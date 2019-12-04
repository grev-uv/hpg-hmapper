#include "consumer.h"


//-----------------------------------------------------

consumer_input_t* consumer_input_init(scheduler_input_t* scheduler, size_t thread_id) {
  char map_path[MAX_FILENAME_LENGTH] = { 0 };
  char map_name[128] = { 0 };

  consumer_input_t* ptr     = calloc(1, sizeof(consumer_input_t));

  ptr->num_chromosomes      = scheduler->num_chromosomes;
  ptr->num_workers          = scheduler->num_workers;
  ptr->thread_id            = thread_id;
  ptr->worker_team_pass     = WORKER_FIRST_PASS;

  ptr->csv_delimiter        = scheduler->csv_delimiter;
  ptr->csv_record_delimiter = scheduler->csv_record_delimiter;
  ptr->coverage             = scheduler->coverage;

  ptr->meth_map_forward_fd  = calloc(ptr->num_chromosomes, sizeof(FILE*));
  ptr->meth_map_reverse_fd  = calloc(ptr->num_chromosomes, sizeof(FILE*));

  for (size_t i = 0; i < ptr->num_chromosomes; ++i) {
    strcpy(map_path, scheduler->output_directory);
    snprintf(map_name, 128, "methylation_map_forward_%lu.csv", i + 1);
    strcat(map_path, map_name);

    ptr->meth_map_forward_fd[i] = fopen(map_path, "w");

    strcpy(map_path, scheduler->output_directory);
    snprintf(map_name, 128, "methylation_map_reverse_%lu.csv", i + 1);
    strcat(map_path, map_name);

    ptr->meth_map_reverse_fd[i] = fopen(map_path, "w");
  }

  ptr->time_sec = 0.0;

  return ptr;
}

//-----------------------------------------------------

void consumer_input_free(consumer_input_t* consumer_input) {
  if (consumer_input) {
    if (consumer_input->meth_map_forward_fd) {
      for (size_t i = 0; i < consumer_input->num_chromosomes; ++i) {
        fclose(consumer_input->meth_map_forward_fd[i]);
      }

      free(consumer_input->meth_map_forward_fd);
    }

    if (consumer_input->meth_map_reverse_fd) {
      for (size_t i = 0; i < consumer_input->num_chromosomes; ++i) {
        fclose(consumer_input->meth_map_reverse_fd[i]);
      }

      free(consumer_input->meth_map_reverse_fd);
    }
  
    free(consumer_input);
  }
}

//-----------------------------------------------------

int consumer_stage_step(void* data, scheduler_input_t* scheduler) {
  int result = STAGE_IN_PROGRESS;

  consumer_input_t* input = (consumer_input_t*)data;

  meth_array_node_t** current_array = NULL;
  size_t* current_array_size = NULL;
  int8_t* thr_schedule = NULL;
  int8_t chromosome = 0;

  // Timing data
  double start_time;

  // Only serialize the maps if all the workers have finished
  if (input->worker_team_pass == WORKER_FINISH_PASS) {
    start_time = omp_get_wtime();

    // Traverse all the existing worker arrays and serialize them
    for (size_t worker = 0; worker < input->num_workers; ++worker) {
      thr_schedule = scheduler->worker_input[worker]->thr_schedule;

      current_array = scheduler->worker_input[worker]->meth_array_forward;
      current_array_size = scheduler->worker_input[worker]->meth_array_forward_size;

      for (size_t chr = 0; chr < input->num_chromosomes && thr_schedule[chr] != -1; ++chr) {
        consumer_meth_array_serialize(current_array[chr], 
                                      current_array_size[chr], 
                                      input->meth_map_forward_fd[thr_schedule[chr]],
                                      input->csv_delimiter, 
                                      input->csv_record_delimiter,
                                      input->coverage);
      }

      current_array = scheduler->worker_input[worker]->meth_array_reverse;
      current_array_size = scheduler->worker_input[worker]->meth_array_reverse_size;

      for (size_t chr = 0; chr < input->num_chromosomes && thr_schedule[chr] != -1; ++chr) {
        consumer_meth_array_serialize(current_array[chr], 
                                      current_array_size[chr], 
                                      input->meth_map_reverse_fd[thr_schedule[chr]],
                                      input->csv_delimiter, 
                                      input->csv_record_delimiter,
                                      input->coverage);
      }
    }

    input->time_sec += omp_get_wtime() - start_time;

    #ifdef DEBUG
    printf("[Cons] Stage completed.\n");
    #endif

    result = STAGE_COMPLETED;
  }

  return result;
}

//-----------------------------------------------------

void consumer_meth_array_serialize(meth_array_node_t* array, size_t length, 
        FILE* fd, char delimiter, char record_delimiter, size_t coverage) {
  uint32_t position = 0;
  uint16_t c = 0, nc = 0, mc = 0, hmc = 0;

  for (size_t i = 0; i < length; ++i) {
    position = array[i].position;
    c = array[i].c_count;
    nc = array[i].nc_count;
    mc = array[i].mc_count;
    hmc = array[i].hmc_count;

    if (c + mc > coverage || c + hmc > coverage) {
      fprintf(fd, "%u%c%u%c%u%c%u%c%u%c", 
        position, delimiter, 
        c, delimiter, 
        nc, delimiter, 
        mc, delimiter, 
        hmc, record_delimiter);
    }
  }
}
