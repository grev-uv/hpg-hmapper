/**
 * scheduler.c
 * - Date: 3 / Feb / 2017
 * - Who: Cesar Gonzalez (Cesar.Gonzalez-Segura@uv.es)
 *
 * The scheduler is in charge of creating the producer, consumer and workers,
 * assigning the roles for each and synchronizing them.
 *
 */

#include "scheduler.h"

#include "producer.h"
#include "worker.h"
#include "consumer.h"


//-----------------------------------------------------

scheduler_input_t* scheduler_input_init(const char* index_path, const char* methyl_read_path,
                      const char* hmc_read_path, const char* mc_bam_path, const char* hmc_bam_path,
                      const char* output_dir, size_t num_threads, size_t treatment,
                      size_t memory_budget, size_t batch_size, size_t output_type,
                      char csv_delimiter, char csv_record_delimiter) {
  scheduler_input_t* scheduler = calloc(1, sizeof(scheduler_input_t));

  scheduler->num_threads = num_threads;

  if (num_threads < 3) {
    scheduler->num_workers = 1;
  } else {
    scheduler->num_workers = num_threads - 2;
  }

  scheduler->output_directory = strdup(output_dir);
  scheduler->memory_budget = memory_budget;
  scheduler->methyl_reads = NULL;

  scheduler->memory_consumption = 0.0f;
  scheduler->alignment_mean_size = 0.0f;
  scheduler->alignment_length = 0.0f;

  scheduler->batch_size = batch_size;
  scheduler->output_type = output_type;

  scheduler->csv_delimiter = csv_delimiter;
  scheduler->csv_record_delimiter = csv_record_delimiter;

  // If the reference genome index directory is provided, load the
  // chromosome information. If it is not, use the built-in homo
  // sapiens GRCh37.68 genome
  if (strlen(index_path)) {
    scheduler_load_chr_from_bwt_index(index_path, scheduler);
  } else {
    scheduler->num_chromosomes = DEFAULT_GENOME_NUM_CHROMOSOMES;
    scheduler->chr_length = calloc(DEFAULT_GENOME_NUM_CHROMOSOMES, sizeof(uint64_t));
    memcpy(scheduler->chr_length, DEFAULT_GENOME_CHR_LENGTH, DEFAULT_GENOME_NUM_CHROMOSOMES * sizeof(uint64_t));
  }

  #ifdef DEBUG
  printf("Potentially meth bases per chr:\n");
  printf("===============================\n");
  for (int kk = 0; kk < scheduler->num_chromosomes; ++kk) {
    printf("Chr %d: %lu\n", kk+1, scheduler->chr_length[kk]);
  }
  printf("\n");
  #endif

  // If the per-chromosome methylation statistics file is provided,
  // create the scheduling matrices using that information. If it is
  // not, create it using the reference genome
  scheduler->thr_schedule = calloc(scheduler->num_workers, sizeof(int8_t*));

  for (size_t i = 0; i < scheduler->num_workers; ++i) {
    scheduler->thr_schedule[i] = calloc(scheduler->num_chromosomes, sizeof(int8_t));
    memset(scheduler->thr_schedule[i], THR_UNUSED, scheduler->num_chromosomes * sizeof(int8_t));
  }

  scheduler->chr_schedule = malloc(scheduler->num_chromosomes * sizeof(int8_t));
  memset(scheduler->chr_schedule, THR_UNUSED, scheduler->num_chromosomes * sizeof(int8_t));

  // Populate the matrices
  if (strlen(methyl_read_path)) {
    scheduler_gen_matrices_from_metyhl(methyl_read_path, scheduler);
  } else {
    scheduler_gen_matrices(scheduler->chr_length, scheduler);
  }

  // Load the hydroximethylation statistics
  if (strlen(hmc_read_path)) {
    FILE* hmc_fd = fopen(hmc_read_path, "rb");

    scheduler->hmc_methyl_reads = calloc(scheduler->num_chromosomes, sizeof(uint32_t));
    fread(scheduler->hmc_methyl_reads, sizeof(uint32_t), scheduler->num_chromosomes, hmc_fd);

    fclose(hmc_fd);
  }

  // Open both BAM files
  scheduler->mc_bam_file = bam_fopen_mode((char*)mc_bam_path, NULL, "r");
  scheduler->hmc_bam_file = bam_fopen_mode((char*)hmc_bam_path, NULL, "r");
  scheduler->mc_bam_path = (char*)mc_bam_path;
  scheduler->hmc_bam_path = (char*)hmc_bam_path;

  // Create the worker queues
  scheduler->worker_in_queue = calloc(scheduler->num_workers, sizeof(static_queue_t*));

  for (size_t i = 0; i < scheduler->num_workers; ++i) {
    scheduler->worker_in_queue[i] = static_queue_new(50000);
  }

  scheduler->producer_status = STAGE_IN_PROGRESS;
  scheduler->consumer_status = STAGE_IN_PROGRESS;
  scheduler->worker_team_status = 0;
  scheduler->worker_status = malloc(scheduler->num_workers * sizeof(size_t));
  memset(scheduler->worker_status, STAGE_IN_PROGRESS, scheduler->num_workers * sizeof(size_t));
  scheduler->worker_input = calloc(scheduler->num_workers, sizeof(worker_input_t*));

  return scheduler;
}

//-----------------------------------------------------

void scheduler_alloc_event(scheduler_input_t* scheduler, int64_t size) {
  #pragma omp critical
  {
    scheduler->memory_consumption += (int64_t)(scheduler->alignment_mean_size) * size;
  }
}

//-----------------------------------------------------

void scheduler_alloc_event_direct(scheduler_input_t* scheduler, int64_t size) {
  #pragma omp critical
  scheduler->memory_consumption += size;
}

//-----------------------------------------------------

void scheduler_free_event(scheduler_input_t* scheduler, int64_t size) {
  #pragma omp critical
  {
    int64_t bytes = ((int64_t)scheduler->alignment_mean_size) * size;

    if (scheduler->memory_consumption - bytes >= 0) {
      scheduler->memory_consumption -= bytes;
    } else {
      scheduler->memory_consumption = 0;
    }
  }
}

//-----------------------------------------------------

void scheduler_free_event_direct(scheduler_input_t* scheduler, int64_t size) {
  #pragma omp critical
  scheduler->memory_consumption -= size;
}

//-----------------------------------------------------

int64_t scheduler_get_memory_usage(scheduler_input_t* scheduler) {
  int64_t mem = 0;

  #pragma omp atomic read
  mem = scheduler->memory_consumption;
  
  return mem;
}

//-----------------------------------------------------

void scheduler_set_producer_status(scheduler_input_t* scheduler, size_t status) {
  #pragma omp critical(producer_state_lock)
  scheduler->producer_status = status;
}

//-----------------------------------------------------

size_t scheduler_set_worker_status(scheduler_input_t* scheduler, size_t worker_id, size_t status) {
  #pragma omp critical(worker_state_lock)
  {
    scheduler->worker_status[worker_id] = status;

    if (status == STAGE_COMPLETED) {
      ++scheduler->worker_team_status;
    }

    if (scheduler->worker_team_status == scheduler->num_workers) {
      scheduler->consumer_input->worker_team_pass = STAGE_COMPLETED;
    }
  }
}

//-----------------------------------------------------

void scheduler_finished_worker_stage_event(scheduler_input_t* scheduler) {
  #pragma omp critical(worker_state_lock)
  {
    scheduler->worker_team_pass_status++;

    if (scheduler->worker_team_pass_status < scheduler->num_workers) {
      scheduler->producer_input->worker_team_status = WORKER_FIRST_PASS;
      scheduler->consumer_input->worker_team_pass = WORKER_FIRST_PASS;
    } else if (scheduler->worker_team_pass_status >= scheduler->num_workers &&
               scheduler->worker_team_pass_status < 2 * scheduler->num_workers) {
      scheduler->producer_input->worker_team_status = WORKER_SECOND_PASS;
      scheduler->consumer_input->worker_team_pass = WORKER_SECOND_PASS;
    } else if (scheduler->worker_team_pass_status >= 2 * scheduler->num_workers) {
      scheduler->producer_input->worker_team_status = WORKER_FINISH_PASS;
      scheduler->consumer_input->worker_team_pass = WORKER_FINISH_PASS;
    }
  }
}

//-----------------------------------------------------

void scheduler_set_producer_pass_status(scheduler_input_t* scheduler, size_t pass_id) {
  #pragma omp critical(worker_state_lock)
  {
    scheduler->producer_pass_status = pass_id;

    if (pass_id == PRODUCER_SECOND_PASS) {
      for (size_t i = 0; i < scheduler->num_workers; ++i) {
        worker_pass_finished_event(scheduler->worker_input[i], WORKER_SECOND_PASS);
      }
    } else if (pass_id == PRODUCER_FINISH_PASS) {
      for (size_t i = 0; i < scheduler->num_workers; ++i) {
        worker_pass_finished_event(scheduler->worker_input[i], WORKER_FINISH_PASS);
      }
    }
  }
}

//-----------------------------------------------------

void scheduler_load_chr_from_bwt_index(const char* index_path, scheduler_input_t* scheduler) {
  char fn[MAX_FILENAME_LENGTH];
  char genome_str[128];
  uint32_t* curr_chr;

  snprintf(fn, MAX_FILENAME_LENGTH, "%s/index", index_path);
  FILE* fd = fopen(fn, "r");

  if (fd) {
    // Count lines in the index file. Each line represents
    // a new chromosome
    do {
      if (fgetc(fd) == '\n') scheduler->num_chromosomes++;
    } while (!feof(fd));

    // Rewind the file to the beggining and load the chromosome
    // lengths
    rewind(fd);

    // Chromosome length has the following format:
    // >[name] [start] [finish]
    scheduler->chr_length = calloc(scheduler->num_chromosomes, sizeof(uint64_t));

    for (int i = 0; i < scheduler->num_chromosomes; ++i) {
      fgets(genome_str, 128, fd);
      curr_chr = &scheduler->chr_length[i];
      
      // Get the chromosome length if the line is correctly
      // formatted
      if (strchr(genome_str, '>')) {
        sscanf(genome_str, "%*s %*u %u", curr_chr);
      }
    }

    fclose(fd);
  }
}

//-----------------------------------------------------

// Lambda function to sort the chromosome lengths
int __scheduler_sort_chr(const void* a, const void* b) {
  const uint32_t* x = a;
  const uint32_t* y = b;
  return *y - *x;
}

//-----------------------------------------------------

void scheduler_gen_matrices_from_metyhl(const char* methyl_read_path, scheduler_input_t* scheduler) {
  uint8_t file_num_chromosomes = 0;
  uint32_t* sz_buf = NULL;
  FILE* fd = fopen(methyl_read_path, "rb");

  if (fd) {
    // Check the number of chromosomes in the methyl stats file to
    // ensure that it is compatible with the current genome
    fread(&file_num_chromosomes, sizeof(uint8_t), 1, fd);

    if (file_num_chromosomes != scheduler->num_chromosomes) {
      printf("Error! The input methylation statistics file is not compatible with \
              the current genome or is badly formed. Reverting to the default \
              genome scheduling.\n");
      
      scheduler_gen_matrices(scheduler->chr_length, scheduler);
    } else {
      // Load the number of reads with methylated C's from the file
      scheduler->methyl_reads = calloc(scheduler->num_chromosomes, sizeof(uint32_t));
      fread(scheduler->methyl_reads, sizeof(uint32_t), scheduler->num_chromosomes, fd);

      // Populate the scheduling matrices
      scheduler_gen_matrices(scheduler->methyl_reads, scheduler);
    }

    fclose(fd);
  }
}

//-----------------------------------------------------

void scheduler_gen_matrices(uint64_t* sz_buffer, scheduler_input_t* scheduler) {
  uint64_t* sz_temp_buffer = NULL;
  uint64_t* sz_orig_buffer = NULL;
  size_t* index_buffer = NULL;
  size_t thr_column = 0;
  size_t thread_id = 0;
  int8_t direction = 1;
  
  sz_temp_buffer = calloc(scheduler->num_chromosomes, sizeof(uint64_t));
  index_buffer = calloc(scheduler->num_chromosomes, sizeof(size_t));

  memcpy(sz_temp_buffer, sz_buffer, sizeof(uint64_t) * scheduler->num_chromosomes);

  // Sort the reads in descending order
  qsort(sz_temp_buffer, scheduler->num_chromosomes, sizeof(uint64_t), __scheduler_sort_chr);

  // Reconstruct the indices after the sorting
  for (size_t i = 0; i < scheduler->num_chromosomes; ++i) {
    for (size_t j = 0; j < scheduler->num_chromosomes; ++j) {
      if (sz_temp_buffer[i] == sz_buffer[j]) {
        index_buffer[i] = j;
        break;
      }
    }
  }

  // Fill the scheduler matrices using the top-down ordering
  size_t ib = 0;

  for (size_t i = 0; i < scheduler->num_chromosomes; ++i) {
    scheduler->thr_schedule[thread_id][thr_column] = index_buffer[i];
    scheduler->chr_schedule[index_buffer[i]] = thread_id;

    if (thread_id == scheduler->num_workers - 1 && direction == 1) {
      direction = -1;
      ++thr_column;
    } else if (thread_id == 0 && direction == -1) {
      direction = 1;
      ++thr_column;
    } else {
      thread_id += direction;
    }
  }

  #ifdef DEBUG
  printf("Thread matrix\n");
  printf("==============\n");
  for (size_t i = 0; i < scheduler->num_workers; ++i) {
    for (size_t j = 0; j < scheduler->num_chromosomes; ++j) {
      printf("%d ", scheduler->thr_schedule[i][j]);
    }

    printf("\n");
  }

  printf("\n\n");

  printf("Chromosome matrix\n");
  printf("==================\n");

  for (size_t i = 0; i < scheduler->num_chromosomes; ++i) {
    printf("%d ", scheduler->chr_schedule[i]);
  }
  printf("\n\n");
  #endif

  free(sz_temp_buffer);
  free(index_buffer);
}

//-----------------------------------------------------

void scheduler_print_statistics(scheduler_input_t* scheduler) {
  producer_input_t* producer = scheduler->producer_input;
  consumer_input_t* consumer = scheduler->consumer_input;

  size_t total_reads_mc = 0, processed_reads_mc = 0;
  size_t total_reads_hmc = 0, processed_reads_hmc = 0;

  for (size_t i = 0; i < scheduler->num_chromosomes; ++i) {
    total_reads_mc += producer->total_reads[MC_QUEUE_INDEX][i];
    total_reads_hmc += producer->total_reads[HMC_QUEUE_INDEX][i];

    processed_reads_mc += producer->num_reads[MC_QUEUE_INDEX][i];
    processed_reads_hmc += producer->num_reads[HMC_QUEUE_INDEX][i];
  }
  
  if (scheduler->output_type == STATS_OUTPUT_FORMAT_TEXT) {
    printf("\n");
    printf("===========================================================\n");
    printf("=               S t a t i s t i c s                       =\n");
    printf("===========================================================\n");
    printf("\n");

    printf("Total processed reads (mC): %lu\n", total_reads_mc);
    printf("Total processed reads (hmC): %lu\n\n", total_reads_hmc);

    printf("Skipped reads (mC): %lu (%.2f%%)\n", total_reads_mc - processed_reads_mc,
          100.0*((double)(total_reads_mc - processed_reads_mc) / (double)(total_reads_mc)));

    printf("Skipped reads (hmC): %lu (%.2f%%)\n\n", total_reads_hmc - processed_reads_hmc,
          100.0*((double)(total_reads_hmc - processed_reads_hmc) / (double)(total_reads_hmc)));
    
    printf("Total processed cytosines: %lu\n", consumer->total_c);
    printf("Total methylated cytosines: %lu (%.2f%%)\n", consumer->total_mc,
          100.0*((double)(consumer->total_mc)/(double)consumer->total_c));

    printf("Total hydroximethylated cytosines: %lu (%.2f%%)\n\n", consumer->total_hmc,
          100.0*((double)(consumer->total_hmc)/(double)consumer->total_c));
  } else if (scheduler->output_type == STATS_OUTPUT_FORMAT_CSV) {
    // First line is the CSV header, second line is the data
    printf("Total 5mC Reads\tTotal 5hmC Reads\tUnprocessed 5mC Reads\tUnprocessed 5hmC Reads\n");
    printf("%lu\t%lu\t", total_reads_mc, total_reads_hmc);
    printf("%lu\t%lu\t\n", total_reads_mc - processed_reads_mc, total_reads_hmc - processed_reads_hmc);
  }

  // Print timing data
  printf("\n\n");

  double producer_time = 0.0, worker_time = 0.0, consumer_time = 0.0, total_time = 0.0;
  double producer_bam_time = 0.0, producer_filter_time = 0.0;
  double worker_first_stage_process = 0.0, worker_bit_to_map = 0.0, 
         worker_second_stage_process = 0.0;
  
  producer_time = scheduler->producer_input->time_sec;
  producer_bam_time = scheduler->producer_input->time_bam_read;
  producer_filter_time = scheduler->producer_input->time_filter_alignment;

  consumer_time = scheduler->consumer_input->time_sec;

  for (size_t i = 0; i < scheduler->num_workers; ++i) {
    worker_time += scheduler->worker_input[i]->time_sec;
    worker_first_stage_process += scheduler->worker_input[i]->first_stage_process;
    worker_second_stage_process += scheduler->worker_input[i]->second_stage_process;
    worker_bit_to_map += scheduler->worker_input[i]->bit_to_map_time;
  }

  worker_time /= (double)scheduler->num_workers;
  worker_first_stage_process /= (double)scheduler->num_workers;
  worker_second_stage_process /= (double)scheduler->num_workers;
  worker_bit_to_map /= (double)scheduler->num_workers;

  worker_time -= (worker_first_stage_process + worker_second_stage_process);

  total_time = fmax(producer_time + producer_bam_time + producer_filter_time, 
                  worker_time + worker_first_stage_process + worker_second_stage_process + 
                   + worker_bit_to_map) + consumer_time;

  if (scheduler->output_type == STATS_OUTPUT_FORMAT_TEXT) {
    printf("\n");
    printf("===========================================================\n");
    printf("=                    T i m i n g                          =\n");
    printf("===========================================================\n");

    printf("[Producer (other)]\t%.6f sec (%.2f%% of total)\n",
              producer_time, 100.0 * (producer_time / total_time));
    printf("[Producer (read BAM)]\t%.6f sec (%.2f%% of total)\n",
              producer_bam_time, 100.0 * (producer_bam_time / total_time));
    printf("[Producer (filter read)]\t%.6f sec (%.2f%% of total)\n",
              producer_filter_time, 100.0 * (producer_filter_time / total_time));

    printf("[Workers (other)]\tElapsed time: %.6f sec (%.2f%% of total)\n",
              worker_time, 100.0 * (worker_time / total_time));
    printf("[Workers (1st stage)]\tElapsed time: %.6f sec (%.2f%% of total)\n",
              worker_first_stage_process, 100.0 * (worker_first_stage_process / total_time));
    printf("[Workers (Bit-array to map)]\tElapsed time: %.6f sec (%.2f%% of total)\n",
              worker_bit_to_map, 100.0 * (worker_bit_to_map / total_time));
    printf("[Workers (2nd stage)]\tElapsed time: %.6f sec (%.2f%% of total)\n",
              worker_second_stage_process, 100.0 * (worker_second_stage_process / total_time));

    printf("[Consumer]\tElapsed time: %.6f sec (%.2f%% of total)\n",
              consumer_time, 100.0 * (consumer_time / total_time));
    
    printf("[TOTAL]\tElapsed time: %.6f sec\n", total_time);
  } else if (scheduler->output_type == STATS_OUTPUT_FORMAT_CSV) {
    // First line is the CSV header, second line is the timing data
    printf("Total Producer Time\tBAM Read Time\tAlignment Filtering Time\tTotal Worker Time\t\
            Worker First Stage Time\tWorker Bit to Map Time\tWorker Second Stage Time\t\
            Consumer Time\tTotal Time\n");
    printf("%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
      producer_time, producer_bam_time, producer_filter_time, worker_time,
      worker_first_stage_process, worker_bit_to_map, worker_second_stage_process,
      consumer_time, total_time);
  }
}

//-----------------------------------------------------

void scheduler_input_free(scheduler_input_t* scheduler) {
  if (scheduler) {
    // Destroy the producer
    if (scheduler->producer_input) {
      producer_input_free(scheduler->producer_input);
    }

    // Destroy the consumer
    if (scheduler->consumer_input) {
      consumer_input_free(scheduler->consumer_input);
    }

    // Destroy the workers
    if (scheduler->worker_input) {
      for (size_t i = 0; i < scheduler->num_workers; ++i) {
        if (scheduler->worker_input[i]) {
          worker_input_free(scheduler->worker_input[i]);
        }
      }

      free(scheduler->worker_input);
    }

    if (scheduler->chr_length) {
      free(scheduler->chr_length);
    }

    if (scheduler->thr_schedule) {
      for (size_t i = 0; i < scheduler->num_workers; ++i) {
        if (scheduler->thr_schedule[i]) {
          free(scheduler->thr_schedule[i]);
        }
      }

      free (scheduler->thr_schedule);
    }

    if (scheduler->chr_schedule) {
      free(scheduler->chr_schedule);
    }

    if (scheduler->methyl_reads) {
      free(scheduler->methyl_reads);
    }

    if (scheduler->hmc_methyl_reads) {
      free(scheduler->hmc_methyl_reads);
    }

    bam_fclose(scheduler->mc_bam_file);
    bam_fclose(scheduler->hmc_bam_file);

    if (scheduler->worker_in_queue) {
      for (size_t i = 0; i < scheduler->num_workers; ++i) {
        static_queue_free(scheduler->worker_in_queue[i]);
      }

      free(scheduler->worker_in_queue);
    }

    if (scheduler->output_directory) {
      free(scheduler->output_directory);
    }

    if (scheduler->worker_status) {
      free(scheduler->worker_status);
    }
    
    free(scheduler);
  }
}

//-----------------------------------------------------

void scheduler_run(scheduler_input_t* scheduler, scheduler_producer_cb_t producer_cb, scheduler_consumer_cb_t consumer_cb,
                        scheduler_worker_cb_t worker_cb) {
  // Set the consumer, producer and workers callback pointers
  scheduler->producer_cb = producer_cb;
  scheduler->consumer_cb = consumer_cb;
  scheduler->worker_cb = worker_cb;

  size_t num_threads = scheduler->num_threads;
  size_t thread_id = 0;
  
  // Launch the pipeline
  if (num_threads < 3) {
    scheduler_run_single_threaded(scheduler);
  } else {
    omp_set_num_threads(num_threads);

    #pragma omp parallel shared(scheduler) private(thread_id, num_threads)
    {
      thread_id = omp_get_thread_num();
    
      if (thread_id == 0) {
        scheduler_producer_start(scheduler);
      } else if (thread_id == 1) {
        scheduler_consumer_start(scheduler);
      } else {
        scheduler_worker_start(scheduler, thread_id);
      }
    }
  }

  scheduler_print_statistics(scheduler);
}

//-----------------------------------------------------

void scheduler_run_single_threaded(scheduler_input_t* scheduler) {
  scheduler->producer_input = producer_input_init(scheduler, 0);
  scheduler->consumer_input = consumer_input_init(scheduler, 0);
  scheduler->worker_input[0] = worker_input_init(scheduler, 0);

  size_t producer_stage = STAGE_IN_PROGRESS;
  size_t worker_stage = STAGE_IN_PROGRESS;
  size_t consumer_stage = STAGE_IN_PROGRESS;

  scheduler_producer_cb_t produce = scheduler->producer_cb;
  scheduler_consumer_cb_t consume = scheduler->consumer_cb;
  scheduler_worker_cb_t work = scheduler->worker_cb;

  while (consumer_stage != STAGE_COMPLETED) {
    // Producer
    if (producer_stage == STAGE_IN_PROGRESS) {
      producer_stage = produce(scheduler->producer_input, scheduler);

      if (producer_stage == STAGE_COMPLETED) {
        scheduler_set_producer_status(scheduler, STAGE_COMPLETED);
      }
    }

    // Worker
    if (worker_stage == STAGE_IN_PROGRESS) {
      worker_stage = work(scheduler->worker_input[0], scheduler);
    }

    // Consumer
    if (producer_stage == STAGE_COMPLETED && worker_stage == STAGE_COMPLETED) {
      consumer_stage = consume(scheduler->consumer_input, scheduler);
    }
  }
}

//-----------------------------------------------------

void scheduler_producer_start(scheduler_input_t* scheduler) {
  scheduler->producer_input = producer_input_init(scheduler, omp_get_thread_num());
  producer_input_t* producer = scheduler->producer_input;

  // Local pointer to the producer callback
  scheduler_producer_cb_t produce = scheduler->producer_cb;

  // Run the producer until it finishes
  while (produce(producer, scheduler) != STAGE_COMPLETED);

  // Mark that the producer has finished
  scheduler_set_producer_status(scheduler, STAGE_COMPLETED);
}

//-----------------------------------------------------

void scheduler_consumer_start(scheduler_input_t* scheduler) {
  scheduler->consumer_input = consumer_input_init(scheduler, omp_get_thread_num());
  consumer_input_t* consumer = scheduler->consumer_input;

  // Local pointer to the consumer callback
  scheduler_consumer_cb_t consume = scheduler->consumer_cb;

  // Run the consumer until it finishes
  while (consume(consumer, scheduler) != STAGE_COMPLETED);
}

//-----------------------------------------------------

void scheduler_worker_start(scheduler_input_t* scheduler, size_t id) {
  scheduler->worker_input[id - 2] = worker_input_init(scheduler, id);
  worker_input_t* worker = scheduler->worker_input[id - 2];

  // Local pointer to the worker callback
  scheduler_worker_cb_t work = scheduler->worker_cb;

  // Run the worker until it finishes
  while (work(worker, scheduler) != STAGE_COMPLETED);

  // Mark that the current worker has finished
  scheduler_set_worker_status(scheduler, worker->worker_id, STAGE_COMPLETED);
}
