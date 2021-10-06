/**
 * worker.c
 * - Date: 15 / Feb / 2017
 * - Who: Cesar Gonzalez (Cesar.Gonzalez-Segura@uv.es)
 *
 * Workers are in charge of retrieving alignments from the input queue, processing them
 * to create the methylation and hydroximethylation map and storing the result into the
 * output queue.
 *
 */

#include "worker.h"

static const uint8_t bit_map_mask[] = {
  0b00000001,
  0b00000010,
  0b00000100,
  0b00001000,
  0b00010000,
  0b00100000,
  0b01000000,
  0b10000000
};

//-----------------------------------------------------

worker_input_t* worker_input_init(scheduler_input_t* scheduler, size_t thread_id) {
  worker_input_t* worker = calloc(1, sizeof(worker_input_t));

  worker->thread_id = thread_id;

  if (scheduler->num_threads < 3) {
    worker->worker_id = 0;
  } else {
    worker->worker_id = thread_id - 2;
  }

  worker->num_chromosomes = scheduler->num_chromosomes;
  
  worker->time_sec = 0.0;
  worker->temp_reads = 0;
  worker->processed_reads = 0;

  worker->batch_size = scheduler->batch_size;
  worker->in_queue = scheduler->worker_in_queue[worker->worker_id];

  worker->thr_schedule = calloc(worker->num_chromosomes, sizeof(int8_t));
  memcpy(worker->thr_schedule, scheduler->thr_schedule[worker->worker_id], worker->num_chromosomes * sizeof(int8_t));

  worker->bit_map_forward = calloc(worker->num_chromosomes, sizeof(uint8_t*));
  worker->bit_map_reverse = calloc(worker->num_chromosomes, sizeof(uint8_t*));

  worker->meth_array_forward = calloc(worker->num_chromosomes, sizeof(meth_array_node_t*));
  worker->meth_array_reverse = calloc(worker->num_chromosomes, sizeof(meth_array_node_t*));

  worker->meth_array_forward_size = calloc(worker->num_chromosomes, sizeof(size_t));
  worker->meth_array_reverse_size = calloc(worker->num_chromosomes, sizeof(size_t));
  worker->bit_map_size = calloc(worker->num_chromosomes, sizeof(size_t));

  for (size_t i = 0; i < worker->num_chromosomes && worker->thr_schedule[i] != -1; ++i) {
    worker->bit_map_size[i] = 1 + (DEFAULT_GENOME_CHR_LENGTH[worker->thr_schedule[i]] / 8);
    worker->bit_map_forward[i] = calloc(worker->bit_map_size[i], sizeof(uint8_t));
    worker->bit_map_reverse[i] = calloc(worker->bit_map_size[i], sizeof(uint8_t));
  }

  worker->pass_id = WORKER_FIRST_PASS;
  worker->next_pass = WORKER_FIRST_PASS;

  return worker;
}

//-----------------------------------------------------

void worker_input_free(worker_input_t* worker_input) {
  if (worker_input) {
    if (worker_input->meth_array_forward) {
      for (size_t i = 0; i < worker_input->num_chromosomes && worker_input->thr_schedule[i] != -1; ++i) {
        free(worker_input->meth_array_forward[i]);
      }

      free(worker_input->meth_array_forward);
    }

    if (worker_input->meth_array_reverse) {
      for (size_t i = 0; i < worker_input->num_chromosomes && worker_input->thr_schedule[i] != -1; ++i) {
        free(worker_input->meth_array_reverse[i]);
      }

      free(worker_input->meth_array_reverse);
    }

    if (worker_input->thr_schedule) {
      free(worker_input->thr_schedule);
    }
	
    if (worker_input->bit_map_size) {
      free(worker_input->bit_map_size);
    }

    free(worker_input);
  }
}

//-----------------------------------------------------

int worker_stage_step(void* data, scheduler_input_t* scheduler) {
  int status = STAGE_IN_PROGRESS;
  worker_input_t* worker_input = (worker_input_t*)data;

  // Input queues
  alignment_batch_t* current_batch = NULL;
  alignment_t* current_alignment = NULL;

  size_t input_count = 0;
  size_t batch_size = 0;

  // Timing data
  double start_time = 0.0;

  // Check if there is data in the input queue and if the producer is still in
  // progress. If the input queue is empty and the producer has finished,
  // mark the worker as finished.
  current_batch = static_queue_pop(worker_input->in_queue);

  if (current_batch) {
    start_time = omp_get_wtime();

    // Process the alignments, and when the processing is over insert them
    // into the second stage queue to be reprocessed later
    current_alignment = static_queue_pop(current_batch->alignments);

    while (current_alignment) {
      worker_process_alignment(worker_input, current_alignment, current_batch->type, scheduler, worker_input->pass_id);

      alignment_free(current_alignment);

      worker_input->temp_reads++;
      worker_input->processed_reads++;

      if (worker_input->temp_reads == SCHEDULER_MEMORY_EVENT_READ_SIZE) {
        worker_input->temp_reads = 0;
        scheduler_free_event(scheduler, SCHEDULER_MEMORY_EVENT_READ_SIZE);
      }
      
      current_alignment = static_queue_pop(current_batch->alignments);
    }

    alignment_batch_free(current_batch);

    worker_input->time_sec += omp_get_wtime() - start_time;
  } else {
    if (worker_input->pass_id == WORKER_FIRST_PASS && worker_input->next_pass == WORKER_SECOND_PASS) {
      #ifdef DEBUG
      printf("[W=%lu] First pass completed.\n", worker_input->worker_id);
      #endif

      worker_input->temp_reads = 0;

      // Convert the methylation trees to linear arrays
      worker_meth_map_to_array(worker_input);

      scheduler_finished_worker_stage_event(scheduler);
      static_queue_clear(worker_input->in_queue);

      worker_input->pass_id = WORKER_SECOND_PASS;
    } else if (worker_input->pass_id == WORKER_SECOND_PASS && worker_input->next_pass == WORKER_FINISH_PASS) {
      #ifdef DEBUG
      printf("[W=%lu] Second pass completed.\n", worker_input->worker_id);
      #endif

      scheduler_finished_worker_stage_event(scheduler);

      status = STAGE_COMPLETED;
    }
  }

  return status;
}

//-----------------------------------------------------

void worker_process_alignment(worker_input_t* worker, alignment_t* alignment, size_t type, 
          scheduler_input_t* scheduler, size_t pass_id) {
  size_t length = alignment->length;
  size_t strand = alignment->seq_strand;

  bam_tag_t* current_tag = NULL;
  bam_string_t* bam_str = NULL;

  char* sequence = NULL, *meth_sequence = NULL;
  double start_time = 0.0, start_time_sub = 0.0;

  if (pass_id == WORKER_FIRST_PASS) {
    uint8_t* current_map = NULL;
    size_t array_size = NULL;

    start_time = omp_get_wtime();

    // Extract the methylation context from the optional tags
    for (size_t i = array_list_size(alignment->optional_tags) - 1; i >= 0; --i) {
      current_tag = array_list_get(i, alignment->optional_tags);

      if (current_tag != NULL && !strcmp(current_tag->tag, XM_TAG_NAME)) {
        bam_str = current_tag->data;
        sequence = bam_str->data;
        break;
      }
    }

    // Select the methylation map for the current chromosome
    for (size_t i = 0; i < worker->num_chromosomes && worker->thr_schedule[i] != -1; ++i) {
      if (worker->thr_schedule[i] == alignment->chromosome) {
        if (strand == STRAND_FORWARD) {
          current_map = worker->bit_map_forward[i];
        } else {
          current_map = worker->bit_map_reverse[i];
        }

        array_size = worker->bit_map_size[i];
        break;
      }
    }

    // Traverse the methylation context and search the methylated C's
    size_t compressed_position = alignment->position / 8;
	  size_t mask = alignment->position % 8;

	  //RICARDO - Incorporando el cigar
	  int cont, pos=0, operations, num;
	  int offset = 0, current_offset=0;
	  char car;
	  int inser = 0;
	  int pos_mod = 0;
	  int pp = 0;
	  char *cigar = strdup(alignment->cigar);
	  for (operations = 0; operations < alignment->num_cigar_operations; operations++) {
	      sscanf(&cigar[offset], "%i%c%n", &num, &car, &current_offset);
	      offset += current_offset;

	      if (car == 'M' || car == '=' || car == 'X') {
	        for (cont = 0; cont < num; cont++, pos++, pos_mod++) {

	        	//alignment->sequence_mod[pos_mod] = alignment->sequence[pos];
	        	//alignment->sequence_met[pos_mod] = sequence[pos];


	        	 if (sequence[pos] != XM_NON_RELEVANT) {
	        	        #ifdef DEBUG
	        	        if (compressed_position >= array_size) {
	        	          printf("Warning! Pc (%lu) < N (%lu)\n", compressed_position, array_size);
	        	          raise(SIGINT);
	        	        }
	        	        #endif
	        	        if (pos_mod<alignment->length){
	        	        	current_map[compressed_position] |= bit_map_mask[mask];
	        	        	}
	        	        /*else{
	        	        	pp = 5;
	        	        }*/

	        		    if (mask == 7) {
	        		      ++compressed_position;
	        			    mask = 0;
	        		    } else {
	        			    ++mask;
	        		    }
	        	 }
	        }

	   }
	    else{
	    	   if (car == 'D' || car == 'N') {
	    		   for (cont = 0; cont < num; cont++, pos_mod++) {

	    			  // alignment->sequence_mod[pos_mod] = 'D';
	    			  // alignment->sequence_met[pos_mod] = '.';

	    			   	   	   if (mask == 7) {
	    			  	             ++compressed_position;
	    			  	       	    mask = 0;
	    			  	       } else {
	    			  	       	    ++mask;
	    			  	       }
	    	   }

	      }
	      else {
	    	       if (car == 'I' || car == 'H' || car == 'S') {
	    	           pos+=num;

	    	       }
	    	   }


	  }
	  }

	  free(cigar);


	 /* for (size_t i = 0; i < length; ++i) {
      if (sequence[i] != XM_NON_RELEVANT) {
        #ifdef DEBUG
        if (compressed_position >= array_size) {
          printf("Warning! Pc (%lu) < N (%lu)\n", compressed_position, array_size);
          raise(SIGINT);
        }
        #endif

        current_map[compressed_position] |= bit_map_mask[mask];
      }
	  
	    if (mask == 7) {
	      ++compressed_position;
		    mask = 0;
	    } else {
		    ++mask;
	    }
    }*/

    worker->first_stage_process += omp_get_wtime() - start_time;
  } else if (pass_id == WORKER_SECOND_PASS) {
    meth_array_node_t* current_array = NULL;
    meth_array_node_t* current_node = NULL;
    size_t array_size = 0;
    size_t current_position = 0;
    int8_t inside_interval = 1;

	  start_time = omp_get_wtime();
	
    // Select the methylation array for the current chromosome
    for (size_t i = 0; i < worker->num_chromosomes && worker->thr_schedule[i] != -1; ++i) {
      if (worker->thr_schedule[i] == alignment->chromosome) {
        if (strand == STRAND_FORWARD) {
          current_array = worker->meth_array_forward[i];
          array_size = worker->meth_array_forward_size[i];
        } else {
          current_array = worker->meth_array_reverse[i];
          array_size = worker->meth_array_reverse_size[i];
        }
      }
    }

    // Extract the methylation context from the optional tags
    for (size_t i = 0; i < array_list_size(alignment->optional_tags); ++i) {
      current_tag = array_list_get(i, alignment->optional_tags);

      if (current_tag != NULL && !strcmp(current_tag->tag, XM_TAG_NAME)) {
        bam_str = current_tag->data;
        meth_sequence = bam_str->data;
        break;
      }
    }

    // Search the left-most closest entry to the interval created by
    // the alignment position and its length
    current_node = worker_meth_array_find(current_array, array_size, alignment->position, &current_position);
    current_position++;	//Por coherencia, avanzamos el siguiente índice
    // If the selected node isn't inside the alignment interval, advance the node
    // until it is inside or until the node gets out of the interval. If it is out
    // of the interval, assume that it doesn't have any relevant hit and skip the
    // alignment.
    if (current_node->position < alignment->position) {
      while (current_node->position < alignment->position && current_position < array_size + 1) {/*+1, no -1*/
        current_node = &current_array[current_position++];
      }

      if (current_node->position >= alignment->position + length) {
        inside_interval = 0;
      }
    }


    //RICARDO - Incorporando el cigar

    char* sequence_mod = (char*) calloc(2*alignment->length + 1, sizeof(char)); //*2 for deleccions
    char* sequence_met = (char*) calloc(2*alignment->length + 1, sizeof(char)); //*2 for deleccions

    	int cont, pos=0, operations, num;
    	  int offset = 0, current_offset=0;
    	  char car;
    	  int inser = 0;
    	  int pos_mod = 0;
    	  char *cigar = strdup(alignment->cigar);
    	  for (operations = 0; operations < alignment->num_cigar_operations; operations++) {
    	      sscanf(&cigar[offset], "%i%c%n", &num, &car, &current_offset);
    	      offset += current_offset;

    	      if (car == 'M' || car == '=' || car == 'X') {
    	        for (cont = 0; cont < num; cont++, pos++, pos_mod++) {

    	        	sequence_mod[pos_mod] = alignment->sequence[pos];
    	        	sequence_met[pos_mod] = meth_sequence[pos];


    	        }

    	   }
    	    else{
    	    	   if (car == 'D' || car == 'N') {
    	    		   for (cont = 0; cont < num; cont++, pos_mod++) {

    	    			   sequence_mod[pos_mod] = 'D';
    	    			   sequence_met[pos_mod] = '.';

    	    	   }

    	      }
    	      else {
    	    	       if (car == 'I' || car == 'H' || car == 'S') {
    	    	           pos+=num;

    	    	       }
    	    	   }


    	  }
    	  }
    	  while(pos_mod < alignment->length)
    	  {
    		  sequence_mod[pos_mod] = 'R';	//Rellenamos secuencia
    		  sequence_met[pos_mod] = '.';
    		  pos_mod++;
    	  }

    	  free(cigar);

    // If there exist hits inside the alignment interval, traverse the hits on the
    // interval and check the sequence for the corresponding nucleotides. If the
    // nucleotides corresponding to a hit are not methylated values (C for + strand or
    // G for - strand) count the base as a non-C.
    char meth_base, seq_base;

    if (inside_interval && meth_sequence) {
      while (current_node->position < alignment->position + length - 1 && current_position < array_size + 1/*- 1*/) {  //OJO, el current_position vale uno de más, originalmente era -1, pero realmente es +1
    /*    seq_base = bam1_seqi(alignment->sequence, current_node->position - alignment->position);
        meth_base = meth_sequence[current_node->position - alignment->position];*/
    	seq_base = sequence_mod[current_node->position - alignment->position];
    	meth_base = sequence_met[current_node->position - alignment->position];

        if (meth_base != XM_NON_RELEVANT) {
          if (meth_base == XM_METHYLATED_CPG ||
              meth_base == XM_METHYLATED_CHG ||
              meth_base == XM_METHYLATED_CHH ||
              meth_base == XM_METHYLATED_MUT ||
			  meth_base == XM_METHYLATED_CUN) {
            if (alignment->methylation_type == MC_QUEUE_INDEX) {
              current_node->mc_count++;
            } else {
              current_node->hmc_count++;
            }
          } else {
            if (alignment->methylation_type == MC_QUEUE_INDEX) {
              current_node->c_count++;
            } else {
              current_node->ch_count++;
            }
          }
        } else {
          if (strand == STRAND_FORWARD) {
            if (seq_base != 'C') {
              if (alignment->methylation_type == MC_QUEUE_INDEX) {
                current_node->nc_count++;
              } else {
                current_node->nch_count++;
              }
            }
          } else {
            if (seq_base != 'G') {
              if (alignment->methylation_type == MC_QUEUE_INDEX) {
                current_node->nc_count++;
              } else {
                current_node->nch_count++;
              }
            }
          }
        }

        current_node = &current_array[current_position++]; //OJO Current position vale 1 más de los datos que tiene el current_node
      }
    }

    free(sequence_mod);
    free(sequence_met);


    worker->second_stage_process += omp_get_wtime() - start_time;
  }
}

//-----------------------------------------------------

void worker_pass_finished_event(worker_input_t* worker, size_t pass_id) {
  #ifdef DEBUG
  printf("[W=%lu] Passing from stage %lu to stage %lu\n", worker->worker_id, worker->pass_id, pass_id);
  #endif
  
  worker->next_pass = pass_id;
}

//-----------------------------------------------------

meth_array_node_t* worker_meth_array_find(meth_array_node_t* array, size_t length, uint64_t position, size_t* current_position) {
  size_t next = length * 0.5;
  size_t start = 0, end = length;
  meth_array_node_t* result = NULL;

  // Handle the edge cases when the methylation maps have length one or two,
  // and the general case when the maps have length greater than two
  if (length == 1) {
    result = &array[0];
    *current_position = 0;
  } else if (length == 2) {
    if (array[0].position <= position) {
      result = &array[0];
      *current_position = 0;
    } else {
      result = &array[1];
      *current_position = 1;
    }
  } else {
    result = &array[next];

    while (result->position != position && (end - start) != 1) {
      if (result->position > position) {
        end = next;
        next = start + ((end - start) * 0.5);
        result = &array[next];
      } else if (result->position < position) {
        start = next;
        next = start + ((end - start) * 0.5);
        result = &array[next];
      }
    }
  }
  
  *current_position = next;
  return result;
}

//-----------------------------------------------------

void worker_meth_map_to_array(worker_input_t* worker) {
  uint8_t* current_map = NULL;
  meth_array_node_t* last_meth_array = NULL;
  size_t last_meth_array_index = 0;

  size_t comp_index = 0, comp_offset = 0, array_index = 0;
  uint8_t current = 0;
  
  double start_time = omp_get_wtime();

  // Convert each methylation map into an array, first the forward strand
  // maps and then the reverse strand maps. After the data has been copied,
  // free the methylation trees.
  for (size_t i = 0; i < worker->num_chromosomes && worker->thr_schedule[i] != -1; ++i) {
	  // Forward map
    current_map = worker->bit_map_forward[i];

    // Traverse the bit array and find how many methylated entries there are
    comp_index = 0;
    comp_offset = 0;
    array_index = 0;
	
    while (comp_index < worker->bit_map_size[i]) {
      current = (current_map[comp_index] >> comp_offset) & 1;

      if (current) {
        worker->meth_array_forward_size[i]++;
      }
      
      if (comp_offset == 7) {
        ++comp_index;
        comp_offset = 0;
      } else {
        ++comp_offset;
      }
    }
    
    // Allocate the methylation map data
    // If there were zero methylated entries in the chromosome, allocate a 1 sized
    // map to avoid problems further in the pipeline
    if (worker->meth_array_forward_size[i] == 0) {
      worker->meth_array_forward_size[i] = 1;
    }

    worker->meth_array_forward[i] = calloc(worker->meth_array_forward_size[i], sizeof(meth_array_node_t));
    
    #ifdef DEBUG
    printf("[Worker %lu] Forward map for chr %d with %lu elements (%lu MB)\n", worker->worker_id, worker->thr_schedule[i],
        worker->meth_array_forward_size[i], worker->meth_array_forward_size[i] * 16 / 1000000);
    #endif

    // Traverse the bit array again and set the initial positions
    comp_index = 0;
    comp_offset = 0;
    array_index = 0;
	
    while (comp_index < worker->bit_map_size[i]) {
      current = (current_map[comp_index] >> comp_offset) & 1;

      if (current) {
        worker->meth_array_forward[i][array_index++].position = comp_offset + (comp_index * 8);
      }
      
      if (comp_offset == 7) {
        ++comp_index;
        comp_offset = 0;
      } else {
        ++comp_offset;
      }
    }

    // Repeat for the reverse strand
    current_map = worker->bit_map_reverse[i];

    // Traverse the bit array and find how many methylated entries there are
    comp_index = 0;
    comp_offset = 0;
    array_index = 0;
	
    while (comp_index < worker->bit_map_size[i]) {
      current = (current_map[comp_index] >> comp_offset) & 1;

      if (current) {
        worker->meth_array_reverse_size[i]++;
      }
      
      if (comp_offset == 7) {
        ++comp_index;
        comp_offset = 0;
      } else {
        ++comp_offset;
      }
    }
    
    // Allocate the methylation map data
    // If there were zero methylated entries in the chromosome, allocate a 1 sized
    // map to avoid problems further in the pipeline
    if (worker->meth_array_reverse_size[i] == 0) {
      worker->meth_array_reverse_size[i] = 1;
    }

    #ifdef DEBUG
    printf("[Worker %lu] Reverse map for chr %d with %lu elements (%lu MB)\n", worker->worker_id, worker->thr_schedule[i],
        worker->meth_array_forward_size[i], worker->meth_array_forward_size[i] * 16 / 1000000);
    #endif

    worker->meth_array_reverse[i] = calloc(worker->meth_array_reverse_size[i], sizeof(meth_array_node_t));
    
    // Traverse the bit array again and set the initial positions
    comp_index = 0;
    comp_offset = 0;
    array_index = 0;
	
    while (comp_index < worker->bit_map_size[i]) {
      current = (current_map[comp_index] >> comp_offset) & 1;

      if (current) {
        worker->meth_array_reverse[i][array_index++].position = comp_offset + (comp_index * 8);
      }
      
      if (comp_offset == 7) {
        ++comp_index;
        comp_offset = 0;
      } else {
        ++comp_offset;
      }
    }
    
    // Free the bit arrays which won't be used anymore
    free(worker->bit_map_forward[i]);
    free(worker->bit_map_reverse[i]);
  }
  
  free(worker->bit_map_forward);
  free(worker->bit_map_reverse);
  
  worker->bit_to_map_time += omp_get_wtime() - start_time;
}
