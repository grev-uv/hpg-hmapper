/**
 * producer.c
 * - Date: 3 / Feb / 2017
 * - Who: Cesar Gonzalez (Cesar.Gonzalez-Segura@uv.es)
 *
 * The producer module is in charge of reading alignments from the BAM
 * file, filtering them if appropriate and possible and finally send
 * the alignments to their designated worker queue to be processed.
 *
 */

#include "producer.h"

#define BAM_ENTRY_DATA_SIZE   16384
static uint8_t BAM_ENTRY_DATA_FIELD[BAM_ENTRY_DATA_SIZE];

//-----------------------------------------------------

producer_input_t* producer_input_init(scheduler_input_t* scheduler, size_t thread_id) {
  const size_t chr_sz = scheduler->num_chromosomes * sizeof(int8_t);
  const size_t mth_sz = scheduler->num_chromosomes * sizeof(uint32_t);

  producer_input_t* producer_input = calloc(1, sizeof(producer_input_t));

  producer_input->num_chromosomes = scheduler->num_chromosomes;
  producer_input->num_threads = scheduler->num_threads;
  producer_input->thread_id = thread_id;
  
  producer_input->valid_reads = 0;
  producer_input->unmapped_reads = 0;
  producer_input->malformed_reads = 0;
  producer_input->temp_reads = 0;
  producer_input->memory_budget = scheduler->memory_budget;

  producer_input->in_queue = scheduler->worker_in_queue;

  producer_input->mc_bam = scheduler->mc_bam_file;
  producer_input->hmc_bam = scheduler->hmc_bam_file;
  producer_input->current_bam = MC_QUEUE_INDEX;

  producer_input->bam_entry = bam_init1();
  producer_input->bam_entry->m_data = BAM_ENTRY_DATA_SIZE;
  producer_input->bam_entry->data = BAM_ENTRY_DATA_FIELD;

  producer_input->num_reads[MC_QUEUE_INDEX] = 0;
  producer_input->num_reads[HMC_QUEUE_INDEX] = 0;

  producer_input->stage_complete[MC_QUEUE_INDEX] = 0;
  producer_input->stage_complete[HMC_QUEUE_INDEX] = 0;

  producer_input->batch_size = scheduler->batch_size;
  producer_input->pending_batches = array_list_new(producer_input->batch_size, 1.25F, COLLECTION_MODE_ASYNCHRONIZED);

  producer_input->chr_schedule = malloc(chr_sz);
  memcpy(producer_input->chr_schedule, scheduler->chr_schedule, chr_sz);

  // Copy the methylation statistics to use them in the filtering preprocess
  if (scheduler->methyl_reads) {
    producer_input->methyl_reads[MC_QUEUE_INDEX] = malloc(mth_sz);
    memcpy(producer_input->methyl_reads[MC_QUEUE_INDEX], scheduler->methyl_reads, mth_sz);
  } else {
    producer_input->methyl_reads[MC_QUEUE_INDEX] = NULL;
  }

  // Copy the hydroximethylation statistics to use them in the filtering preprocess
  if (scheduler->hmc_methyl_reads) {
    producer_input->methyl_reads[HMC_QUEUE_INDEX] = malloc(mth_sz);
    memcpy(producer_input->methyl_reads[HMC_QUEUE_INDEX], scheduler->hmc_methyl_reads, mth_sz);
  } else {
    producer_input->methyl_reads[HMC_QUEUE_INDEX] = NULL;
  }

  producer_input->num_reads[MC_QUEUE_INDEX] = calloc(producer_input->num_chromosomes, sizeof(size_t));
  producer_input->num_reads[HMC_QUEUE_INDEX] = calloc(producer_input->num_chromosomes, sizeof(size_t));

  producer_input->total_reads[MC_QUEUE_INDEX] = calloc(producer_input->num_chromosomes, sizeof(size_t));
  producer_input->total_reads[HMC_QUEUE_INDEX] = calloc(producer_input->num_chromosomes, sizeof(size_t));

  producer_input->time_sec = 0.0;

  return producer_input;
}

//-----------------------------------------------------

void producer_input_free(producer_input_t* producer_input) {
  if (producer_input) {
    bam_destroy1_nodealloc(producer_input->bam_entry);

    if (producer_input->chr_schedule) {
      free(producer_input->chr_schedule);
    }

    if (producer_input->pending_batches) {
      array_list_free(producer_input->pending_batches, NULL);
    }

    if (producer_input->methyl_reads[MC_QUEUE_INDEX]) {
      free(producer_input->methyl_reads[MC_QUEUE_INDEX]);
    }

    if (producer_input->methyl_reads[HMC_QUEUE_INDEX]) {
      free(producer_input->methyl_reads[HMC_QUEUE_INDEX]);
    }

    if (producer_input->num_reads[MC_QUEUE_INDEX]) {
      free(producer_input->num_reads[MC_QUEUE_INDEX]);
    }

    if (producer_input->num_reads[HMC_QUEUE_INDEX]) {
      free(producer_input->num_reads[HMC_QUEUE_INDEX]);
    }

    if (producer_input->total_reads[MC_QUEUE_INDEX]) {
      free(producer_input->total_reads[MC_QUEUE_INDEX]);
    }

    if (producer_input->total_reads[HMC_QUEUE_INDEX]) {
      free(producer_input->total_reads[HMC_QUEUE_INDEX]);
    }

    free(producer_input);
  }
} 

//-----------------------------------------------------

int producer_stage_step(void* data, scheduler_input_t* scheduler) {
  const size_t mem_budget = scheduler->memory_budget;
  size_t result = STAGE_IN_PROGRESS;
  int bam_result = 0;
  producer_input_t* input = (producer_input_t*)data;

  // Alignment and BAM file data
  bamFile fd;
  alignment_t* alignment = NULL;

  alignment_batch_t* current_batch = NULL;
  int8_t worker_id = 0;
  size_t type = input->current_bam;

  // Timing data
  double start_time = 0.0, start_time_read = 0.0;

  if (input->pass_id == PRODUCER_FIRST_PASS ||
     (input->pass_id == PRODUCER_SECOND_PASS && input->worker_team_status == WORKER_SECOND_PASS)) {
    // Read a single alignment from the current BAM file
    // Skip if the current BAM file has been processed or
    // if the memory budget has been reached
    if (!input->stage_complete[type] /*&& scheduler_get_memory_usage(scheduler) < input->memory_budget*/) {
      // Allocate space for a single alignment
      start_time_read = omp_get_wtime();

      // Read the alignment from the appropriate file and
      // convert it into an OpenCB alignment
      if (type == MC_QUEUE_INDEX) {
        fd = input->mc_bam->bam_fd;
      } else if (type == HMC_QUEUE_INDEX) {
        fd = input->hmc_bam->bam_fd;
      }

      bam_result = bam_read1_nonalloc(fd, input->bam_entry);

      // Check if the complete BAM file has been read
      if (bam_result >= 0) {
        input->total_bytes[type] += bam_result;

        // There still are alignments left to read.
        // Convert the BAM entry into a OpenCB alignment
        alignment = alignment_new_by_bam(input->bam_entry, PHRED_BASE_QUALITY);

        // Update the average alignment length
        if (scheduler->alignment_length == 0.0f) {
          scheduler->alignment_length = alignment->length;
          scheduler->alignment_mean_size = alignment->data_length;
        }

        scheduler->alignment_length = 0.5f * (scheduler->alignment_length + (float)alignment->length);
        scheduler->alignment_mean_size = 0.5f * (scheduler->alignment_mean_size + (float)alignment->data_length);

        input->time_bam_read += omp_get_wtime() - start_time_read;
        
        // Process the alignment
        if (alignment->chromosome >= 0) {
          ++input->total_reads[type][alignment->chromosome];
        }
        
        alignment = producer_process_alignment(alignment, input, type, scheduler, input->pass_id);
        start_time = omp_get_wtime();

        // If the alignment wasn't filtered, find a batch for this
        // chromosome and insert the alignment
        if (alignment) {
          size_t found = 0;
          size_t batch_count = array_list_size(input->pending_batches);

          if (batch_count) {
            for (int i = batch_count - 1; i >= 0; --i) {
              current_batch = array_list_get(i, input->pending_batches);

              if (current_batch->chromosome == alignment->chromosome &&
                  current_batch->type == type) {
                found = 1;
                break;
              }
            }

            if (found) {
              static_queue_push(alignment, current_batch->alignments);
            } else {
              // There isn't any batch matching the criteria, create
              // a new one
              current_batch = alignment_batch_new(input->batch_size, type, alignment->chromosome);
              static_queue_push(alignment, current_batch->alignments);
              array_list_insert(current_batch, input->pending_batches);
            }
          } else {
            // There are no pending batches, create a new one.
            current_batch = alignment_batch_new(input->batch_size, type, alignment->chromosome);
            static_queue_push(alignment, current_batch->alignments);
            array_list_insert(current_batch, input->pending_batches);
          }
        }

        input->time_sec += omp_get_wtime() - start_time;
      } else {
        // The end of file has been reached. The producer
        // stage for this BAM file is complete
        input->stage_complete[type] = 1;
        
        if (type == MC_QUEUE_INDEX) {
          input->current_bam = HMC_QUEUE_INDEX;
        }
      }
    }
  }

  // Check if there are pending batches, enqueue them into
  // the worker queues if it's the case
  start_time = omp_get_wtime();

  for (int i = array_list_size(input->pending_batches) - 1; i >= 0; --i) {
    current_batch = array_list_get(i, input->pending_batches);
    worker_id = input->chr_schedule[current_batch->chromosome];

    if (static_queue_is_full(current_batch->alignments)) {
      static_queue_push(current_batch, input->in_queue[worker_id]);
      array_list_remove_at(i, input->pending_batches);
    }
  }

  // Check if the completion flag is set for both
  // BAM files. If it is, flush all the pending 
  // batches to the workers.
  if (input->stage_complete[MC_QUEUE_INDEX] &&
      input->stage_complete[HMC_QUEUE_INDEX]) {
    for (int16_t i = array_list_size(input->pending_batches) - 1; i >= 0; --i) {
      current_batch = array_list_get(i, input->pending_batches);
      worker_id = input->chr_schedule[current_batch->chromosome];

      static_queue_push(current_batch, input->in_queue[worker_id]);
      array_list_remove_at(i, input->pending_batches);
    }

    printf("[Prod] Completed pass stats:\n");
    printf("       - Valid:     %lu reads (%.02f%%)\n", input->valid_reads,
        100.0f*((float)input->valid_reads/(float)(input->valid_reads + input->malformed_reads + input->unmapped_reads)));
    printf("       - Unmapped:  %lu reads (%.02f%%)\n", input->unmapped_reads, 
        100.0f*((float)input->unmapped_reads/(float)(input->valid_reads + input->malformed_reads + input->unmapped_reads)));
    printf("       - Malformed: %lu reads (%.02f%%)\n", input->malformed_reads,
        100.0f*((float)input->malformed_reads/(float)(input->valid_reads + input->malformed_reads + input->unmapped_reads)));

    input->valid_reads = 0;
    input->unmapped_reads = 0;
    input->malformed_reads = 0;

    if (input->pass_id == PRODUCER_FIRST_PASS) {
      #ifdef DEBUG
      printf("[Prod] Producer first pass completed.\n");
      #endif

      input->stage_complete[MC_QUEUE_INDEX] = 0;
      input->stage_complete[HMC_QUEUE_INDEX] = 0;
      input->current_bam = MC_QUEUE_INDEX;

      producer_reopen_bam_files(input, scheduler);

      input->pass_id = PRODUCER_SECOND_PASS;
      scheduler_set_producer_pass_status(scheduler, PRODUCER_SECOND_PASS);
    } else if (input->pass_id == PRODUCER_SECOND_PASS) {
      #ifdef DEBUG
      printf("[Prod] Producer second pass completed.\n");
      #endif

      scheduler_set_producer_pass_status(scheduler, PRODUCER_FINISH_PASS);
      result = STAGE_COMPLETED;
    }
  }

  input->time_sec += omp_get_wtime() - start_time;
  return result;
}

//-----------------------------------------------------

// Lambda function to free the tag objects
void __producer_remove_tag_lambda(void* data) {
  if (data) bam_tag_free((bam_tag_t*) data);
}

//-----------------------------------------------------

alignment_t* producer_process_alignment(alignment_t* alignment, producer_input_t* input, size_t type, 
        scheduler_input_t* scheduler, size_t pass_id) {
  size_t sq_length = 0;
  size_t alig_size = alignment->data_length;
  uint8_t meth_test_passed = 1, found_xm = 0;

  const int chromosome = alignment->chromosome;
  const uint32_t *methyl_stats;
  
  bam_tag_t* current_tag = NULL;
  size_t *num_reads;

  double start_time_filter = omp_get_wtime();

  if (type == MC_QUEUE_INDEX) {
    methyl_stats = input->methyl_reads[MC_QUEUE_INDEX];
    num_reads = input->num_reads[MC_QUEUE_INDEX];
  } else if (type == HMC_QUEUE_INDEX) {
    methyl_stats = input->methyl_reads[HMC_QUEUE_INDEX];
    num_reads = input->num_reads[HMC_QUEUE_INDEX];
  }

  // Skip the alignment if it comes from an unmapped read or if it is malformed
  if (chromosome >= 0 && alignment->position >= 0) {
    const uint8_t worker_id = input->chr_schedule[chromosome];

    // Check if the first 2 nucleotides are valid
    // values (AGCTN)
    if (alignment->sequence &&
        (alignment->sequence[0] == 'A' || alignment->sequence[0] == 'C' || alignment->sequence[0] == 'G' || alignment->sequence[0] == 'T' || alignment->sequence[0] == 'N') &&
        (alignment->sequence[1] == 'A' || alignment->sequence[1] == 'C' || alignment->sequence[1] == 'G' || alignment->sequence[1] == 'T' || alignment->sequence[1] == 'N'))
    {
      // If the alignment has BAM tags, extract them to check if
      // the ZM tag is present and if it is, check if there are
      // any methylated C's in the alignment. If there aren't any,
      // skip the alignment
      if (alignment->optional_tags) {
        // Skip alignments without all the optional mcontext tags
        if (array_list_size(alignment->optional_tags) == 7) {
          // Search for the ZM tag
          for (size_t i = 0; i < array_list_size(alignment->optional_tags); ++i) {
            current_tag = array_list_get(i, alignment->optional_tags);

            // If it is present and there is not any methylated C in the
            // alignment, mark the methylation test as failed. Also, if the
            // XM tag is not present, skip the alignment.
            if (current_tag) {
              if (!strcmp(current_tag->tag, XM_TAG_NAME)) {
                found_xm = 1;
              } else if (!strcmp(current_tag->tag, ZM_TAG_NAME)) {
                if (bam_tag_get_int(current_tag) == 0) {
                  meth_test_passed = 0;
                  break;
                }
              }
            }
          }
        }
      }

      if (!found_xm) {
        meth_test_passed = 0;
      }
    } else {
      // The tag data is badly formed. Skip the alignment
      meth_test_passed = 0;
    }

    // If the per-chromosome methylation statistics are available, skip
    // the alignment if all the methylated or hydroximethylated alignments 
    // have already been processed for this chromosome. Only process if the 
    // alignment wasn't discarded by the previous filter
    if (meth_test_passed && methyl_stats) {
      if (num_reads[chromosome] > methyl_stats[chromosome]) {
        meth_test_passed = 0;
      }
    }

    // If the alignment wasn't filtered, add it to the
    // input queue of the appropriate thread.
    if (meth_test_passed) {
      // Update the processed alignment counter
      ++num_reads[chromosome];
      alignment->methylation_type = type;
      ++input->valid_reads;
      ++input->temp_reads;

      if (input->temp_reads == SCHEDULER_MEMORY_EVENT_READ_SIZE) {
        input->temp_reads = 0;
        scheduler_alloc_event(scheduler, SCHEDULER_MEMORY_EVENT_READ_SIZE);
      }
    } else {
      if (chromosome == -1) {
        ++input->unmapped_reads;
      } else {
        ++input->malformed_reads;
      }
      
      alignment_free(alignment);
      alignment = NULL;
    }
  } else {
    if (chromosome == -1) {
      ++input->unmapped_reads;
    } else {
      ++input->malformed_reads;
    }

    alignment_free(alignment);
    alignment = NULL;
  }

  input->time_filter_alignment += omp_get_wtime() - start_time_filter;

  return alignment;
}

//-----------------------------------------------------

void producer_reopen_bam_files(producer_input_t* input, scheduler_input_t* scheduler) {
  bam_fclose(scheduler->mc_bam_file);
  bam_fclose(scheduler->hmc_bam_file);

  scheduler->mc_bam_file = bam_fopen_mode((char*)scheduler->mc_bam_path, NULL, "r");
  scheduler->hmc_bam_file = bam_fopen_mode((char*)scheduler->hmc_bam_path, NULL, "r");

  input->mc_bam = scheduler->mc_bam_file;
  input->hmc_bam = scheduler->hmc_bam_file;
}
