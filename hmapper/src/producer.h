/**
 * producer.h
 * - Date: 3 / Feb / 2017
 * - Who: Cesar Gonzalez (Cesar.Gonzalez-Segura@uv.es)
 *
 * The producer module is in charge of reading alignments from the BAM
 * file, filtering them if appropriate and possible and finally send
 * the alignments to their designated worker queue to be processed.
 *
 */

#ifndef PRODUCER_H
#define PRODUCER_H

#include <ctype.h>

#include "bioformats/bam/bam_tags.h"
#include "commons/commons.h"
#include "scheduler.h"


/*
 * producer_input struct
 * - Object containing all the state for the producer.
 */
typedef struct producer_input {
  bam_file_t* mc_bam;     /**< OpenCB BAM file pointer to the dataset containing the 5mC (BS-seq) reads */
  bam_file_t* hmc_bam;    /**< OpenCB BAM file pointer to the dataset containing the 5hmC (TAB-seq) reads */
  bam1_t* bam_entry;      /**< Temporary BAM entry used to hold the current BAM alignment data */
  size_t current_bam;     /**< Size of the current BAM entry */

  static_queue_t** in_queue;        /**< Input queues for all the workers */
  array_list_t* pending_batches;    /**< Alignment batches pending to be sent through the pipeline*/

  size_t num_chromosomes;   /**< Number of chromosomes in the current genome */
  size_t num_threads;       /**< Number of threads used to execute HPG-Hmapper */
  size_t thread_id;         /**< Global thread ID for the producer */

  size_t batch_size;          /**< Batch size for the worker input queues */
  size_t memory_budget;       /**< Global memory budget in bytes, set by the memory manager */
  size_t pass_id;             /**< Current stage of the producer*/
  size_t worker_team_status;  /**< Aggregate stage ID of the worker team */

  int8_t* chr_schedule;                       /**< Chromosome scheduling matrix */
  uint32_t* methyl_reads[METH_QUEUE_COUNT];   /**< Used for statistics purposes */

  size_t stage_complete[METH_QUEUE_COUNT];    /**< Used to track the current stage for both 5mC and 5hmC processing queues */
  size_t* num_reads[METH_QUEUE_COUNT];        /**< Used for statistic purposes */
  size_t* total_reads[METH_QUEUE_COUNT];      /**< Used for statistic purposes */
  size_t total_bytes[METH_QUEUE_COUNT];       /**< Used for statistic purposes */

  size_t valid_reads;       /**< Used for statistic purposes */
  size_t unmapped_reads;    /**< Used for statistic purposes */
  size_t malformed_reads;   /**< Used for statistic purposes */
  size_t temp_reads;        /**< Used for statistic purposes */

  double time_sec;                /**< Used for statistic purposes */
  double time_bam_read;           /**< Used for statistic purposes */
  double time_filter_alignment;   /**< Used for statistic purposes */
} producer_input_t;



/* ============================================================
 *                      PUBLIC MEMBERS
 * ============================================================*/

/**
*  @brief Initialize the producer object.
*  @param[in] scheduler     Pointer to the scheduler object.
*  @param     thread_id     Global thread ID for the producer.
*  @return                  Pointer to the new producer object.
*  
*/
producer_input_t* producer_input_init(scheduler_input_t* scheduler, size_t thread_id);

/**
*  @brief Free the resources used by the producer object.
*  @param[in,out]   consumer_input   Pointer to the producer object.
*  
*/
void producer_input_free(producer_input_t* scheduler);

/**
*  @brief Perform a step of the producer object.
*  @param[in,out]   data        Pointer to the input data (producer_input_t).
*  @param[in]       scheduler   Pointer to the scheduler object.
*  @return                      State of the producer after the processing (check scheduler.h).
*  
*/
int producer_stage_step(void* data, scheduler_input_t* scheduler);



/* ============================================================
 *                      PRIVATE MEMBERS
 * ============================================================*/


/**
*  @brief Process a single alignment. This includes filtering unmapped and unmethylated reads,
*         as well as detecting malformed reads.
*  @param[in,out]   alignment   Pointer to the input read.
*  @param[in]       input       Pointer to the producer object.
*  @param           type        Read type (5mC or 5hmC, check hmc_common.h).
*  @param[in]       scheduler   Pointer to the scheduler object.
*  @param           pass_id     ID of the current producer stage
*  @return                      The pointer to the alignment object if it was deemed as valid, NULL if
*                               the alignment was filtered.
*  
*/
alignment_t* producer_process_alignment(alignment_t* alignment, producer_input_t* input, size_t type, scheduler_input_t* scheduler, size_t pass_id);

/**
*  @brief Reset the BAM file descriptors and set up the auxiliary structures for the next stage.
*  @param[in]       input       Pointer to the producer object.
*  @param[in]       scheduler   Pointer to the scheduler object.
*  
*/
void producer_reopen_bam_files(producer_input_t* input, scheduler_input_t* scheduler);

#endif // PRODUCER_H
