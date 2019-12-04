/**
 * scheduler.c
 * - Date: 3 / Feb / 2017
 * - Who: Cesar Gonzalez (Cesar.Gonzalez-Segura@uv.es)
 *
 * The scheduler is in charge of creating the producer, consumer and workers,
 * assigning the roles for each and synchronizing them.
 *
 */

#ifndef SCHEDULER_H
#define SCHEDULER_H

#include <omp.h>
#include <signal.h>

#include "bioformats/bam/bam_file.h"
#include "containers/static_queue.h"
#include "options.h"
#include "hmc_common.h"

#define STAGE_IN_PROGRESS           0    /**< The pipeline stage is working   */
#define STAGE_COMPLETED             1    /**< The pipeline stage has finished */

#define PRODUCER_FIRST_PASS         0    /**< ID used to signal the producer to pass to the first stage */
#define PRODUCER_SECOND_PASS        1    /**< ID used to signal the producer to pass to the second stage */
#define PRODUCER_FINISH_PASS        2    /**< ID used to signal the producer to finish */

#define WORKER_FIRST_PASS           0    /**< ID used to signal the workers to pass to the first stage */
#define WORKER_SECOND_PASS          1    /**< ID used to signal the workers to pass to the second stage */
#define WORKER_FINISH_PASS          2    /**< ID used to signal the workers to finish */

typedef int (*scheduler_producer_cb_t)(void*, void*);   /**< Callback signature for the producer thread */
typedef int (*scheduler_consumer_cb_t)(void*, void*);   /**< Callback signature for the consumer thread */
typedef int (*scheduler_worker_cb_t)(void*, void*);     /**< Callback signature for the worker thread   */

/**
 * @brief Global structure to store the scheduler state.
 * 
 */

typedef struct scheduler_input {
  size_t num_chromosomes;                 /**< Number of chromosomes in the reference genome */
  size_t num_threads;                     /**< Number of threads to use                      */
  size_t num_workers;                     /**< Number of worker processes to use             */
  uint64_t* chr_length;                   /**< Chrosome lengths from the reference genome    */
  int8_t** thr_schedule;                  /**< Per-thread chromosome scheduling matrix       */
  int8_t* chr_schedule;                   /**< Per-chromosome thread index                   */
  uint32_t* methyl_reads;                 /**< Number of reads with methylated C's per
                                               chromosome. NULL if not present               */
  uint32_t* hmc_methyl_reads;             /**< Number of reads with hydroximethylated C's per
                                                chromosome. NULL if not present              */
  bam_file_t* mc_bam_file;                /**< Input methylation BAM file                    */
  bam_file_t* hmc_bam_file;               /**< Input hydroximethylation BAM file             */
  char* mc_bam_path;
  char* hmc_bam_path;
  size_t treatment;                       /**< Hydroximethylation amplification treatment    */

  char* output_directory;                 /**< Output directory path                         */
  size_t memory_budget;                   /**< Maximum memory that the application should use
                                               in bytes. The application should never exceed
                                               this budget when allocating memory. Each module
                                               is responsible for doing so                   */

  int64_t memory_consumption;             /**< Memory being used by alignment data           */
  float alignment_mean_size;              /**< Mean size of an alignment in bytes            */
  float alignment_length;                 /**< Alignment string length in characters         */

  scheduler_producer_cb_t producer_cb;    /**< Pointer to the workflow producer callback     */
  scheduler_consumer_cb_t consumer_cb;    /**< Pointer to the workflow consumer callback     */
  scheduler_worker_cb_t worker_cb;        /**< Pointer to the workflow producer callback     */

  struct producer_input* producer_input;  /**< Pointer to the data input structure used by the producer thread  */
  struct consumer_input* consumer_input;  /**< Pointer to the data input structure used by the consumer thread  */
  struct worker_input** worker_input;     /**< Pointers to the data input structures used by the worker threads */

  size_t producer_status;                 /**< Status of the producer (STAGE_IN_PROGRESS or STAGE_COMPLETED)    */
  size_t producer_pass_status;
  size_t consumer_status;                 /**< Status of the consumer (STAGE_IN_PROGRESS or STAGE_COMPLETED)    */
  size_t* worker_status;                  /**< Status of the workers (STAGE_IN_PROGRESS or STAGE_COMPLETED)     */
  size_t worker_team_status;              /**< Joined status of all the workers (in progress if smaller than the
                                               number of workers in the team.                                   */
  size_t worker_team_pass_status;         /**< Joint stage status of the whole worker team   */

  static_queue_t** worker_in_queue;       /**< Worker input queues, one per worker thread. Workers fetch
                                               data from this queue                          */

  size_t mc_processed_reads;              /**< Used for global statistics                    */
  size_t hmc_processed_reads;             /**< Used for global statistics                    */

  size_t batch_size;                      /**< Batch size for the inter-stage queues         */
  size_t output_type;                     /**< Output format for the global statistics 
                                               (text or CSV, check hmc_common.h)             */

  char csv_delimiter;                     /**< User-defined delimiter for the CSV columns    */
  char csv_record_delimiter;              /**< User-defined delimiter for the CSV rows       */
  
  size_t quality_cutoff;

  size_t coverage;                        /**< User-defined minimum coverage for each methylated position       */
} scheduler_input_t;



/* ============================================================
 *                      PUBLIC MEMBERS
 * ============================================================*/

/**
*  @brief Initialize a new scheduler object
*  @param[in] index_path         Path to the directory containing the BWT index, NULL if missing.
*  @param[in] methyl_read_path   Path to the global methylation statistics file, NULL if missing.
*  @param[in] hmc_read_path      Path to the global hydroximethylation statistics file, NULL if missing.
*  @param[in] mc_bam_path        Path to the methylation BAM file.
*  @param[in] hmc_bam_path       Path to the hydroximethylation BAM file.
*  @param[in] output_dir         Path to the output directory.
*  @param[in] num_threads        Total number of threads to use in the pool (combined workers + producer/consumer).
*  @param[in] treatment          Treatment used for the hmC amplification (HMC_MARKER_TAB_SEQ or HMC_MARKER_OX_SEQ).
*  @param[in] timing             Pointer to the timing module, NULL if missing.
*  @param[in] memory_budget      Maximum memory to be used by the thread pool, in bytes.
*  @return                       Pointer to the new scheduler object.
*  
*/
scheduler_input_t* scheduler_input_init(const char* index_path,
                                        const char* methyl_read_path,
                                        const char* hmc_read_path, 
                                        const char* mc_bam_path,
                                        const char* hmc_bam_path,
                                        const char* output_dir,
                                        size_t num_threads,
                                        size_t treatment,
                                        size_t memory_budget, 
                                        size_t batch_size, 
                                        size_t output_type, 
                                        char csv_delimiter, 
                                        char csv_record_delimiter, 
                                        size_t quality_cutoff, 
                                        size_t coverage);

/**
*  @brief Deallocate a scheduler object.
*  @param[in] scheduler   Pointer to the scheduler object.
*  
*/
void scheduler_input_free(scheduler_input_t* scheduler);

/**
*  @brief Run the scheduler with the given workflow stage callbacks. The number of threads is
*         specified in the scheduler constructor.
*  @param[in]    scheduler     Pointer to the scheduler object.
*  @param[in]    producer_cb   Callback to the producer stage function.
*  @param[in]    consumer_cb   Callback to the consumer stage function.
*  @param[in]    worker_cb     Callback to the worker stage function.
*  
*/
void scheduler_run(scheduler_input_t* scheduler, scheduler_producer_cb_t producer_cb, scheduler_consumer_cb_t consumer_cb,
                        scheduler_worker_cb_t worker_cb);

/**
*  @brief Notice the memory manager that a read batch has been allocated
*  @param[in] scheduler     Pointer to the scheduler object.
*  @param[in] size          Number of allocated reads.
*  
*/
void scheduler_alloc_event(scheduler_input_t* scheduler, int64_t size);

/**
*  @brief Notice the memory manager that a data structure has been allocated
*  @param[in] scheduler     Pointer to the scheduler object.
*  @param[in] size          Allocated size in bytes.
*  
*/
void scheduler_alloc_event_direct(scheduler_input_t* scheduler, int64_t size);

/**
*  @brief Notice the memory manager that a read batch has been freed
*  @param[in] scheduler     Pointer to the scheduler object.
*  @param[in] size          Number of freed reads.
*  
*/
void scheduler_free_event(scheduler_input_t* scheduler, int64_t size);

/**
*  @brief Notice the memory manager that a data structure has been freed
*  @param[in] scheduler     Pointer to the scheduler object.
*  @param[in] size          Freed size in bytes.
*  
*/
void scheduler_free_event_direct(scheduler_input_t* scheduler, int64_t size);

/**
*  @brief Retrieve the memory used by user data in bytes
*  @param[in] scheduler     Pointer to the scheduler object.
*  @return                  Memory used by user data in bytes
*  
*/
int64_t scheduler_get_memory_usage(scheduler_input_t* scheduler);


/**
*  @brief Set the producer status 
*  @param[in] scheduler     Pointer to the scheduler object.
*  @param     status        Producer stage status (STAGE_COMPLETED or STAGE_IN_PROGRESS)
*  
*/
void scheduler_set_producer_status(scheduler_input_t* scheduler, size_t status);

/**
*  @brief Set the status for a given worker
*  @param[in] scheduler     Pointer to the scheduler object.
*  @param     worker_id     Worker stage status (STAGE_COMPLETED or STAGE_IN_PROGRESS)
*  
*/
size_t scheduler_set_worker_status(scheduler_input_t* scheduler, size_t worker_id, size_t status);

/**
*  @brief Print the execution statistics through stdout.
*  @param[in] scheduler     Pointer to the scheduler object.
*  
*/
void scheduler_print_statistics(scheduler_input_t* scheduler);

/**
*  @brief Notify the scheduler that a worker has completed its current stage.
*  @param[in] scheduler     Pointer to the scheduler object.
*  
*/
void scheduler_finished_worker_stage_event(scheduler_input_t* scheduler);

/**
*  @brief Set the consumer status 
*  @param[in] scheduler     Pointer to the scheduler object.
*  @param     status        Consumer stage status (STAGE_COMPLETED or STAGE_IN_PROGRESS)
*  
*/
void scheduler_set_producer_pass_status(scheduler_input_t* scheduler, size_t pass_id);



/* ============================================================
 *                      PRIVATE MEMBERS
 * ============================================================*/

/**
*  @brief Load the chromosome count and lengths from the BWT index of the reference genome.
*  @param[in]     index_path         Path to the directory containing the BWT index.
*  @param[in,out] scheduler          Pointer to the scheduler object.
*  
*  This is an internal function and should not be called by user code.
*/
void scheduler_load_chr_from_bwt_index(const char* index_path, scheduler_input_t* scheduler);

/**
*  @brief Generate the scheduling matrices using the global methylation statistics.
*  @param[in]     methyl_read_path   Path to the global methylation statistics file.
*  @param[in,out] scheduler          Pointer to the scheduler object.
*  
*  This is an internal function and should not be called by user code.
*/
void scheduler_gen_matrices_from_metyhl(const char* methyl_read_path, scheduler_input_t* scheduler);

/**
*  @brief Generate the scheduling matrices using a set of chromosome lengths.
*  @param[in]     sz_buffer   Buffer containing the target chromosome lengths.
*  @param[in,out] scheduler   Pointer to the scheduler object.
*  
*  This is an internal function and should not be called by user code.
*/
void scheduler_gen_matrices(uint64_t* sz_buffer, scheduler_input_t* scheduler);

/**
*  @brief Start the producer stage.
*  @param[in] scheduler   Pointer to the scheduler object.
*  
*  This is an internal function and should not be called by user code.
*/
void scheduler_producer_start(scheduler_input_t* scheduler);

/**
*  @brief Start the consumer stage.
*  @param[in] scheduler   Pointer to the scheduler object.
*  
*  This is an internal function and should not be called by user code.
*/
void scheduler_consumer_start(scheduler_input_t* scheduler);

/**
*  @brief Start a worker stage.
*  @param[in] scheduler   Pointer to the scheduler object.
*  @param[in] id          Worker identifier number.
*  
*  This is an internal function and should not be called by user code.
*/
void scheduler_worker_start(scheduler_input_t* scheduler, size_t id);

/**
*  @brief Run the pipeline using a single thread
*  @param[in] scheduler     Pointer to the scheduler object.
*  
*/
void scheduler_run_single_threaded(scheduler_input_t* scheduler);

#endif // SCHEDULER_H
