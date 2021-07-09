/**
 * worker.h
 * - Date: 15 / Feb / 2017
 * - Who: Cesar Gonzalez (Cesar.Gonzalez-Segura@uv.es)
 *
 * Workers are in charge of retrieving alignments from the input queue, processing them
 * to create the methylation and hydroximethylation map and storing the result into the
 * output queue.
 *
 */

#ifndef WORKER_H
#define WORKER_H

#include "scheduler.h"

/*
 * meth_array_node struct
 * - Used to store the methylation data in the second stage bit arrays.
 */
typedef struct __attribute__((__packed__)) meth_array_node {
  uint32_t position;    /* Position in the chromosome */
  uint16_t c_count;      /* Number of non-methylated cytosines */
  uint16_t nc_count;     /* Number of overlaping non-cytosines (5mC)*/
  uint16_t mc_count;     /* Number of 5-mC cytosines   */
  uint16_t ch_count;     /* Number of non-hydroxymethylated cytosines */
  uint16_t nch_count;    /* Number of overlaping non-cytosines (5hmC)*/
  uint16_t hmc_count;    /* Number of 5-hmC cytosines  */
} meth_array_node_t;


/*
 * worker_input struct
 * - Object containing all the state for a single worker. Direct access
 *   to the object members is only acceptable for the public members.
 *   Atomic access is not needed for reading. Write access to the members
 *   should be done only through the object's public methods.
 */
typedef struct worker_input {
  static_queue_t* in_queue;   /**< Input queue for the incoming reads. Used by the scheduler object */
  int8_t* thr_schedule;       /**< Thread scheduling matrix for the local worker */

  size_t num_chromosomes;   /**< Total number of chromosomes in the current genome */
  size_t thread_id;         /**< Global thread ID for this worker */
  size_t worker_id;         /**< Local worker ID */
  
  size_t batch_size;        /**< PRIVATE Batch size for the input queues */
  size_t temp_reads;        /**< PRIVATE Used for debug purposes */
  size_t processed_reads;   /**< PRIVATE Used for debug purposes */
  size_t pass_id;           /**< PRIVATE Index of the stage being executed currently */
  size_t next_pass;         /**< PRIVATE Index of the stage that must be executed when the current stage is completed */

  uint8_t** bit_map_forward;  /**< PRIVATE Pointers to the forward bit arrays, used in the first stage */
  uint8_t** bit_map_reverse;  /**< PRIVATE Pointers to the reverse bit arrays, used in the first stage */
  size_t* bit_map_size;       /**< PRIVATE Length of the bit array pointers, used in the first stage */

  meth_array_node_t** meth_array_forward;   /**< PRIVATE Pointers to the forward methylation maps, used in the second stage */
  meth_array_node_t** meth_array_reverse;   /**< PRIVATE Pointers to the reverse methylation maps, used in the second stage */
  size_t* meth_array_forward_size;  /**< PRIVATE Length of the forward methylation maps, used in the second stage */
  size_t* meth_array_reverse_size;  /**< PRIVATE Length of the reverse methylation maps, used in the second stage */

  double time_sec;              /**< PRIVATE Used for timing purposes */ 
  double first_stage_process;   /**< PRIVATE Used for timing purposes */ 
  double second_stage_process;  /**< PRIVATE Used for timing purposes */ 
  double bit_to_map_time;       /**< PRIVATE Used for timing purposes */ 
} worker_input_t;



/* ============================================================
 *                      PUBLIC MEMBERS
 * ============================================================*/

/**
*  @brief Initialize a worker object.
*  @param[in] scheduler     Pointer to the scheduler object.
*  @param     thread_id     Global thread ID for this worker.
*  @return                  Pointer to the new worker object.
*  
*/
worker_input_t* worker_input_init(scheduler_input_t* scheduler, size_t thread_id);

/**
*  @brief Free the resources used by the worker object.
*  @param[in,out]   worker_input   Pointer to the worker object.
*  
*/
void worker_input_free(worker_input_t* worker_input);

/**
*  @brief Perform a step of the current worker object.
*  @param[in,out]   data        Pointer to the input data (worker_input_t).
*  @param[in]       scheduler   Pointer to the scheduler object.
*  @return                      State of the worker after the processing (check scheduler.h).
*  
*/
int worker_stage_step(void* data, scheduler_input_t* scheduler);

/**
*  @brief Report the worker that the producer has completed a stage, to
*         let it modify its state accordingly. Must be called from a critical
*         section.
*  @param[in,out]   worker    Pointer to the worker object.
*  @param           pass_id   Target pass for the worker.
*/
void worker_pass_finished_event(worker_input_t* worker, size_t pass_id);

/**
*  @brief Perform a blocking query to find a methylation node in the given methylation
*         map, located in the passed position in the genome.
*  @param[in]      array              Pointer to the methylation map object.
*  @param          length             Length of the methylation map.
*  @param          position           Target position in the reference genome.
*  @param[out]     current_position   Position in the methylation map where the genome position is located.
*/
meth_array_node_t* worker_meth_array_find(meth_array_node_t* array, size_t length, uint64_t position, size_t* current_position);



/* ============================================================
 *                      PRIVATE MEMBERS
 * ============================================================*/

/**
*  @brief Process a single alignment. Called by worker_stage_step.
*  @param[in]    worker         Pointer to the worker object.
*  @param[in]    alignment      Pointer to the aligned read.
*  @param        type           Methylation type (5mC or 5hmC, see hmc_common.h).
*  @param[in]    scheduler      Pointer to the scheduler object.
*  @param        pass_id        ID of the current stage (see scheduler.h).
*/
void worker_process_alignment(worker_input_t* worker, alignment_t* alignment, size_t type, scheduler_input_t* scheduler, size_t pass_id);

/**
*  @brief Start a blocking call to convert the methylation bit arrays of a worker to
*         methylation maps. The bit array data will be deleted. Called by worker_stage_step
*         when passing from the first to the second stage (check scheduler.h).
*  @param[in]    worker         Pointer to the worker object.
*/
void worker_meth_map_to_array(worker_input_t* worker);

#endif // WORKER_H
