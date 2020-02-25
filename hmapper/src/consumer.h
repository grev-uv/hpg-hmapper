/**
 * consumer.h
 * - Date: 17 / Jul / 2017
 * - Who: Cesar Gonzalez (Cesar.Gonzalez-Segura@uv.es)
 *
 * The consumer module is in charge of writing the methylation maps
 * to the disk using the format specified in the options.
 *
 */

#ifndef CONSUMER_H
#define CONSUMER_H

#include "scheduler.h"
#include "worker.h"

/*
 * consumer_input struct
 * - Object containing all the state for the consumer. Direct access
 *   to the object members is only acceptable for the public members.
 *   Atomic access is not needed for reading. Write access to the members
 *   should be done only through the object's public methods.
 */
typedef struct consumer_input {
  size_t num_chromosomes;       /**< Number of chromosomes in the current genome. */
  size_t num_workers;           /**< Number of active workers */
  size_t thread_id;             /**< Global thread ID of the consumer */

  FILE** meth_map_forward_fd;   /**< PRIVATE File descriptors to the output forward maps */
  FILE** meth_map_reverse_fd;   /**< PRIVATE File descriptors to the output reverse maps */
  FILE** meth_map_mix_fd;       /**< PRIVATE File descriptors to the output forward reverse mixed maps */

  size_t total_c;               /**< PRIVATE Used for statistics purposes by the scheduler object */
  size_t total_mc;              /**< PRIVATE Used for statistics purposes by the scheduler object */
  size_t total_hmc;             /**< PRIVATE Used for statistics purposes by the scheduler object */

  size_t worker_team_pass;      /**< PRIVATE Aggregated stage ID of the worker team. Used to signal
                                     when to start logging to the output files (check scheduler.h). */
  char csv_delimiter;           /**< PRIVATE User-specified delimiter for each column in the CSV files */
  char csv_record_delimiter;    /**< PRIVATE User-specified delimiter for each row in the CSV files */

  double time_sec;              /**< PRIVATE Used for timing purposes */
  
  size_t coverage;              /**< Minimum number of reads for each position to write in a file*/
} consumer_input_t;



/* ============================================================
 *                      PUBLIC MEMBERS
 * ============================================================*/

/**
*  @brief Initialize the consumer object.
*  @param[in] scheduler     Pointer to the scheduler object.
*  @param     thread_id     Global thread ID for the consumer.
*  @return                  Pointer to the new consumer object.
*  
*/
consumer_input_t* consumer_input_init(scheduler_input_t* scheduler, size_t thread_id);

/**
*  @brief Free the resources used by the consumer object.
*  @param[in,out]   consumer_input   Pointer to the consumer object.
*  
*/
void consumer_input_free(consumer_input_t* consumer_input);

/**
*  @brief Perform a step of the consumer object.
*  @param[in,out]   data        Pointer to the input data (consumer_input_t).
*  @param[in]       scheduler   Pointer to the scheduler object.
*  @return                      State of the consumer after the processing (check scheduler.h).
*  
*/
int consumer_stage_step(void* data, scheduler_input_t* scheduler);



/* ============================================================
 *                      PRIVATE MEMBERS
 * ============================================================*/

/**
*  @brief Serialize a methylation map to a CSV file. Called by consumer_stage_step.
*  @param[in]   array             Pointer to the methylation map.
*  @param       length            Length of the methylation map.
*  @param[in]   fd                Descriptor of the output file.
*  @param       delimiter         Delimiter for the CSV columns.
*  @param       record_delimiter  Delimiter for the CSV rows.
*  
*/
void consumer_meth_array_serialize(meth_array_node_t* array, 
                                   size_t length, 
                                   FILE* fd, 
                                   char delimiter, 
                                   char record_delimiter,
                                   size_t coverage);

void consumer_meth_array_serialize_mix(meth_array_node_t* array_f, 
                                       size_t length_f, 
                                       meth_array_node_t* array_r, 
                                       size_t length_r, 
                                       FILE* fd, 
                                       char delimiter, 
                                       char record_delimiter,
                                       size_t coverage);

#endif // CONSUMER_H
