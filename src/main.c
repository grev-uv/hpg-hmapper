/**
 * main.c
 * - Date: 1 / Dec / 2016
 * - Who: Cesar Gonzalez (Cesar.Gonzalez-Segura@uv.es)
 *
 * Entry point for the hmc-generator application.
 *
 */
#include "options.h"
#include "scheduler.h"
#include "producer.h"
#include "consumer.h"
#include "worker.h"

int main(int argc, char** argv) {
  options_t options;
  scheduler_input_t* scheduler = NULL;

  // Read the command line options
  if (options_read_cmd(argc, argv, &options)) {
    exit(1);
  }

  // Create the scheduler and launch the pipeline
  scheduler = scheduler_input_init(options.index_directory, options.mc_methyl_stats_file,
                options.hmc_methyl_stats_file, options.mc_bam_file, options.hmc_bam_file, 
                options.output_directory, options.num_threads, options.treatment,
                options.memory_budget, options.batch_size, options.stats_output_format, 
                options.csv_delimiter, options.csv_record_delimiter);

  // Start the scheduler
  scheduler_run(scheduler, (scheduler_producer_cb_t)producer_stage_step, 
                  (scheduler_consumer_cb_t)consumer_stage_step, 
                  (scheduler_worker_cb_t)worker_stage_step);

  // Clean-up
  scheduler_input_free(scheduler);
  exit(0);
}