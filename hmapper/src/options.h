/**
 * options.h
 * - Date: 1 / Feb / 2017
 * - Who: Cesar Gonzalez (Cesar.Gonzalez-Segura@uv.es)
 *
 * The options module is in charge of parsing the command line
 * options and generating the set-up for running the application.
 *
 */

#ifndef OPTIONS_H
#define OPTIONS_H

#include <omp.h>
#include <sys/sysinfo.h>

#include "commons/commons.h"
#include "commons/file_utils.h"
#include "hmc_common.h"

/**
 * @brief Command line input errors
 *
 */
#define OPTIONS_ENOFILE   0   /**< Argument is not a file */
#define OPTIONS_ENOTBAM   1   /**< File is not a BAM file */
#define OPTIONS_ENOTCSV   2   /**< File is not a CSV file */
#define OPTIONS_ENOTCMD   3   /**< Unknown command line argument */
#define OPTIONS_ENOTBIN   4   /**< File is not a BIN file */
#define OPTIONS_EBADCMD   5   /**< Invalid argument option */
#define OPTIONS_ENOTDIR   6   /**< Argument is not a directory */
#define OPTIONS_ETOOFEW   7   /**< Mandatory options have not been set*/
#define OPTIONS_ETHREAD   8   /**< Tried to use less than three threads */
#define OPTIONS_EBUDGET   9   /**< Memory budget is over the system memory */

/**
 * @brief Command line argument strings
 *
 */
#define TREATMENT_CMD_STRING              "--treatment"              /**< Read treatment (TAB-Seq, OX-Seq) */
#define MC_FILE_CMD_STRING                "--mc-bam-input-file"      /**< Methylation input file */
#define MC_FILE_CMD_STRING_SHORT          "-mc"                      /**< Methylation input file */
#define HMC_FILE_CMD_STRING               "--hmc-bam-input-file"     /**< Hydroximethylation input file */
#define HMC_FILE_CMD_STRING_SHORT         "-hmc"                     /**< Hydroximethylation input file */
#define TIMING_CMD_STRING                 "--timing"                 /**< Profile timing */
#define TIMING_CMD_STRING_SHORT           "-t"                       /**< Profile timing */
#define NUM_THREADS_CMD_STRING            "--cpu-threads"            /**< Thread count */
#define INDEX_CMD_STRING                  "--bwt-index"              /**< Path to the genome index */
#define INDEX_CMD_STRING_SHORT            "-i"                       /**< Path to the genome index */
#define MC_METH_STATS_CMD_STRING          "--mc-methyl-stats-file"   /**< Methylation per-chromosome stats file */
#define HMC_METH_STATS_CMD_STRING         "--hmc-methyl-stats-file"  /**< Hydroximethylation per-chromosome stats file */
#define OUTPUT_PATH_CMD_STRING            "--output"                 /**< Output directory */
#define OUTPUT_PATH_CMD_STRING_SHORT      "-o"                       /**< Output directory */
#define MEMORY_BUDGET_CMD_STRING          "--memory"                 /**< Memory budget string (#.#M or #.#G) */
#define BATCH_SIZE_CMD_STRING             "--batch-size"
#define STATS_OUTPUT_CMD_STRING           "--output-format"
#define CSV_DELIMITER_CMD_STRING          "--csv-delimiter"
#define CSV_RECORD_DELIMITER_CMD_STRING   "--csv-record-delimiter"
#define QUALITY_CUTOFF_CMD_STRING         "--quality"

#define OUTPUT_TYPE_TEXT_STR              "text"
#define OUTPUT_TYPE_CSV_STR               "csv"

/**
 * @brief Structure to store the application options.
 * 
 */
typedef struct options {
  char hmc_bam_file[MAX_FILENAME_LENGTH];           /**< Input hydroximethylated BAM file */
  char mc_bam_file[MAX_FILENAME_LENGTH];            /**< Input methylated BAM file */
  char output_directory[MAX_FILENAME_LENGTH];       /**< Output directory path */
  char index_directory[MAX_FILENAME_LENGTH];        /**< Genome index path */
  char mc_methyl_stats_file[MAX_FILENAME_LENGTH];   /**< Input per-chromosome methylation stats file */
  char hmc_methyl_stats_file[MAX_FILENAME_LENGTH];  /**< Input per-chromosome hydroximethylation stats file */

  size_t treatment;       /**< Marker used for the hydroximethylation treatment */
  size_t timing;          /**< Is timing enabled? */
  size_t num_threads;     /**< Thread count */
  size_t memory_budget;   /**< Maximum ammount of memory in bytes the application can use */

  size_t batch_size;            /**< Batch size for the inter-stage queues */
  size_t stats_output_format;   /**< Output format for the global statistics (check hmc_common.h) */

  char csv_delimiter;         /**< User-specified delimiter for the output CSV columns */
  char csv_record_delimiter;  /**< User-specified delimiter for the output CSV rows */

  size_t quality_cutoff;
} options_t;

#define MIN_OPTION_COUNT    2   /**< Minimum number of command line arguments */

/**
*  @brief Reads and validates the command line options.
*  @param         argc    Command line argument count.
*  @param         argv    Pointer to the command line argument list.
*  @param[in,out] options Pointer to the options list.
*  @return 0 on success, -1 on failure.
*  
*/
int options_read_cmd(const int argc, char** argv, options_t* options);

/**
*  @brief Prints the error message associated through stdout
*         with a given internal error.
*  @param err Error code.
*  
*/
void options_error(const int err);

/**
*  @brief Prints the application usage through stdout.
*  
*/
void options_print_usage();

#endif // OPTIONS_H
