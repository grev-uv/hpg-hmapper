/**
 * options.h
 * - Date: 1 / Feb / 2017
 * - Who: Cesar Gonzalez (Cesar.Gonzalez-Segura@uv.es)
 *
 * The options module is in charge of parsing the command line
 * options and generating the set-up for running the application.
 *
 */

#include "options.h"

//------------------------------------------------------------------------------------

void options_error(const int err) {
  switch (err) {
    case OPTIONS_ENOFILE:
      printf("Error: argument must be a file.\n\n");
      break;

    case OPTIONS_ENOTBAM:
      printf("Error: argument must be a BAM file.\n\n");
      break;
    
    case OPTIONS_ENOTCSV:
      printf("Error: argument must be a CSV file.\n\n");
      break;
    
    case OPTIONS_ENOTCMD:
      printf("Error: an invalid argument was introduced.\n");
      break;
    
    case OPTIONS_ENOTBIN:
      printf("Error: argument must be a BIN file.\n");
      break;

    case OPTIONS_ENOTDIR:
      printf("Error: argument must be an existing directory.\n");
      break;

    case OPTIONS_EBADCMD:
      printf("Error: the specified argument is not valid.\n");
      break;

    case OPTIONS_ETOOFEW:
      printf("Error: mandatory options have not been set.\n");
      break;

    case OPTIONS_ETHREAD:
      printf("Error: A minimum of three threads are required to run the application.\n");
      break;

    case OPTIONS_EBUDGET:
      printf("Error: The entered memory budget is bigger than the available system memory.\n");
      break;
  }
}

//------------------------------------------------------------------------------------

void options_print_usage() {
  printf("hmc-report\n");
  printf("usage:\n");
  printf("hmc-report -mc|--mc-bam-input-file=<file> -hmc|--hmc-bam-input-file=<file>");
  printf(" [--treatment=<int>] [-t|--timing] [--num-threads=<int>] [-i|--bwt-index=<path>] ");
  printf("[--mc-methyl-stats-file=<file>] [--hmc-methyl-stats-file=<file>] [-o|--output=<path>]\n");

  printf("-mc, --mc-bam-input-file\t\tBisulphite treated BAM input\n");
  printf("-hmc, --hmc-bam-input-file\t\tTAB-Seq / OX-Seq treated BAM input\n");
  printf("--num-threads\t\t\t\tNumber of threads to use (default is the num. of available threads in the system\n");
  printf("-i, --bwt-index\t\t\t\tPath to the HPG-Methyl BWT index (if not specified human GRCh37.68 is used as default)\n");
  printf("--mc-methyl-stats-file\t\t\tPer-chromosome methylation stats file (if not specified genome ordering is used as default and mC filtering is disabled)\n");
  printf("--hmc-methyl-stats-file\t\t\tPer-chromosome hydroximethylation stats file (if not specified hmC filtering is disabled)\n");
  printf("-o, --output\t\t\t\tOutput directory (default is working directory)\n");
  printf("--memory\t\t\t\tMemory budget to assign to the program (X.XXM/G). Default is 80%% of system memory.\n");
  printf("--output-format [text|csv]\t\t\t\tStatistics output format (plain text or CSV)\n");
  printf("--batch-size\t\t\t\tBatch size for the worker queues.\n");
  printf("--csv-delimiter\t\t\t\tDelimiter for the CSV columns\n");
  printf("--csv-record-delimiter\t\t\t\tDelimiter for the CSV rows\n");
  printf("--quality\t\t\tMinimum quality threshold (QUAL) for the input reads\n");
}

//------------------------------------------------------------------------------------

int options_read_cmd(const int argc, char** argv, options_t* options) {
  struct sysinfo sinfo;
  int error = 0;
  
  memset(options, 0, sizeof(options_t));

  // Initialize options with default values
  options->treatment = HMC_MARKER_TAB_SEQ;
  options->num_threads = omp_get_num_procs();

  options->batch_size = 1000000;
  options->quality_cutoff = 20;
  options->stats_output_format = STATS_OUTPUT_FORMAT_TEXT;

  options->csv_delimiter = ' ';
  options->csv_record_delimiter = '\n';

  memset(options->hmc_bam_file, 0, MAX_FILENAME_LENGTH);
  memset(options->mc_bam_file, 0, MAX_FILENAME_LENGTH);
  memset(options->index_directory, 0, MAX_FILENAME_LENGTH);
  memset(options->mc_methyl_stats_file, 0, MAX_FILENAME_LENGTH);
  memset(options->hmc_methyl_stats_file, 0, MAX_FILENAME_LENGTH);
  getcwd(options->output_directory, MAX_FILENAME_LENGTH);

  // Get the available system memory and set the default
  // memory budget to the 80% of the total memory
  sysinfo(&sinfo);
  options->memory_budget = 0.8 * sinfo.totalram * sinfo.mem_unit;

  // Check if the user has introduced at least
  // one command line argument
  if (argc < 1 + MIN_OPTION_COUNT) {
    error = 1;
    options_error(OPTIONS_ENOFILE);
  } else {
    for (size_t i = 1; i < argc; ++i) {
      // Compare each input argument with all possible
      // commands
      //
      // Memory budget
      if (!strcmp(argv[i], MEMORY_BUDGET_CMD_STRING)) {
        if (i + 1 < argc) {
          // Try to parse the memory budget
          // The first argument should be a decimal number
          // with the budget size, the second argument the
          // units (M(egabytes) or G(igabytes))
          double budget = 0.0;
          char units[3];

          sscanf(argv[i + 1], "%lf%s", &budget, units);

          // Check if the selected units are valid and that
          // the entered memory budget is inside the machine
          // maximum memory
          if (units[0] == 'M') {
            options->memory_budget = 1048576 * (size_t)budget;
            ++i;
          } else if (units[0] == 'G') {
            options->memory_budget = 1073741824 * (size_t)budget;
            ++i;
          } else {
            error = 1;
            options_error(OPTIONS_EBADCMD);
            break;
          }

          if (options->memory_budget > sinfo.totalram) {
            error = 1;
            options_error(OPTIONS_EBUDGET);
            break;
          }
        } else {
          error = 1;
          options_error(OPTIONS_ENOTCMD);
          break;
        }
      }

      // Statistics output type
      if (!strcmp(argv[i], STATS_OUTPUT_CMD_STRING)) {
        if (i + 1 < argc) {
          if (!strcmp(OUTPUT_TYPE_TEXT_STR, argv[i + 1])) {
            options->stats_output_format = STATS_OUTPUT_FORMAT_TEXT;
          } else if (!strcmp(OUTPUT_TYPE_CSV_STR, argv[i + 1])) {
            options->stats_output_format = STATS_OUTPUT_FORMAT_CSV;
          } else {
            options_error(OPTIONS_EBADCMD);
            break;
          }

          ++i;
        } else {
          options_error(OPTIONS_EBADCMD);
          break;
        }
      }

      // Batch size
      if (!strcmp(argv[i], BATCH_SIZE_CMD_STRING)) {
        if (i + 1 < argc) {
          options->batch_size = atoi(argv[i + 1]);
          ++i;
        }
      }

      // Quality cutoff
      if (!strcmp(argv[i], QUALITY_CUTOFF_CMD_STRING)) {
        if (i + 1 < argc) {
          options->quality_cutoff = atoi(argv[i + 1]);
          ++i;
        }
      }

      // HMC treatment
      if (!strcmp(argv[i], TREATMENT_CMD_STRING)) {
        if (i + 1 < argc) {
          size_t tmp = atoi(argv[i + 1]);

          if (tmp == HMC_MARKER_TAB_SEQ || tmp == HMC_MARKER_OX_SEQ) {
            options->treatment = tmp;
            ++i;
          } else {
            error = 1;
            options_error(OPTIONS_EBADCMD);
            break;
          }
        } else {
          error = 1;
          options_error(OPTIONS_ENOTCMD);
          break;
        }
      }

      // Methylation file
      if (!strcmp(argv[i], MC_FILE_CMD_STRING) || !strcmp(argv[i], MC_FILE_CMD_STRING_SHORT)) {
        if (i + 1 < argc) {
          if (is_file(argv[i + 1])) {
            strcpy(options->mc_bam_file, argv[i + 1]);
            ++i;
          } else {
            error = 1;
            options_error(OPTIONS_ENOFILE);
            break;
          }
        } else {
          error = 1;
          options_error(OPTIONS_ENOTCMD);
          break;
        }
      }

      // Hydroximethylation file
      if (!strcmp(argv[i], HMC_FILE_CMD_STRING) || !strcmp(argv[i], HMC_FILE_CMD_STRING_SHORT)) {
        if (i + 1 < argc) {
          if (is_file(argv[i + 1])) {
            strcpy(options->hmc_bam_file, argv[i + 1]);
            ++i;
          } else {
            error = 1;
            options_error(OPTIONS_ENOFILE);
            break;
          }
        } else {
          error = 1;
          options_error(OPTIONS_ENOTCMD);
          break;
        }
      }

      // Output directory
      if (!strcmp(argv[i], OUTPUT_PATH_CMD_STRING) || !strcmp(argv[i], OUTPUT_PATH_CMD_STRING_SHORT)) {
        if (i + 1 < argc) {
          // Check if the target output path is an existing
          // directory. If it isn't, create it
          struct stat st;
          stat(argv[i + 1], &st);

          if (!S_ISDIR(st.st_mode)) {
            create_directory(argv[i + 1]);
          }

          strcpy(options->output_directory, argv[i + 1]);
          ++i;
        } else {
          error = 1;
          options_error(OPTIONS_ENOTCMD);
          break;
        }
      }

      // Genome index directory
      if (!strcmp(argv[i], INDEX_CMD_STRING) || !strcmp(argv[i], INDEX_CMD_STRING_SHORT)) {
        if (i + 1 < argc) {
          // Check if the target output path is an existing
          // directory
          struct stat st;
          stat(argv[i + 1], &st);

          if (S_ISDIR(st.st_mode)) {
            strcpy(options->index_directory, argv[i + 1]);
            ++i;
          } else {
            error = 1;
            options_error(OPTIONS_ENOTDIR);
            break;
          }
        } else {
          error = 1;
          options_error(OPTIONS_ENOTCMD);
          break;
        }
      }

      // Methylation stats file
      if (!strcmp(argv[i], MC_METH_STATS_CMD_STRING)) {
        if (i + 1 < argc) {
          if (is_file(argv[i + 1])) {
            strcpy(options->mc_methyl_stats_file, argv[i + 1]);
            ++i;
          } else {
            error = 1;
            options_error(OPTIONS_ENOFILE);
            break;
          }
        } else {
          error = 1;
          options_error(OPTIONS_ENOTCMD);
          break;
        }
      }

       // Hydroximethylation stats file
      if (!strcmp(argv[i], HMC_METH_STATS_CMD_STRING)) {
        if (i + 1 < argc) {
          if (is_file(argv[i + 1])) {
            strcpy(options->hmc_methyl_stats_file, argv[i + 1]);
            ++i;
          } else {
            error = 1;
            options_error(OPTIONS_ENOFILE);
            break;
          }
        } else {
          error = 1;
          options_error(OPTIONS_ENOTCMD);
          break;
        }
      }

      // CSV delimiter character
      if (!strcmp(argv[i], CSV_DELIMITER_CMD_STRING)) {
        if (i + 1 < argc) {
          if (strlen(argv[i+1]) == 1) {
            options->csv_delimiter = argv[i+1][0];
          } else {
            if (!strcmp(argv[i+1], "\\t")) {
              options->csv_delimiter = '\t';
            }
          }
        } else {
          error = 1;
          options_error(OPTIONS_ENOTCMD);
          break;
        }
      }

      // CSV record delimiter character
      if (!strcmp(argv[i], CSV_RECORD_DELIMITER_CMD_STRING)) {
        if (i + 1 < argc) {
          if (strlen(argv[i+1]) == 1) {
            options->csv_record_delimiter = argv[i+1][0];
          } else {
            if (!strcmp(argv[i+1], "\\n")) {
              options->csv_record_delimiter = '\n';
            }
          }
        } else {
          error = 1;
          options_error(OPTIONS_ENOTCMD);
          break;
        }
      }

      // Number of threads to use
      if (!strcmp(argv[i], NUM_THREADS_CMD_STRING)) {
        if (i + 1 < argc) {
          size_t tmp = atoi(argv[i + 1]);

          // Check the thread count. It must be bigger than three and ideally
          // lower than the number of CPU's available in the system
          if (tmp > 0) {
            if (tmp > options->num_threads) {
              printf("Warning! Assigning more threads than available processors!\n");
            }

            options->num_threads = tmp;
            ++i;
          } else {
            error = 1;
            options_error(OPTIONS_EBADCMD);
            break;
          }
        } else {
          error = 1;
          options_error(OPTIONS_ENOTCMD);
          break;
        }
      }
    }
  }

  // Check if the mandatory options have been set
  if (!strlen(options->mc_bam_file) || !strlen(options->hmc_bam_file)) {
    error = 1;
    options_error(OPTIONS_ETOOFEW);
  }

  // Print the command line options if the input
  // was badly formed
  if (error) {
    options_print_usage();
  }

  return error;
}
