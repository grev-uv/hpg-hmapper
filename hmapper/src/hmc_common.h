/**
 * hmc_common.h
 * - Date: 1 / Feb / 2017
 * - Who: Cesar Gonzalez (Cesar.Gonzalez-Segura@uv.es)
 *
 * Common definitions and constants for the hmC pipeline.
 *
 */

#ifndef HMC_COMMON_H
#define HMC_COMMON_H

#include <stdint.h>

/**
 * Hydroximethylation common constants
 */
#define HMC_MARKER_TAB_SEQ    0   /**< Amplified using Tab-SEQ */
#define HMC_MARKER_OX_SEQ     1   /**< Amplified using Ox-SEQ */

#define WF_READ_PRODUCER      0   /**< BAM reader task identifier */
#define WF_MAP_WORKER         1   /**< Mapper task identifier */
#define WF_WRITE_CONSUMER     2   /**< Writter task identifier */

#define THR_UNUSED            -1  /**< Unused slot indicator for the scheduling matrix */

#define PHRED_BASE_QUALITY    33  /**< Base quality in Phred scale using in the alignments */

#define MC_QUEUE_INDEX         0    /**< The pipeline stage is processing a mC batch  */
#define HMC_QUEUE_INDEX        1    /**< The pipeline stage is processing a hmC batch */
#define METH_QUEUE_COUNT       2

#define STATS_OUTPUT_FORMAT_TEXT    0   /**< Output statistics and timing as plain text */
#define STATS_OUTPUT_FORMAT_CSV     1   /**< Output statistics and timing as a CSV */

#define SCHEDULER_MEMORY_EVENT_READ_SIZE    100000     /**< Number of bytes to track before notifying a memory alloc or free
                                                            event to the scheduler memory manager */

/**
 * Default homo sapiens GRCh37.68 genome chromosome data
 */
#define DEFAULT_GENOME_NUM_CHROMOSOMES    24
static const uint64_t DEFAULT_GENOME_CHR_LENGTH[] = {
  249250620,    /**< Chromosome 1 */
  243199372,    /**< Chromosome 2 */
  198022429,    /**< Chromosome 3 */
  191154275,    /**< Chromosome 4 */
  180915259,    /**< Chromosome 5 */
  171115066,    /**< Chromosome 6 */
  159138662,    /**< Chromosome 7 */
  146364021,    /**< Chromosome 8 */
  141213430,    /**< Chromosome 9 */
  135534746,    /**< Chromosome 10 */
  135006515,    /**< Chromosome 11 */
  133851894,    /**< Chromosome 12 */
  115169877,    /**< Chromosome 13 */
  107349539,    /**< Chromosome 14 */
  102531391,    /**< Chromosome 15 */
  90354752,     /**< Chromosome 16 */
  81195209,     /**< Chromosome 17 */
  78077247,     /**< Chromosome 18 */
  59128982,     /**< Chromosome 19 */
  63025519,     /**< Chromosome 20 */
  48129894,     /**< Chromosome 21 */
  51304565,     /**< Chromosome 22 */
  155270559,    /**< Chromosome X */
  59373565      /**< Chromosome Y */
};

#endif /* HMC_COMMON_H */
