// This script is part of HPG-Hmapper
//
// HPG-Hmapper is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// HPG-Hmapper is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with HPG-Hmapper. If not, see <http://www.gnu.org/licenses/>.

// Creator: César González
// Date: March 2017

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define SEQ_LENGTH  16384

int main (int argc, char* argv[]) {
  if (argc < 3) {
    printf("Usage: dwgsid_fastq_hmc_convert [hmc_fastq_file] [hmc_ratio]\n");
    exit(1);
  }

  double hmc_ratio = atof(argv[2]) / 100.0;
  double non_hmc_ratio = 1.0 - hmc_ratio;

  char out_filename[FILENAME_MAX];
  snprintf(out_filename, FILENAME_MAX, "%s_hmc_convert.fastq", argv[1]);

  char id[SEQ_LENGTH];
  char seq[SEQ_LENGTH];
  char strand[SEQ_LENGTH];
  char quality[SEQ_LENGTH];
  char *result = NULL;

  int seq_length = 0;

  size_t converted_c = 0;
  size_t total_c = 0;

  double inv_rm = 1.0 / (double)RAND_MAX;

  FILE* fastq_in = fopen(argv[1], "r");
  FILE* fastq_out = fopen(out_filename, "w");

  if (fastq_in && fastq_out) {
    while (!feof(fastq_in)) {
      result = fgets(id, SEQ_LENGTH, fastq_in);
      result = fgets(seq, SEQ_LENGTH, fastq_in);
      result = fgets(strand, SEQ_LENGTH, fastq_in);
      result = fgets(quality, SEQ_LENGTH, fastq_in);

      if (id && seq && strand && quality) {
        if (seq_length == 0) seq_length = strlen(seq);

        for (int i = 0; i < seq_length; ++i) {
          if (seq[i] == 'C') {
            double r = (double)rand() * inv_rm;
            ++total_c;

            if (r > non_hmc_ratio) {
              seq[i] = 'T';
              ++converted_c;
            }
          }
        }

        fputs(id, fastq_out);
        fputs(seq, fastq_out);
        fputs(strand, fastq_out);
        fputs(quality, fastq_out);
      } else {
        break;
      }
    }
  } else {
    printf("Error opening FASTQ files\n");
  }

  printf("Converted %lu of total %lu (%.03f%%) C's to T's.\n", converted_c, total_c, (float)converted_c / (float)total_c * 100.0f);

  return 0;
}