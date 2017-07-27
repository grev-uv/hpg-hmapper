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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define SEQUENCE  1
#define SEQUENCE_LEN 16384
#define tam SEQUENCE_LEN
#define MAX_ISLES 80000

#define CPG_CHR     0
#define CPG_START   1
#define CPG_END     2
#define CPG_ELEM    3

void replace_seqs(char * input_name, char * zone_name, float limit_reg, float limit_zone, int tissue_count, int samples_per_tissue) {
  srand(time(NULL));

  size_t out_name_length = strlen(input_name) + 50;
  char *output_name = (char *)malloc(out_name_length * sizeof(char));
  snprintf(output_name, out_name_length, "%s_convert.fastq", input_name);

  FILE *input  = fopen(input_name, "r");
  FILE *output = fopen(output_name, "w");
  FILE *zone   = fopen(zone_name, "r");

  size_t count_C  = 0, count_G  = 0;
  size_t count_C_isle  = 0, count_G_isle  = 0;
  size_t count_CT = 0, count_GA = 0;
  size_t count_CT_isle = 0, count_GA_isle = 0;

  const float nr = 1.0f / (float)RAND_MAX;

  char sequence[SEQUENCE_LEN];
  char caracter11, caracter12;
  char caracter22, caracter21;
  char src, dst;

  size_t cpg_zones[MAX_ISLES][CPG_ELEM];
  size_t num_zones = 0;

  float inv_limit_reg = 1.0f - limit_reg;
  size_t chromo, position, strand;
  float limit = limit_reg;
  char *p;
  size_t current_sequence = 0;
  size_t total_sequences  = 0;
  size_t current_tissue   = 0;

  // Generate N random methylation percentages for the CpG islands
  float* tissue_percentages = calloc(tissue_count, sizeof(float));
  size_t* tissue_samples    = calloc(tissue_count, sizeof(size_t));

  for (size_t i = 0; i < tissue_count; ++i) {
    tissue_percentages[i] = rand() * nr;
    if (tissue_percentages[i] > limit_zone) tissue_percentages[i] > limit_zone;
  }

  // Read the CpG islands from the CSV file
  while (!feof(zone)) {
    fscanf(zone, "%lu", &cpg_zones[num_zones][CPG_CHR]);
    fscanf(zone, "%lu", &cpg_zones[num_zones][CPG_START]);
    fscanf(zone, "%lu", &cpg_zones[num_zones][CPG_END]);
    num_zones++;
  }

  num_zones--;

  // For every sequence in the FASTQ
  while(!feof(input)) {
    fgets(sequence, SEQUENCE_LEN, input);

    if (sequence[0] != '@') {
      printf("The FASTQ sequence %lu is malformed. Terminating or EOF reached.\n", total_sequences);
      break;
    }

    // Write the sequence name to the output FASTQ
    fputs(sequence, output);

    // Extract chromosome, strand and start position
    sequence[0] = ' ';
    p = strtok(sequence, "_");
    chromo = atoi(p);
    p = strtok(NULL, "_");
    position = atoi(p);
    p = strtok(NULL, "_");
    p = strtok(NULL, "_");
    strand = atoi(p);

    // Simulate random strand amplification mutation 
    fgets(sequence, SEQUENCE_LEN, input);

    if (rand() % 2 == 0) {
      src = 'C';
      dst = 'T';
    } else {
      src = 'G';
      dst = 'A';
    }

    // Detect if the read is inside a CpG island or not
    limit = inv_limit_reg;

    for (size_t i = 0; i < num_zones; i++) {
      if (cpg_zones[i][CPG_CHR] == chromo) {
	if ((position >= cpg_zones[i][CPG_START] && position <= cpg_zones[i][CPG_END]) ||
	    (position + strlen(sequence) >= cpg_zones[i][CPG_START] && position + strlen(sequence) <= cpg_zones[i][CPG_END])) {
	  limit = 1.0f - tissue_percentages[current_tissue];
	  break;
	}
      }
    }

    // Handle tissue switching
    ++current_sequence;
    ++total_sequences;

    if (current_sequence > samples_per_tissue) {
      if (current_tissue == tissue_count - 1) {
        current_tissue = 0;
      } else {

        ++current_tissue;
      }

      current_sequence = 0;
    }

    // Transform the read
    for (size_t i = 0; i < strlen(sequence); i++) {
      if (sequence[i] == src) {
        if (limit == inv_limit_reg) {
          if (src == 'C') ++count_C;
          if (src == 'G') ++count_G;
        } else {
          if (src == 'C') ++count_C_isle;
          if (src == 'G') ++count_G_isle;
        }

	if ((float)rand() * nr > limit) {
	  sequence[i] = dst;

          if (limit == inv_limit_reg) {
            if (dst == 'T') ++count_CT;
	    if (dst == 'A') ++count_GA;
          } else {
            if (dst == 'T') ++count_CT_isle;
            if (dst == 'A') ++count_GA_isle;
            ++tissue_samples[current_tissue];
          }
	}
      }
    }

    // Write the sequence string to the output FASTQ
    fputs(sequence, output);

    // Write the strand to the output FASTQ
    fgets(sequence, SEQUENCE_LEN, input);
    fputs(sequence, output);

    // Write the quality string to the output FASTQ
    fgets(sequence, SEQUENCE_LEN, input);
    fputs(sequence, output);
  }

  fclose(input);
  fclose(output);
  fclose(zone);

  printf("Processed sequences: %lu\n", total_sequences);

  printf("Number of C: %lu, methylated :%lu\n", count_C, count_CT);
  printf("Number of G: %lu, methylated :%lu\n", count_G, count_GA);
  printf("Number of C in CpG: %lu, methylated :%lu\n", count_C_isle, count_CT_isle);
  printf("Number of G in CpG: %lu, methylated :%lu\n", count_G_isle, count_GA_isle);

  printf("Simulating %d tissues in CpG isles, rates and methylated CG:\n", tissue_count);

  size_t tissue_total = 0;

  for (size_t i = 0; i < tissue_count; ++i) {
    printf("[Tissue %lu] CpG meth rate: %.02f%%, methylated CG: %lu\n", i + 1, tissue_percentages[i] * 100.0f, tissue_samples[i]);
    tissue_total += tissue_samples[i];
  }

  printf("[All tissues] %lu nt\n", tissue_total);

  free(tissue_percentages);
  free(tissue_samples);
}


int main(int argc, char *argv[]){
  if (argc < 2){
    printf("error:\n");
    printf("%s option input_file [...]\n", argv[0]);
    printf("option:\n%d\tSequence Bisulfite treatment\n", SEQUENCE);
    return 0;
  }

  int opt = atoi(argv[1]);
  
  switch(opt) {
  case SEQUENCE:
    if (argc < 8){
      printf("error:\n");
      printf("%s option [fastq_input_file] [cpg_zone_file] [methilated_percent] [max_methilated_percent_in_cpg_island] [tissue_count] [samples_per_tissue]\n", argv[0]);
      return 0;
    }

    printf("Replace sequence nucleotides with Bisulphite\n");
    replace_seqs(argv[2], argv[3], 0.01f * atof(argv[4]), 0.01f * atof(argv[5]), atoi(argv[6]), atoi(argv[7]));
    break;
  default:
    printf("Incorrect value to replace nucleotides\n");
  }
  
  return 0;
}