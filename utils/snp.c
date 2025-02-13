/**
 * @file snp.c
 * @brief SNP Generator for Simulated Individuals with Large Genome Support
 *
 * This program reads a genome sequence from a reference file in chunks and generates 10 individual
 * genome sequences with artificially introduced Single Nucleotide Polymorphisms (SNPs).
 * Each individual has a random number of unique SNPs and a random number of group-wide SNPs.
 * The group-wide SNPs are shared by a subset of individuals, with different probabilities for each
 * group. The reference genome is read in chunks to support large genome sequences that may not fit
 * in stack memory.
 *
 * @author Namir Garib
 * @date 2025-01-18
 * @version 1.1
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NUM_INDIVIDUALS 10
#define GROUP_1_SIZE 4      // individuals 0..3
#define GROUP_2_SIZE 4      // individuals 4..7
#define GROUP_3_SIZE 2      // individuals 8..9
#define CHUNK_SIZE 1000000  // Process 1 million bases at a time

#define GROUP_SHARED_PROB 0.3  // Probability that a chosen SNP is group-wide (instead of unique)
#define GROUP_1_PROB 0.4       // Among group-wide events, 40% for group 1
#define GROUP_2_PROB 0.4       // 40% for group 2
#define GROUP_3_PROB 0.2       // 20% for group 3

/**
 * @brief Randomly pick a different nucleotide from the reference sequence.
 */
char random_base_different_from(char original) {
  char bases[] = "ACGT";
  char new_base;
  do {
    new_base = bases[rand() % 4];
  } while (new_base == original);
  return new_base;
}

/**
 * @brief Read the reference genome into a streaming buffer, introducing SNPs for groups
 * and individuals.
 *
 * @param reference_file Path to the reference genome file.
 * @param min_snps Minimum number of SNPs to introduce per chunk.
 * @param max_snps Maximum number of SNPs to introduce per chunk.
 */
void process_large_genome(const char *reference_file, int min_snps, int max_snps) {
  FILE *ref_file = fopen(reference_file, "r");
  if (!ref_file) {
    fprintf(stderr, "Error opening reference file: %s\n", reference_file);
    exit(EXIT_FAILURE);
  }

  // Create output files for each individual
  FILE *out_files[NUM_INDIVIDUALS];
  char  filename[32];
  for (int i = 0; i < NUM_INDIVIDUALS; i++) {
    sprintf(filename, "ind%d.txt", i + 1);
    out_files[i] = fopen(filename, "w");
    if (!out_files[i]) {
      fprintf(stderr, "Error creating file: %s\n", filename);
      exit(EXIT_FAILURE);
    }
  }

  srand((unsigned)time(NULL));

  // Allocate buffers on the HEAP instead of the stack.
  char *ref_chunk = (char *)malloc(CHUNK_SIZE * sizeof(char));
  if (!ref_chunk) {
    fprintf(stderr, "Failed to allocate ref_chunk\n");
    exit(EXIT_FAILURE);
  }

  char **ind_chunk = (char **)malloc(NUM_INDIVIDUALS * sizeof(char *));
  if (!ind_chunk) {
    fprintf(stderr, "Failed to allocate ind_chunk pointers\n");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < NUM_INDIVIDUALS; i++) {
    ind_chunk[i] = (char *)malloc(CHUNK_SIZE * sizeof(char));
    if (!ind_chunk[i]) {
      fprintf(stderr, "Failed to allocate ind_chunk[%d]\n", i);
      exit(EXIT_FAILURE);
    }
  }

  // We'll keep track of the global position so we can print where SNPs occur.
  size_t global_position = 0;
  size_t total_length = 0;
  size_t read_size;

  // Read reference in a streaming fashion
  while ((read_size = fread(ref_chunk, sizeof(char), CHUNK_SIZE, ref_file)) > 0) {
    // Copy reference chunk into each individual's chunk
    for (int ind = 0; ind < NUM_INDIVIDUALS; ind++) {
      memcpy(ind_chunk[ind], ref_chunk, read_size);
    }

    // Decide how many SNPs we want to introduce in this chunk
    int snps_this_chunk = min_snps + rand() % (max_snps - min_snps + 1);

    // Now let's choose random positions in this chunk to mutate
    for (int s = 0; s < snps_this_chunk; s++) {
      size_t pos_in_chunk = rand() % read_size;
      char   original_base = ref_chunk[pos_in_chunk];
      if (original_base != 'A' && original_base != 'C' && original_base != 'G' &&
          original_base != 'T') {
        // If reference isn't a standard base, skip it
        continue;
      }

      // Decide if this SNP is group-wide or unique
      double r = (double)rand() / RAND_MAX;
      char   new_base = random_base_different_from(original_base);

      if (r < GROUP_SHARED_PROB) {
        // It's a group-wide SNP; decide which group
        double rg = (double)rand() / RAND_MAX;
        if (rg < GROUP_1_PROB) {
          // Group 1
          for (int i = 0; i < GROUP_1_SIZE; i++) {
            ind_chunk[i][pos_in_chunk] = new_base;
          }
        } else if (rg < GROUP_1_PROB + GROUP_2_PROB) {
          // Group 2
          for (int i = GROUP_1_SIZE; i < GROUP_1_SIZE + GROUP_2_SIZE; i++) {
            ind_chunk[i][pos_in_chunk] = new_base;
          }
        } else {
          // Group 3
          for (int i = GROUP_1_SIZE + GROUP_2_SIZE; i < NUM_INDIVIDUALS; i++) {
            ind_chunk[i][pos_in_chunk] = new_base;
          }
        }
      } else {
        // Unique SNP for exactly one individual
        int chosen_ind = rand() % NUM_INDIVIDUALS;
        ind_chunk[chosen_ind][pos_in_chunk] = new_base;
      }
    }

    // Write out each individual's chunk to its file
    for (int ind = 0; ind < NUM_INDIVIDUALS; ind++) {
      fwrite(ind_chunk[ind], sizeof(char), read_size, out_files[ind]);
    }

    total_length += read_size;
    global_position += read_size;
  }

  // Cleanup
  fclose(ref_file);
  for (int i = 0; i < NUM_INDIVIDUALS; i++) {
    fclose(out_files[i]);
    free(ind_chunk[i]);
  }
  free(ind_chunk);
  free(ref_chunk);

  printf("\nProcessed %zu bases from %s.\n", total_length, reference_file);
  printf("SNP-modified genome sequences generated for %d individuals.\n", NUM_INDIVIDUALS);
}

int main(int argc, char *argv[]) {
  if (argc < 4) {
    fprintf(stderr, "Usage: %s <reference_file> <min_snps> <max_snps>\n", argv[0]);
    return EXIT_FAILURE;
  }

  const char *reference_file = argv[1];
  int         min_snps = atoi(argv[2]);
  int         max_snps = atoi(argv[3]);
  if (min_snps <= 0 || max_snps <= 0 || min_snps > max_snps) {
    fprintf(stderr, "Invalid SNP range\n");
    return EXIT_FAILURE;
  }

  process_large_genome(reference_file, min_snps, max_snps);

  return EXIT_SUCCESS;
}
