/**
 * @file main.c
 * @brief Main entry point for PCA-based population stratification
 *
 * Author: Namir Garib
 * Created: January 2025
 */

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "analysis.h"
#include "utils.h"

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <ref_genome> <ind_folder>\n", argv[0]);
    return 1;
  }

  const char *ref_file = argv[1];
  const char *ind_folder = argv[2];

  /* read individuals from the folder (all .txt files) */
  DIR *dir = opendir(ind_folder);
  if (!dir) {
    fprintf(stderr, "Cannot open directory: %s\n", ind_folder);
    return 1;
  }

  struct dirent *entry;
  int            num_individuals = 0;
  while ((entry = readdir(dir)) != NULL) {
    if (entry->d_type == DT_REG) {
      size_t len = strlen(entry->d_name);
      if (len > 4 && strcmp(entry->d_name + (len - 4), ".txt") == 0) {
        num_individuals++;
      }
    }
  }
  rewinddir(dir);

  if (num_individuals == 0) {
    fprintf(stderr, "No .txt files found in %s\n", ind_folder);
    closedir(dir);
    return 1;
  }

  const char **individuals_files = secure_alloc(num_individuals * sizeof(char *));
  int          idx = 0;
  while ((entry = readdir(dir)) != NULL && idx < num_individuals) {
    if (entry->d_type == DT_REG) {
      size_t len = strlen(entry->d_name);
      if (len > 4 && strcmp(entry->d_name + (len - 4), ".txt") == 0) {
        size_t path_len = strlen(ind_folder) + 1 + strlen(entry->d_name) + 1;
        char  *fullp = secure_alloc(path_len);
        snprintf(fullp, path_len, "%s/%s", ind_folder, entry->d_name);
        individuals_files[idx++] = fullp;
      }
    }
  }
  closedir(dir);

  /* run analysis */
  perform_full_analysis(ref_file, num_individuals, individuals_files);

  /* free allocated file paths */
  for (int i = 0; i < num_individuals; i++) {
    free((void *)individuals_files[i]);
  }
  free(individuals_files);

  return 0;
}
