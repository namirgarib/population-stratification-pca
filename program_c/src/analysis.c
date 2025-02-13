/**
 * @file analysis.c
 * @brief High-level pipeline for variant calling + partial PCA
 *
 * Author: Namir Garib
 * Created: January 2025
 */

#include "analysis.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#include "pca.h"
#include "utils.h"
#include "variant_calling.h"

/* Default variant calling parameters */
VariantParams params = {.transition_weight = 0.28,
                        .transversion_weight = 1.1,
                        .cpg_multiplier = 1.8,
                        .cluster_factor = 0.12,
                        .logistic_scale = 0.6};

/**
 * @brief Return file length in bytes (minus trailing newline if present).
 */
static size_t get_file_length(const char *filepath) {
  FILE *f = fopen(filepath, "rb");
  if (!f) {
    fprintf(stderr, "get_file_length: Cannot open %s\n", filepath);
    return 0;
  }
  fseeko(f, 0, SEEK_END);
  off_t sz = ftello(f);
  if (sz < 0) {
    fprintf(stderr, "get_file_length: ftello error on %s\n", filepath);
    fclose(f);
    return 0;
  }
  /* if last char is newline => ignore */
  if (sz > 0) {
    fseeko(f, -1, SEEK_END);
    char last_ch;
    if (fread(&last_ch, 1, 1, f) == 1) {
      if (last_ch == '\n') {
        sz -= 1;
      }
    }
  }
  fclose(f);
  return (size_t)sz;
}

/**
 * @brief We'll read ref_file + each individual's file in chunk_size blocks,
 *        call call_variants_chunk(...) for each block,
 *        and store results in a sparse array for each individual.
 */
static int gather_variants_sparse(const char *ref_file, const char **individuals, int n_inds,
                                  size_t chunk_size, IndividualVariants *ivars_out,
                                  size_t *total_len_out) {
  FILE *fref = fopen(ref_file, "rb");
  if (!fref) {
    fprintf(stderr, "Failed to open reference %s\n", ref_file);
    return 1;
  }

  FILE **f_ind = secure_alloc(n_inds * sizeof(FILE *));
  for (int i = 0; i < n_inds; i++) {
    f_ind[i] = fopen(individuals[i], "rb");
    if (!f_ind[i]) {
      fprintf(stderr, "Failed to open individual %s\n", individuals[i]);
      for (int j = 0; j < i; j++) {
        fclose(f_ind[j]);
      }
      free(f_ind);
      fclose(fref);
      return 1;
    }
  }

  char  *ref_chunk = secure_alloc(chunk_size);
  char **ind_chunk = secure_alloc(n_inds * sizeof(char *));
  for (int i = 0; i < n_inds; i++) {
    ind_chunk[i] = secure_alloc(chunk_size);
  }

  size_t global_offset = 0;
  while (!feof(fref)) {
    size_t nref = fread(ref_chunk, 1, chunk_size, fref);
    if (nref == 0) break; /* done */

    /* read same amount from each individual's file */
    for (int i = 0; i < n_inds; i++) {
      size_t nind = fread(ind_chunk[i], 1, nref, f_ind[i]);
      if (nind < nref) {
        /* mismatch or partial read => possible error */
        fprintf(stderr, "Warning: partial read for %s\n", individuals[i]);
      }
    }

    /* call variants for each individual */
    for (int i = 0; i < n_inds; i++) {
      call_variants_chunk(ref_chunk, ind_chunk[i], nref, global_offset, &ivars_out[i], &params);
    }
    global_offset += nref;

    if (nref < chunk_size) {
      break; /* near EOF */
    }
  }

  *total_len_out = global_offset;

  /* cleanup */
  fclose(fref);
  for (int i = 0; i < n_inds; i++) {
    fclose(f_ind[i]);
    free(ind_chunk[i]);
  }
  free(ind_chunk);
  free(f_ind);
  free(ref_chunk);

  return 0;
}

static const char *shorten_bases(size_t bases) {
  static char buf[32];
  if (bases >= 1000000) {
    snprintf(buf, sizeof(buf), "%zuM", bases / 1000000);
  } else if (bases >= 1000) {
    snprintf(buf, sizeof(buf), "%zuk", bases / 1000);
  } else {
    snprintf(buf, sizeof(buf), "%zu", bases);
  }
  return buf;
}

void perform_full_analysis(const char *ref_file, int num_individuals,
                           const char **individuals_files) {
  if (!ref_file || num_individuals <= 0 || !individuals_files) {
    fprintf(stderr, "perform_full_analysis: invalid args.\n");
    return;
  }

  /* Initialize sparse arrays for each individual */
  IndividualVariants *all_ivars = secure_alloc(num_individuals * sizeof(IndividualVariants));
  for (int i = 0; i < num_individuals; i++) {
    init_individual_variants(&all_ivars[i]);
  }

  /* Stream-based reading in chunks of 1 million bytes */
  size_t chunk_size = 1000000;
  size_t total_len = 0;
  int    err = gather_variants_sparse(ref_file, individuals_files, num_individuals, chunk_size,
                                      all_ivars, &total_len);
  if (err) {
    fprintf(stderr, "perform_full_analysis: gather_variants_sparse failed.\n");
    for (int i = 0; i < num_individuals; i++) {
      free_individual_variants(&all_ivars[i]);
    }
    free(all_ivars);
    return;
  }

  printf("Reference genome size (streamed) = %zu bases\n", total_len);

  /* 2) Perform partial PCA => we ask for top_k=4 or however many you want. */
  int               top_k = 4;
  PartialPCA_Result pca_res = partial_pca_sparse(all_ivars, num_individuals, total_len, top_k);

  /* 3) Build an output folder for the results */
  mkdir("./results", 0777);

  time_t     now = time(NULL);
  struct tm *t_info = localtime(&now);
  char       time_buf[32];
  strftime(time_buf, sizeof(time_buf), "%Y%m%d%H%M%S", t_info);

  char out_folder[128];
  snprintf(out_folder, sizeof(out_folder), "./results/%s_%s", time_buf, shorten_bases(total_len));
  mkdir(out_folder, 0777);

  /* 4) Write scores => shape (num_individuals x top_k) */
  char results_file[256];
  snprintf(results_file, sizeof(results_file), "%s/results.csv", out_folder);
  FILE *f_scores = fopen(results_file, "w");
  if (!f_scores) {
    fprintf(stderr, "Cannot open %s for writing.\n", results_file);
  } else {
    for (int i = 0; i < num_individuals; i++) {
      for (int comp = 0; comp < top_k; comp++) {
        fprintf(f_scores, "%.6f", pca_res.scores[i * top_k + comp]);
        if (comp < top_k - 1) fprintf(f_scores, ",");
      }
      fprintf(f_scores, "\n");
    }
    fclose(f_scores);
  }

  /* 5) Write eigenvalues => shape top_k */
  char eval_file[256];
  snprintf(eval_file, sizeof(eval_file), "%s/eigenvalues.csv", out_folder);
  FILE *f_eval = fopen(eval_file, "w");
  if (f_eval) {
    for (int k = 0; k < top_k; k++) {
      fprintf(f_eval, "%d,%.6f\n", k + 1, pca_res.eigenvalues[k]);
    }
    fclose(f_eval);
  }

  printf("Results written to %s/\n", out_folder);

  /* cleanup */
  for (int i = 0; i < num_individuals; i++) {
    free_individual_variants(&all_ivars[i]);
  }
  free(all_ivars);

  free(pca_res.eigenvalues);
  free(pca_res.pc_vectors);
  free(pca_res.scores);
}
