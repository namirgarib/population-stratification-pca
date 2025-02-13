/**
 * @file variant_calling.c
 * @brief Variant calling implementation with sparse representation
 *
 * Author: Namir Garib
 * Created: January 2025
 */

#include "variant_calling.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

static double _calculate_substitution_score(char ref, char alt, const VariantParams *params) {
  /* transitions = { A<->G, C<->T } */
  const char *transitions[] = {"AG", "GA", "CT", "TC", NULL};
  for (const char **t = transitions; *t; t++) {
    if (ref == (*t)[0] && alt == (*t)[1]) {
      return params->transition_weight;
    }
  }
  return params->transversion_weight;
}

static bool _is_cpg_site(const char *ref_chunk, size_t pos, size_t chunk_len) {
  if (pos > 0 && pos < chunk_len - 1) {
    return (ref_chunk[pos] == 'C' && ref_chunk[pos + 1] == 'G') ||
           (ref_chunk[pos] == 'G' && ref_chunk[pos - 1] == 'C');
  }
  return false;
}

static size_t _count_context_variants(const char *ref_chunk, const char *sample_chunk, size_t pos,
                                      size_t chunk_len) {
  const size_t radius = 2;
  size_t       start = (pos > radius) ? pos - radius : 0;
  size_t       end = (pos + radius < chunk_len) ? pos + radius : (chunk_len - 1);
  size_t       count = 0;
  for (size_t i = start; i <= end; i++) {
    if (ref_chunk[i] != sample_chunk[i]) {
      count++;
    }
  }
  return count;
}

/* -------------- Sparse Structures -------------- */

void init_individual_variants(IndividualVariants *iv) {
  iv->entries = NULL;
  iv->num_entries = 0;
  iv->capacity = 0;
}

void free_individual_variants(IndividualVariants *iv) {
  if (iv->entries) {
    free(iv->entries);
  }
  iv->entries = NULL;
  iv->num_entries = 0;
  iv->capacity = 0;
}

static void add_variant(IndividualVariants *iv, size_t col, double score) {
  if (iv->num_entries == iv->capacity) {
    size_t        new_cap = (iv->capacity == 0) ? 1024 : (iv->capacity * 2);
    VariantEntry *tmp = realloc(iv->entries, new_cap * sizeof(VariantEntry));
    if (!tmp) {
      fprintf(stderr, "add_variant: Out of memory.\n");
      exit(EXIT_FAILURE);
    }
    iv->entries = tmp;
    iv->capacity = new_cap;
  }
  iv->entries[iv->num_entries].col = col;
  iv->entries[iv->num_entries].score = score;
  iv->num_entries++;
}

/* ------------------- Public call_variants_chunk() ------------------- */

void call_variants_chunk(const char *ref_chunk, const char *sample_chunk, size_t chunk_len,
                         size_t global_offset, IndividualVariants *ivar,
                         const VariantParams *params) {
  if (!ref_chunk || !sample_chunk || !ivar || !params) return;

  for (size_t i = 0; i < chunk_len; i++) {
    if (ref_chunk[i] == sample_chunk[i]) {
      /* no variant => skip */
      continue;
    }
    char ref_base = ref_chunk[i];
    char alt_base = sample_chunk[i];
    /* Basic validity check */
    if (ref_base != 'A' && ref_base != 'C' && ref_base != 'G' && ref_base != 'T') {
      /* skip if non-standard base in ref or handle 'N' etc. */
      continue;
    }
    if (alt_base != 'A' && alt_base != 'C' && alt_base != 'G' && alt_base != 'T') {
      continue;
    }

    double score = _calculate_substitution_score(ref_base, alt_base, params);

    if (_is_cpg_site(ref_chunk, i, chunk_len)) {
      score *= params->cpg_multiplier;
    }
    size_t ctx_vars = _count_context_variants(ref_chunk, sample_chunk, i, chunk_len);
    score *= 1.0 + (params->cluster_factor * ctx_vars);

    /* logistic transform => [0,1] range */
    double logistic_score = 1.0 / (1.0 + exp(-params->logistic_scale * score));
    if (logistic_score > 0.0) {
      add_variant(ivar, global_offset + i, logistic_score);
    }
  }
}
