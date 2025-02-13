#ifndef VARIANT_CALLING_H
#define VARIANT_CALLING_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double transition_weight;
  double transversion_weight;
  double cpg_multiplier;
  double cluster_factor;
  double logistic_scale;
} VariantParams;

typedef struct {
  size_t col;
  double score;
} VariantEntry;

typedef struct {
  VariantEntry *entries;
  size_t        num_entries;
  size_t        capacity;
} IndividualVariants;

/**
 * @brief Initializes an IndividualVariants struct.
 */
void init_individual_variants(IndividualVariants *iv);

/**
 * @brief Frees all memory in an IndividualVariants struct.
 */
void free_individual_variants(IndividualVariants *iv);

/**
 * @brief Processes a chunk of data (ref_chunk vs. sample_chunk),
 *        storing only nonzero variant scores in a sparse structure.
 *
 * @param ref_chunk     reference chunk
 * @param sample_chunk  individual's chunk
 * @param chunk_len     length of chunk read
 * @param global_offset position in the overall genome
 * @param ivar          where to store sparse entries
 * @param params        variant-calling parameters
 */
void call_variants_chunk(const char *ref_chunk, const char *sample_chunk, size_t chunk_len,
                         size_t global_offset, IndividualVariants *ivar,
                         const VariantParams *params);

#ifdef __cplusplus
}
#endif

#endif
