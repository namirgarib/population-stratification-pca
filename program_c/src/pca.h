#ifndef PCA_H
#define PCA_H

#include <stddef.h>

#include "variant_calling.h"

/**
 * @brief PartialPCA_Result: structure to hold partial PCA output.
 *
 *  - num_components: how many principal components found
 *  - eigenvalues: length = num_components
 *  - pc_vectors: length = num_components * d (row-major).
 *  - scores: shape (n_inds x num_components).
 */
typedef struct {
  int     num_components;
  double *eigenvalues;
  double *pc_vectors;
  double *scores;
} PartialPCA_Result;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief partial_pca_sparse: naive partial PCA with power method
 *        on sparse variant data. Final scores are z-scored
 *        so each PC dimension is mean=0 and stdev=1 across individuals.
 *
 * @param ivars   array of n_inds, each with (col, score)
 * @param n_inds  number of individuals
 * @param d       dimension (# of columns or max base index + 1)
 * @param top_k   # of principal components to extract
 * @return        PartialPCA_Result (caller frees .eigenvalues, .pc_vectors, .scores)
 */
PartialPCA_Result partial_pca_sparse(const IndividualVariants *ivars, int n_inds, size_t d,
                                     int top_k);

#ifdef __cplusplus
}
#endif

#endif
