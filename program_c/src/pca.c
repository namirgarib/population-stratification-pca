/**
 * @file pca.c
 * @brief Sparse partial PCA with final z-score normalization of PC scores
 *
 * Author: Namir Garib
 * Created: January 2025
 */

#include "pca.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

/* ---------------------------------- *
 *   multiply_X_vec: y = X * v        *
 * ---------------------------------- *
 * X is n_inds x d in sparse form,
 * y is length n_inds. For row i:
 *   y[i] = sum_{(col,score)}(score * v[col]).
 */
static void multiply_X_vec(const IndividualVariants *ivars, int n_inds, const double *v,
                           double *y) {
  /* zero out y */
  for (int i = 0; i < n_inds; i++) {
    y[i] = 0.0;
  }
  /* accumulate */
  for (int i = 0; i < n_inds; i++) {
    const IndividualVariants *indv = &ivars[i];
    for (size_t e = 0; e < indv->num_entries; e++) {
      size_t c = indv->entries[e].col;
      double sc = indv->entries[e].score;
      y[i] += sc * v[c];
    }
  }
}

/* ----------------------------------- *
 *   multiply_Xt_vec: z = X^T * y      *
 * ----------------------------------- *
 * X is n_inds x d. z is length d.
 */
static void multiply_Xt_vec(const IndividualVariants *ivars, int n_inds, size_t d, const double *y,
                            double *z) {
  memset(z, 0, d * sizeof(double));
  for (int i = 0; i < n_inds; i++) {
    double                    y_i = y[i];
    const IndividualVariants *indv = &ivars[i];
    for (size_t e = 0; e < indv->num_entries; e++) {
      size_t c = indv->entries[e].col;
      double sc = indv->entries[e].score;
      z[c] += sc * y_i;
    }
  }
}

/* ------------------------------------------------------- *
 *   remove_component: v <- v - (dot(v, pc)*pc)            *
 * ------------------------------------------------------- */
static void remove_component(double *vec, const double *pc_vec, size_t length) {
  double dot = 0.0;
  for (size_t i = 0; i < length; i++) {
    dot += vec[i] * pc_vec[i];
  }
  for (size_t i = 0; i < length; i++) {
    vec[i] -= dot * pc_vec[i];
  }
}

/* -------------------------------------------- *
 *  z_score_pca_scores: for each PC dimension,  *
 *  subtract the mean and scale to unit stdev.  *
 * -------------------------------------------- */
static void z_score_pca_scores(double *scores, int n_inds, int top_k) {
  /* For comp=0..top_k-1: compute mean & stdev of that column. */
  for (int comp_i = 0; comp_i < top_k; comp_i++) {
    /* 1) mean */
    double sum = 0.0;
    for (int i = 0; i < n_inds; i++) {
      sum += scores[i * top_k + comp_i];
    }
    double mean = sum / (double)n_inds;

    /* subtract mean */
    for (int i = 0; i < n_inds; i++) {
      scores[i * top_k + comp_i] -= mean;
    }

    /* 2) stdev */
    double sum_sq = 0.0;
    for (int i = 0; i < n_inds; i++) {
      double val = scores[i * top_k + comp_i];
      sum_sq += val * val;
    }
    double variance = (n_inds > 1) ? (sum_sq / (double)(n_inds - 1)) : 0.0;
    double stddev = (variance > 0.0) ? sqrt(variance) : 1.0;

    /* divide by stddev */
    for (int i = 0; i < n_inds; i++) {
      scores[i * top_k + comp_i] /= stddev;
    }
  }
}

/* --------------------------------------------------- *
 *  partial_pca_sparse                                 *
 *  1) naive power method for top_k PCs in sparse data  *
 *  2) final z-score of PC scores to produce a typical  *
 *     PCA-like scale across individuals               *
 * --------------------------------------------------- */
PartialPCA_Result partial_pca_sparse(const IndividualVariants *ivars, int n_inds, size_t d,
                                     int top_k) {
  PartialPCA_Result result;
  result.num_components = top_k;

  result.eigenvalues = secure_alloc(top_k * sizeof(double));
  result.pc_vectors = secure_alloc((size_t)top_k * d * sizeof(double));
  result.scores = secure_alloc((size_t)n_inds * top_k * sizeof(double));

  double *v = secure_alloc(d * sizeof(double));
  double *z = secure_alloc(d * sizeof(double));
  double *y = secure_alloc(n_inds * sizeof(double));

  /* For comp_i in [0..top_k-1], do a naive power iteration. */
  for (int comp_i = 0; comp_i < top_k; comp_i++) {
    /* 1) Initialize v randomly. */
    for (size_t c = 0; c < d; c++) {
      v[c] = (double)rand() / RAND_MAX - 0.5;
    }

    /* remove subspace of previously found comps */
    for (int prev = 0; prev < comp_i; prev++) {
      remove_component(v, &result.pc_vectors[prev * d], d);
    }

    /* normalize v */
    {
      double norm = 0.0;
      for (size_t c = 0; c < d; c++) {
        norm += v[c] * v[c];
      }
      norm = sqrt(norm);
      if (norm < 1e-15) norm = 1.0;
      for (size_t c = 0; c < d; c++) {
        v[c] /= norm;
      }
    }

    /* 2) do ~20 power iterations */
    const int MAX_ITERS = 20;
    for (int iter = 0; iter < MAX_ITERS; iter++) {
      /* y = X*v => length n_inds */
      multiply_X_vec(ivars, n_inds, v, y);

      /* z = X^T*y => length d */
      multiply_Xt_vec(ivars, n_inds, d, y, z);

      /* remove subspace again */
      for (int prev = 0; prev < comp_i; prev++) {
        remove_component(z, &result.pc_vectors[prev * d], d);
      }

      /* normalize z -> new v */
      double normz = 0.0;
      for (size_t c = 0; c < d; c++) {
        normz += z[c] * z[c];
      }
      normz = sqrt(normz);
      if (normz < 1e-15) normz = 1.0;
      for (size_t c = 0; c < d; c++) {
        v[c] = z[c] / normz;
      }
    }

    /* store v as PC comp_i */
    memcpy(&result.pc_vectors[comp_i * d], v, d * sizeof(double));

    /* approximate eigenvalue => ||Xv||^2 / (n_inds - 1) */
    multiply_X_vec(ivars, n_inds, v, y);
    double sum_sq = 0.0;
    for (int i = 0; i < n_inds; i++) {
      sum_sq += y[i] * y[i];
    }
    result.eigenvalues[comp_i] = sum_sq / (double)(n_inds - 1);
  }

  /* 3) Build final scores => scores[i, comp_i] = dot( row i, pc_vectors[comp_i] ). */
  for (int comp_i = 0; comp_i < top_k; comp_i++) {
    const double *pcv = &result.pc_vectors[comp_i * d];
    multiply_X_vec(ivars, n_inds, pcv, y); /* y = X*pcv => length n_inds */
    for (int i = 0; i < n_inds; i++) {
      result.scores[i * top_k + comp_i] = y[i];
    }
  }

  /* 4) z-score each PC dimension => mean=0, stdev=1 => helps spread them out. */
  z_score_pca_scores(result.scores, n_inds, top_k);

  /* Cleanup temps */
  free(v);
  free(z);
  free(y);

  return result;
}
