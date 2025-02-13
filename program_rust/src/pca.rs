/*!
 * @file pca.rs
 * @brief Principal Component Analysis (PCA) routines in Rust
 *
 * Author: Namir Garib
 * Created: January 2025
 */

use std::f64;

/**
 * @struct PCAResult
 * @brief Container for PCA results (eigenvalues, eigenvectors).
 */
pub struct PCAResult {
    pub eigenvalues: Vec<f64>,
    pub eigenvectors: Vec<f64>,
    pub num_components: usize,
    pub dimension: usize,
}

/**
 * @brief Center data column-wise (subtract mean from each column).
 *
 * @param data   Data in row-major format: n x d
 * @param n      Number of samples
 * @param d      Dimension (number of SNP positions)
 * @return A new Vec<f64> containing the centered data.
 */
pub fn center_data(data: &Vec<f64>, n: usize, d: usize) -> Vec<f64> {
    let mut centered = vec![0.0; n * d];
    let mut means = vec![0.0; d];

    // Compute column means
    for col in 0..d {
        let mut sum = 0.0;
        for row in 0..n {
            sum += data[row * d + col];
        }
        means[col] = sum / (n as f64);
    }

    // Subtract means
    for row in 0..n {
        for col in 0..d {
            centered[row * d + col] = data[row * d + col] - means[col];
        }
    }

    centered
}

/**
 * @brief Compute covariance matrix (d x d).
 *
 * @param centered_data The centered data (n x d).
 * @param n             Number of samples.
 * @param d             Dimension.
 * @return Vec<f64>     A new vector storing the covariance matrix in row-major order.
 */
pub fn compute_covariance_matrix(centered_data: &Vec<f64>, n: usize, d: usize) -> Vec<f64> {
    let mut cov = vec![0.0; d * d];
    for i in 0..d {
        for j in 0..d {
            let mut sum = 0.0;
            for k in 0..n {
                sum += centered_data[k * d + i] * centered_data[k * d + j];
            }
            cov[i * d + j] = sum / ((n - 1) as f64);
        }
    }
    cov
}

/**
 * @brief Naive eigen decomposition for demonstration.
 *        In real usage, use a numeric library for large d.
 */
pub fn eigen_decomposition(cov_matrix: &Vec<f64>, d: usize) -> PCAResult {
    let mut eigenvectors = cov_matrix.clone();
    let mut eigenvalues = vec![0.0; d];

    // Pretend diagonal are eigenvalues
    for i in 0..d {
        eigenvalues[i] = eigenvectors[i * d + i];
    }

    // Sort in descending order (bubble sort demonstration)
    for i in 0..(d - 1) {
        for j in 0..(d - i - 1) {
            if eigenvalues[j] < eigenvalues[j + 1] {
                eigenvalues.swap(j, j + 1);
                // swap the entire row in eigenvectors
                for col in 0..d {
                    let idx1 = j * d + col;
                    let idx2 = (j + 1) * d + col;
                    eigenvectors.swap(idx1, idx2);
                }
            }
        }
    }

    PCAResult {
        eigenvalues,
        eigenvectors,
        num_components: d,
        dimension: d,
    }
}

/**
 * @brief Project data onto the principal components.
 *
 * @param centered_data The centered data (n x d).
 * @param n             Number of samples.
 * @param d             Dimension.
 * @param pca_result    Contains eigenvectors and eigenvalues.
 * @return Vec<f64>     The projected data (n x num_components).
 */
pub fn project_data(
    centered_data: &Vec<f64>,
    n: usize,
    d: usize,
    pca_result: &PCAResult,
) -> Vec<f64> {
    let dims_to_use = pca_result.num_components;
    let mut projections = vec![0.0; n * dims_to_use];

    for row in 0..n {
        for comp in 0..dims_to_use {
            let mut sum = 0.0;
            for col in 0..d {
                sum += centered_data[row * d + col] * pca_result.eigenvectors[comp * d + col];
            }
            projections[row * dims_to_use + comp] = sum;
        }
    }

    projections
}
