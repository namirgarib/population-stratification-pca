/*!
 * @file analysis.rs
 * @brief High-level pipeline: read reference, read individuals, variant call, PCA.
 *
 * Author: Namir Garib
 * Created: January 2025
 */

use crate::pca::{center_data, compute_covariance_matrix, eigen_decomposition, project_data};
use crate::utils::{get_file_length, read_file_in_chunks};
use crate::variant_calling::call_variants;
use std::fs::File;
use std::io::Write;

/**
 * @brief Orchestrates the pipeline for multiple individuals.
 *
 * @param ref_file        Path to the reference genome.
 * @param individuals_files A slice of paths to individuals' genomes.
 */
pub fn perform_full_analysis(ref_file: &str, individuals_files: &[String]) {
    let ref_length = get_file_length(ref_file);
    if ref_length == 0 {
        eprintln!(
            "Reference genome size is 0 or error reading file: {}",
            ref_file
        );
        return;
    }
    println!("Reference genome length: {}", ref_length);

    let ref_data = match read_file_in_chunks(ref_file, ref_length) {
        Ok(buf) => buf,
        Err(e) => {
            eprintln!("Error reading reference file: {}", e);
            return;
        }
    };

    let n = individuals_files.len();
    let d = ref_length;

    // Prepare data matrix for variant calls: n x d
    let mut data_matrix: Vec<f64> = Vec::with_capacity(n * d);

    // For each individual, call variants
    for (i, indiv) in individuals_files.iter().enumerate() {
        let indiv_length = get_file_length(indiv);
        if indiv_length != d {
            eprintln!(
                "Individual {} length {} != reference length {}",
                i, indiv_length, d
            );
            return;
        }

        let indiv_data = match read_file_in_chunks(indiv, indiv_length) {
            Ok(buf) => buf,
            Err(e) => {
                eprintln!("Error reading individual file {}: {}", indiv, e);
                return;
            }
        };

        let variants = call_variants(&ref_data, &indiv_data);
        data_matrix.extend_from_slice(&variants);
    }

    // Perform PCA (n = number of individuals, d = length of genome)
    let centered = center_data(&data_matrix, n, d);
    let cov = compute_covariance_matrix(&centered, n, d);
    let pca_res = eigen_decomposition(&cov, d);
    let scores = project_data(&centered, n, d, &pca_res);

    // Write results
    {
        let mut f_scores = match File::create("results.csv") {
            Ok(file) => file,
            Err(e) => {
                eprintln!("Failed to create results.csv: {}", e);
                return;
            }
        };
        for row in 0..n {
            for comp in 0..d {
                let val = scores[row * d + comp];
                if comp < (d - 1) {
                    write!(f_scores, "{:.6},", val).unwrap();
                } else {
                    write!(f_scores, "{:.6}\n", val).unwrap();
                }
            }
        }
    }

    {
        let mut f_evals = match File::create("eigenvalues.csv") {
            Ok(file) => file,
            Err(e) => {
                eprintln!("Failed to create eigenvalues.csv: {}", e);
                return;
            }
        };
        for (i, &val) in pca_res.eigenvalues.iter().enumerate() {
            writeln!(f_evals, "{},{}", i + 1, val).unwrap();
        }
    }

    println!("PCA analysis completed. See results.csv and eigenvalues.csv");
}
