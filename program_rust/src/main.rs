/*!
 * @file main.rs
 * @brief Main entry point for PCA-based population stratification (Rust version)
 *
 * Author: Namir Garib
 * Created: January 2025
 */

mod analysis;
mod pca;
mod utils;
mod variant_calling;

use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!(
            "Usage: {} <ref_genome> <num_individuals> <indiv1> [indiv2 ...]",
            args[0]
        );
        process::exit(1);
    }

    let ref_file = &args[1];
    let num_individuals: usize = args[2].parse().unwrap_or(0);
    if num_individuals < 1 || (args.len() - 3) < num_individuals {
        eprintln!("Invalid number of individuals or not enough file paths.");
        process::exit(1);
    }

    let individuals_files = &args[3..(3 + num_individuals)];

    analysis::perform_full_analysis(ref_file, individuals_files);

    println!("Analysis complete. Check results.csv and eigenvalues.csv.");
}
