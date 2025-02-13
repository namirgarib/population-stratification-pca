/*!
 * @file variant_calling.rs
 * @brief Variant calling logic (naive: 0 if same base, 1 if different).
 *
 * Author: Namir Garib
 * Created: January 2025
 */

/**
 * @brief Compare a reference genome and individual's genome to produce a variant vector.
 *
 * @param ref_genome    A slice of bytes for the reference genome.
 * @param indiv_genome  A slice of bytes for the individual's genome.
 * @return Vec<f64>     0.0 if same base, 1.0 if different base (naive).
 */
pub fn call_variants(ref_genome: &[u8], indiv_genome: &[u8]) -> Vec<f64> {
    let length = ref_genome.len();
    let mut variants = Vec::with_capacity(length);

    for i in 0..length {
        if ref_genome[i] == indiv_genome[i] {
            variants.push(0.0);
        } else {
            variants.push(1.0);
        }
    }

    variants
}
