/*!
 * @file utils.rs
 * @brief Utility functions for secure memory handling and file operations.
 *
 * Author: Namir Garib
 * Created: January 2025
 */

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

/**
 * @brief Get the file length in bytes.
 *
 * @param path The path to the file
 * @return The file length in bytes, or 0 if an error occurs.
 */
pub fn get_file_length(path: &str) -> usize {
    let file_path = Path::new(path);
    match File::open(file_path) {
        Ok(mut file) => {
            if let Ok(metadata) = file.metadata() {
                return metadata.len() as usize;
            }
        }
        Err(e) => {
            eprintln!("get_file_length: Failed to open file {}: {}", path, e);
            return 0;
        }
    }
    0
}

/**
 * @brief Read the entire file into a Vec<u8> in chunks.
 *
 * @param path   The file path.
 * @param length Number of bytes to read (assumes we know the file size).
 * @return A Result<Vec<u8>, String> containing the file data or an error message.
 *
 * For extremely large files, consider memory mapping or streaming approach.
 */
pub fn read_file_in_chunks(path: &str, length: usize) -> Result<Vec<u8>, String> {
    let mut file = match File::open(path) {
        Ok(f) => f,
        Err(e) => return Err(format!("Failed to open file {}: {}", path, e)),
    };

    let mut buffer = vec![0u8; length];
    let chunk_size = 1024 * 1024; // 1MB
    let mut total_read = 0;

    while total_read < length {
        let to_read = std::cmp::min(chunk_size, length - total_read);
        match file.read(&mut buffer[total_read..total_read + to_read]) {
            Ok(n) => {
                if n == 0 {
                    return Err("Unexpected EOF".to_string());
                }
                total_read += n;
            }
            Err(e) => return Err(format!("Error reading file: {}", e)),
        }
    }

    Ok(buffer)
}
