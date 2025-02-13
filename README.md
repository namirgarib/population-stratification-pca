# population-stratification-pca

C and Rust implementation of principal component analysis (PCA) used for population stratification

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Introduction

This project provides implementations of Principal Component Analysis (PCA) in both C and Rust, specifically designed for population stratification. PCA is a statistical procedure that transforms possibly correlated variables into a set of linearly uncorrelated variables called principal components.

## Features

- Efficient PCA computation in C and Rust
- Designed for large-scale genomic data
- Supports various input formats
- High performance and low memory usage

## Installation

### Prerequisites

- C compiler (e.g., GCC)
- Rust toolchain (e.g., rustc, cargo)

### Building from Source

#### C Implementation

```sh
cd program_c
make
```

#### Rust Implementation

```sh
cd program_rust
cargo build --release
```

## Usage

### C Implementation

```sh
./program_c/pca input_data.txt output_data.txt
```

### Rust Implementation

```sh
./program_rust/target/release/pca input_data.txt output_data.txt
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any improvements or bug fixes.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
