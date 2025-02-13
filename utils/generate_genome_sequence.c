/**
 * @file generate_genome_sequence.c
 * @brief Generates a random DNA sequence with specified GC content
 * @author Namir Garib
 * @date 2025-01-12
 *
 * This program generates large DNA sequences with user-specified GC content
 * using a high-quality pseudorandom number generator.
 * Output is written in chunks for optimal I/O performance.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

/** Xorshift128+ PRNG state */
static uint64_t s[2];

/**
 * @brief Xorshift128+ pseudorandom number generator
 * @return 64-bit pseudorandom number
 */
static uint64_t xorshift128plus(void) {
  uint64_t       x = s[0];
  uint64_t const y = s[1];
  s[0] = y;
  x ^= x << 23;
  s[1] = x ^ y ^ (x >> 17) ^ (y >> 26);
  return s[1] + y;
}

#define BUFFER_SIZE (16 * 1024 * 1024)  // 16 MB buffer

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "Usage: %s filename num_bases gc_content\n", argv[0]);
    return EXIT_FAILURE;
  }

  const char    *filename = argv[1];
  const uint64_t num_bases = strtoull(argv[2], NULL, 10);
  const double   gc_content = strtod(argv[3], NULL);

  if (num_bases == 0) {
    fprintf(stderr, "Error: Number of bases must be positive\n");
    return EXIT_FAILURE;
  }
  if (gc_content < 0.0 || gc_content > 1.0) {
    fprintf(stderr, "Error: GC content must be between 0 and 1\n");
    return EXIT_FAILURE;
  }

  // Initialize PRNG with entropy from /dev/urandom or fallback
  FILE *urandom = fopen("/dev/urandom", "rb");
  if (urandom) {
    if (fread(s, sizeof(uint64_t), 2, urandom) != 2) {
      fprintf(stderr, "Warning: Using fallback seed\n");
      s[0] = (uint64_t)time(NULL) ^ ((uint64_t)getpid() << 32);
      s[1] = (uint64_t)time(NULL) << 32 | getpid();
    }
    fclose(urandom);
  } else {
    fprintf(stderr, "Warning: /dev/urandom unavailable. Using fallback seed\n");
    s[0] = (uint64_t)time(NULL) ^ ((uint64_t)getpid() << 32);
    s[1] = (uint64_t)time(NULL) << 32 | getpid();
  }

  uint32_t threshold = (uint32_t)(gc_content * 65536.0);
  threshold = threshold > 65535 ? 65535 : threshold;  // GC threshold (16-bit precision)

  FILE *output = fopen(filename, "w");
  if (!output) {
    perror("Error opening output file");
    return EXIT_FAILURE;
  }

  char *buffer = malloc(BUFFER_SIZE);  // Initialize buffer
  if (!buffer) {
    perror("Memory allocation failed");
    fclose(output);
    return EXIT_FAILURE;
  }
  size_t   buf_idx = 0;
  uint64_t remaining = num_bases;

  // Generate DNA sequence in 4-base chunks
  while (remaining >= 4) {
    uint64_t rand_val = xorshift128plus();
    uint16_t parts[4] = {(uint16_t)(rand_val >> 48), (uint16_t)(rand_val >> 32),
                         (uint16_t)(rand_val >> 16), (uint16_t)rand_val};

    for (int i = 0; i < 4; i++) {
      if (parts[i] < threshold) {
        buffer[buf_idx++] = (parts[i] & 1) ? 'C' : 'G';
      } else {
        buffer[buf_idx++] = (parts[i] & 1) ? 'T' : 'A';
      }
    }
    remaining -= 4;

    if (buf_idx >= BUFFER_SIZE) {
      fwrite(buffer, 1, buf_idx, output);
      buf_idx = 0;
    }
  }

  // Process remaining bases (0-3)
  if (remaining > 0) {
    uint64_t rand_val = xorshift128plus();
    for (int i = 0; i < 4 && remaining > 0; i++) {
      uint16_t part = (rand_val >> (48 - i * 16)) & 0xFFFF;
      if (part < threshold) {
        buffer[buf_idx++] = (part & 1) ? 'C' : 'G';
      } else {
        buffer[buf_idx++] = (part & 1) ? 'T' : 'A';
      }
      remaining--;
    }
  }

  if (buf_idx > 0) {
    fwrite(buffer, 1, buf_idx, output);
  }

  free(buffer);
  fclose(output);
  return EXIT_SUCCESS;
}