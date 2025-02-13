/**
 * @file utils.c
 * @brief Utility functions for secure memory allocation
 *
 * Author: Namir Garib
 * Created: January 2025
 */

#include "utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void *secure_alloc(size_t size) {
  if (size == 0) {
    fprintf(stderr, "secure_alloc: Requested size is 0.\n");
    exit(EXIT_FAILURE);
  }
  void *ptr = calloc(1, size);
  if (!ptr) {
    fprintf(stderr, "secure_alloc: Out of memory requesting %zu bytes.\n", size);
    exit(EXIT_FAILURE);
  }
  return ptr;
}
