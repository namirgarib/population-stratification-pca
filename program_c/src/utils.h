#ifndef UTILS_H
#define UTILS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Allocates `size` bytes (zero-initialized). Exits on failure.
 */
void *secure_alloc(size_t size);

#ifdef __cplusplus
}
#endif

#endif
