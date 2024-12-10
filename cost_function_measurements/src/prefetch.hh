#pragma once

/*
 * Helper constants and macros for AMAC prefetch operations.
 *
 * void __builtin_prefetch (const void *addr, ...)
 *
 * arg1: addr:     address of the memory to prefetch
 * arg2: rw:       0 = read or 1 = write
 * arg3: locality: must be a compile-time constant integer between zero and three.
 *                 0 means that the data has no temporal locality, so it need not be left in the cache after the access.
 *                 3 means that the data has a high degree of temporal locality and should be left in all levels of cache possible.
 *                 1 or 2 mean, respectively, a low or moderate degree of temporal locality.
 *                 Default: 3.
 *
 */

static constexpr int PP_R = 0; // prefetch parameter read
static constexpr int PP_W = 1; // prefetch parameter write
static constexpr int PP_S = 0; // prefetch parameter single   access
static constexpr int PP_M = 1; // prefetch parameter multiple accesses
static constexpr int PP_B = PP_S; // build
static constexpr int PP_P = PP_M; // probe

#define PREFETCH_BUILD(ADDR) __builtin_prefetch(ADDR, PP_W, PP_B)
#define PREFETCH_PROBE(ADDR) __builtin_prefetch(ADDR, PP_R, PP_P)
