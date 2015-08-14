// This file is a part of Julia. License is MIT: http://julialang.org/license

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include "dtypes.h"
#include "julia.h"

#ifdef __cplusplus
extern "C" {
#endif

jl_bytebuffer_t *jl_bytebuffer_new(size_t size)
{
    jl_bytebuffer_t *b = (jl_bytebuffer_t *)LLT_ALLOC(B_N_INLINE
        + sizeof(uint8_t *) + 2 * sizeof(size_t));
    b->len = 0;
    if (size <= B_N_INLINE) {
        b->items = &b->_space[0];
        b->max = B_N_INLINE;
    }
    else {
        b->items = (uint8_t*)LLT_ALLOC(size);
        b->max = size;
    }
    if (b->items == NULL) return NULL;
    return b;
}

void jl_bytebuffer_free(jl_bytebuffer_t *b)
{
    if (b->items != &b->_space[0])
        LLT_FREE(b->items);
    b->len = 0;
    b->max = B_N_INLINE;
    b->items = &b->_space[0];
}

#ifdef __cplusplus
}
#endif
