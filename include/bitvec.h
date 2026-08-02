#pragma once
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define BV_ALIGNEMENT 64

typedef uint64_t bv_limb;
#define BYTES_PER_LIMB ((size_t) sizeof(bv_limb))
#define BITS_PER_LIMB  (BYTES_PER_LIMB * 8)
#define BV_LIMBS_PER_BLOCK (BV_ALIGNEMENT / BYTES_PER_LIMB)

typedef struct
{
    size_t len;
    size_t nbits;
    bv_limb* rep;
} bitvec;

typedef bitvec bitvec_t[1];

// zero every bit at index >= nbits
static inline void bv_trim(bitvec_t bv)
{
    size_t full = bv->nbits / BITS_PER_LIMB;
    size_t rem  = bv->nbits % BITS_PER_LIMB;

    if (rem)
    {
        bv->rep[full] &= (((bv_limb) 1) << rem) - 1;
        full++;
    }

    if (full < bv->len)
    {
        memset(bv->rep + full, 0, (bv->len - full) * BYTES_PER_LIMB);
    }
}

static inline void bv_init(bitvec_t bv, size_t nbits, uint8_t byte_val)
{
    size_t limbs  = (nbits + BITS_PER_LIMB - 1) / BITS_PER_LIMB;
    size_t blocks = (limbs + BV_LIMBS_PER_BLOCK - 1) / BV_LIMBS_PER_BLOCK;

    if (blocks == 0)
    {
        blocks = 1;
    }

    bv->nbits = nbits;
    bv->len   = blocks * BV_LIMBS_PER_BLOCK;
    bv->rep = aligned_alloc(BV_ALIGNEMENT, bv->len * BYTES_PER_LIMB);
    memset(bv->rep, byte_val, bv->len * BYTES_PER_LIMB);

    if (byte_val)
    {
        bv_trim(bv);
    }
}

static inline void bv_clear(bitvec_t bv)
{
    free(bv->rep);
    bv->rep = NULL;
    bv->len = 0;
    bv->nbits = 0;
}

static inline void bv_zero(bitvec_t bv)
{
    memset(bv->rep, 0, bv->len * BYTES_PER_LIMB);
}

static inline void bv_set_bit(bitvec_t bv, size_t bit)
{
    assert(bit < bv->nbits);
    bv->rep[bit / BITS_PER_LIMB] |= ((bv_limb) 1) << (bit % BITS_PER_LIMB);
}

static inline void bv_clr_bit(bitvec_t bv, size_t bit)
{
    assert(bit < bv->nbits);
    bv->rep[bit / BITS_PER_LIMB] &= ~(((bv_limb) 1) << (bit % BITS_PER_LIMB));
}

static inline void bv_set(bitvec_t dest, const bitvec_t src)
{
    assert(dest->len >= src->len);
    memcpy(dest->rep, src->rep, src->len * BYTES_PER_LIMB);
    if (dest->len > src->len)
    {
        memset(dest->rep + src->len, 0, (dest->len - src->len) * BYTES_PER_LIMB);
    }
}

void bv_and(bitvec_t dest, const bitvec_t src1, const bitvec_t src2);
void bv_or(bitvec_t dest, const bitvec_t src1, const bitvec_t src2);

// dest &= src, returning 1 if any bit of the result is set.
int bv_and_into(bitvec_t dest, const bitvec_t src);

// dest |= src
void bv_or_into(bitvec_t dest, const bitvec_t src);

// dest = src << shift, with bit 0 the least significant
void bv_shift(bitvec_t dest, const bitvec_t src, size_t shift);

// dest |= src << shift, which is what replicating a periodic mask needs
void bv_or_shift_into(bitvec_t dest, const bitvec_t src, size_t shift);

// index of the lowest set bit, or BV_NO_BIT if there is none
#define BV_NO_BIT ((size_t) -1)
size_t bv_ctz(const bitvec_t src);

int bv_is_zero(const bitvec_t src);

// which kernel the build selected: "avx512", "avx2" or "scalar"
const char* bv_simd_name(void);
