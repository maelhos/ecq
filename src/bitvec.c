#include "bitvec.h"

#if defined(__AVX512F__)

#include <immintrin.h>
#define BV_SIMD "avx512"
typedef __m512i bv_vec;
#define BV_VEC_LIMBS   8
#define BV_LOAD(p)     _mm512_load_si512((const void*) (p))
#define BV_LOADU(p)    _mm512_loadu_si512((const void*) (p))
#define BV_STORE(p, v) _mm512_store_si512((void*) (p), (v))
#define BV_STOREU(p,v) _mm512_storeu_si512((void*) (p), (v))
#define BV_AND(a, b)   _mm512_and_si512((a), (b))
#define BV_OR(a, b)    _mm512_or_si512((a), (b))
#define BV_SLLI(a, n)  _mm512_slli_epi64((a), (n))
#define BV_SRLI(a, n)  _mm512_srli_epi64((a), (n))
#define BV_ZERO()      _mm512_setzero_si512()
#define BV_IS_ZERO(a)  (_mm512_test_epi64_mask((a), (a)) == 0)

#elif defined(__AVX2__)

#include <immintrin.h>
#define BV_SIMD "avx2"
typedef __m256i bv_vec;
#define BV_VEC_LIMBS   4
#define BV_LOAD(p)     _mm256_load_si256((const __m256i*) (p))
#define BV_LOADU(p)    _mm256_loadu_si256((const __m256i*) (p))
#define BV_STORE(p, v) _mm256_store_si256((__m256i*) (p), (v))
#define BV_STOREU(p,v) _mm256_storeu_si256((__m256i*) (p), (v))
#define BV_AND(a, b)   _mm256_and_si256((a), (b))
#define BV_OR(a, b)    _mm256_or_si256((a), (b))
#define BV_SLLI(a, n)  _mm256_slli_epi64((a), (n))
#define BV_SRLI(a, n)  _mm256_srli_epi64((a), (n))
#define BV_ZERO()      _mm256_setzero_si256()
#define BV_IS_ZERO(a)  _mm256_testz_si256((a), (a))

#else

#define BV_SIMD "scalar"
typedef bv_limb bv_vec;
#define BV_VEC_LIMBS   1
#define BV_LOAD(p)     (*(const bv_limb*) (p))
#define BV_LOADU(p)    (*(const bv_limb*) (p))
#define BV_STORE(p, v) (*(bv_limb*) (p) = (v))
#define BV_STOREU(p,v) (*(bv_limb*) (p) = (v))
#define BV_AND(a, b)   ((a) & (b))
#define BV_OR(a, b)    ((a) | (b))
#define BV_SLLI(a, n)  ((a) << (n))
#define BV_SRLI(a, n)  ((a) >> (n))
#define BV_ZERO()      ((bv_limb) 0)
#define BV_IS_ZERO(a)  ((a) == 0)

#endif

const char* bv_simd_name(void)
{
    return BV_SIMD;
}

void bv_and(bitvec_t dest, const bitvec_t src1, const bitvec_t src2)
{
    size_t n = src1->len < src2->len ? src1->len : src2->len;
    size_t i;

    assert(dest->len >= n);

    for (i = 0; i < n; i += BV_VEC_LIMBS)
    {
        BV_STORE(dest->rep + i, BV_AND(BV_LOAD(src1->rep + i), BV_LOAD(src2->rep + i)));
    }

    // past the shorter operand the AND is zero
    if (dest->len > n)
    {
        memset(dest->rep + n, 0, (dest->len - n) * BYTES_PER_LIMB);
    }
}

void bv_or(bitvec_t dest, const bitvec_t src1, const bitvec_t src2)
{
    const bitvec *hi = src1->len < src2->len ? src2 : src1;
    size_t n = src1->len < src2->len ? src1->len : src2->len;
    size_t i;

    assert(dest->len >= hi->len);

    for (i = 0; i < n; i += BV_VEC_LIMBS)
    {
        BV_STORE(dest->rep + i, BV_OR(BV_LOAD(src1->rep + i), BV_LOAD(src2->rep + i)));
    }

    for (; i < hi->len; i += BV_VEC_LIMBS)
    {
        BV_STORE(dest->rep + i, BV_LOAD(hi->rep + i));
    }

    if (dest->len > hi->len)
    {
        memset(dest->rep + hi->len, 0, (dest->len - hi->len) * BYTES_PER_LIMB);
    }
}

int bv_and_into(bitvec_t dest, const bitvec_t src)
{
    size_t n = dest->len < src->len ? dest->len : src->len;
    bv_vec acc = BV_ZERO();
    size_t i;

    // this is the sieve's hot loop: the running OR keeps the emptiness test free
    for (i = 0; i < n; i += BV_VEC_LIMBS)
    {
        bv_vec v = BV_AND(BV_LOAD(dest->rep + i), BV_LOAD(src->rep + i));
        BV_STORE(dest->rep + i, v);
        acc = BV_OR(acc, v);
    }

    if (dest->len > n)
    {
        memset(dest->rep + n, 0, (dest->len - n) * BYTES_PER_LIMB);
    }

    return !BV_IS_ZERO(acc);
}

void bv_or_into(bitvec_t dest, const bitvec_t src)
{
    size_t n = dest->len < src->len ? dest->len : src->len;
    size_t i;

    for (i = 0; i < n; i += BV_VEC_LIMBS)
    {
        BV_STORE(dest->rep + i, BV_OR(BV_LOAD(dest->rep + i), BV_LOAD(src->rep + i)));
    }
}

int bv_is_zero(const bitvec_t src)
{
    bv_vec acc = BV_ZERO();
    size_t i;

    for (i = 0; i < src->len; i += BV_VEC_LIMBS)
    {
        acc = BV_OR(acc, BV_LOAD(src->rep + i));
    }

    return BV_IS_ZERO(acc);
}

static inline void bv_shift_impl(bitvec_t dest, const bitvec_t src, size_t shift, int accumulate)
{
    size_t off = shift / BITS_PER_LIMB;
    unsigned bit = (unsigned) (shift % BITS_PER_LIMB);
    size_t lo, hi, i;

    lo = off + 1;
    hi = src->len + off;
    if (hi > dest->len)
    {
        hi = dest->len;
    }

    i = dest->len;
    while (i-- > hi)
    {
        if (!accumulate)
        {
            dest->rep[i] = 0;
        }
    }

    if (hi > lo)
    {
        size_t nvec = (hi - lo) / BV_VEC_LIMBS;
        size_t vec_lo = hi - nvec * BV_VEC_LIMBS;

        i = hi;
        while (i > vec_lo)
        {
            i -= BV_VEC_LIMBS;

            if (bit == 0)
            {
                bv_vec v = BV_LOADU(src->rep + (i - off));
                BV_STOREU(dest->rep + i, accumulate ? BV_OR(BV_LOADU(dest->rep + i), v) : v);
            }
            else
            {
                bv_vec v_hi = BV_LOADU(src->rep + (i - off));
                bv_vec v_lo = BV_LOADU(src->rep + (i - off - 1));
                bv_vec v = BV_OR(BV_SLLI(v_hi, bit), BV_SRLI(v_lo, BITS_PER_LIMB - bit));
                BV_STOREU(dest->rep + i, accumulate ? BV_OR(BV_LOADU(dest->rep + i), v) : v);
            }
        }

        hi = vec_lo;
    }

    // scalar remainder, still descending
    i = hi;
    while (i-- > 0)
    {
        bv_limb v = 0;

        if (i >= off)
        {
            size_t j = i - off;

            if (j < src->len)
            {
                v = src->rep[j] << bit;
            }

            if (bit != 0 && j > 0 && j - 1 < src->len)
            {
                v |= src->rep[j - 1] >> (BITS_PER_LIMB - bit);
            }
        }

        if (accumulate)
        {
            dest->rep[i] |= v;
        }
        else
        {
            dest->rep[i] = v;
        }
    }

    bv_trim(dest);
}

void bv_shift(bitvec_t dest, const bitvec_t src, size_t shift)
{
    bv_shift_impl(dest, src, shift, 0);
}

// dest |= src << shift, which is exactly what replicating a periodic mask needs
void bv_or_shift_into(bitvec_t dest, const bitvec_t src, size_t shift)
{
    bv_shift_impl(dest, src, shift, 1);
}

size_t bv_ctz(const bitvec_t src)
{
    size_t i;

    for (i = 0; i < src->len; i++)
    {
        if (src->rep[i] != 0)
        {
            // the limb offset is added to the trailing zero count, not to the limb
            return (size_t) __builtin_ctzll(src->rep[i]) + i * BITS_PER_LIMB;
        }
    }

    return BV_NO_BIT;
}
