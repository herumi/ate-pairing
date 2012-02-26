/* -*- mode: c++ -*- */
#ifndef MIE_V128_H_
#define MIE_V128_H_

#include <string>

#ifdef _WIN32
	#include <intrin.h>
#else
	#include <xmmintrin.h>
	#include <tmmintrin.h>
#endif

struct V128 {
	__m128i x_;
	V128()
	{
	}
	V128(const uint32_t *p)
		: x_(*(const __m128i*)p)
	{
	}
	V128(__m128i x)
		: x_(x)
	{
	}
	V128(uint32_t x0, uint32_t x1, uint32_t x2, uint32_t x3)
		: x_(_mm_set_epi32(x0, x1, x2, x3))
	{
	}
	V128(uint32_t x)
		: x_(_mm_cvtsi32_si128(x))
	{
	}
	V128(const V128& rhs)
		: x_(rhs.x_)
	{
	}
	void clear()
	{
		*this = _mm_setzero_si128();
	}
	void set(uint32_t x)
	{
		x_ = _mm_set1_epi32(x);
	}
	V128& operator|=(const V128& rhs)
	{
		x_ = _mm_or_si128(x_, rhs.x_);
		return *this;
	}
	V128& operator&=(const V128& rhs)
	{
		x_ = _mm_and_si128(x_, rhs.x_);
		return *this;
	}
	V128& operator^=(const V128& rhs)
	{
		x_ = _mm_xor_si128(x_, rhs.x_);
		return *this;
	}
	V128& andn(const V128& rhs)
	{
		x_ = _mm_andnot_si128(x_, rhs.x_);
		return *this;
	}
	MIE_FORCE_INLINE void store(uint32_t *p) const
	{
		*(__m128i*)p = x_;
	}
	uint32_t get32bit() const
	{
		return _mm_cvtsi128_si32(x_);
	}
	void load(const uint32_t *p)
	{
		x_ = *(const __m128i*)p;
	}
	/*
		*this >>= n
	*/
	template<int n>
	void shrBit();
	/*
		*this <<= n
	*/
	template<int n>
	void shlBit();
	uint32_t pmovmskb() const { return _mm_movemask_epi8(x_); }
	V128& cmpeq(const V128& rhs)
	{
		x_ = _mm_cmpeq_epi32(x_, rhs.x_);
		return *this;
	}
};

static MIE_FORCE_INLINE V128 Zero()
{
	return _mm_setzero_si128();
}

template<int n>
static MIE_FORCE_INLINE V128 psrldq(const V128& a)
{
	return _mm_srli_si128(a.x_, n);
}

template<int n>
static MIE_FORCE_INLINE V128 pslldq(const V128& a)
{
	return _mm_slli_si128(a.x_, n);
}

template<int n>
static MIE_FORCE_INLINE V128 psrlq(const V128& a)
{
	return _mm_srli_epi64(a.x_, n);
}

static MIE_FORCE_INLINE V128 psrlq(const V128& a, const V128& n)
{
	return _mm_srl_epi64(a.x_, n.x_);
}

template<int n>
static MIE_FORCE_INLINE V128 psllq(const V128& a)
{
	return _mm_slli_epi64(a.x_, n);
}

static MIE_FORCE_INLINE V128 psllq(const V128& a, const V128& n)
{
	return _mm_sll_epi64(a.x_, n.x_);
}

template<int n>
static MIE_FORCE_INLINE V128 palignr(const V128& a, const V128& b)
{
#ifdef MIE_USE_SSSE3
	return _mm_alignr_epi8(a.x_, b.x_, n);
#else
	return pslldq<16 - n>(a) | psrldq<n>(b);
#endif
}

static MIE_FORCE_INLINE V128 punpckhqdq(const V128& a, const V128& b)
{
	return _mm_unpackhi_epi64(a.x_, b.x_);
}

static MIE_FORCE_INLINE V128 paddd(const V128& a, const V128& b)
{
	return _mm_add_epi32(a.x_, b.x_);
}

static MIE_FORCE_INLINE V128 andn(const V128& lhs, const V128& rhs)
{
	return V128(lhs).andn(rhs);
}

static MIE_FORCE_INLINE V128 operator|(const V128& lhs, const V128& rhs)
{
	return V128(lhs) |= rhs;
}

static MIE_FORCE_INLINE V128 operator&(const V128& lhs, const V128& rhs)
{
	return V128(lhs) &= rhs;
}

static MIE_FORCE_INLINE V128 operator^(const V128& lhs, const V128& rhs)
{
	return V128(lhs) ^= rhs;
}

static MIE_FORCE_INLINE void swap128(uint32_t *p, uint32_t *q)
{
	V128 t(p);
	V128(q).store(p);
	t.store(q);
}

static MIE_FORCE_INLINE void copy128(uint32_t *dest, const uint32_t *src)
{
	V128(src).store(dest);
}

template<int n>
void V128::shrBit()
{
	assert(n < 64);
	*this = psrlq<n>(*this) | psllq<64 - n>(psrldq<8>(*this));
#if 0
	if (n == 64) {
		*this = psrldq<8>(*this);
	} else if (n > 64) {
		*this = psrlq<n - 64>(psrldq<8>(*this));
	}
#endif
}
template<int n>
void V128::shlBit()
{
	assert(n < 64);
	*this = psllq<n>(*this) | psrlq<64 - n>(pslldq<8>(*this));
#if 0
	if (n == 64) {
		*this = pslldq<8>(*this);
	} else if (n > 64) {
		*this = psllq<n - 64>(pslldq<8>(*this));
	}
#endif
}

/*
	byte rotr [x2:x1:x0]
*/
template<int n>
static MIE_FORCE_INLINE void rotrByte(V128& x0, V128& x1, V128& x2)
{
	V128 t(x0);
	x0 = palignr<n>(x1, x0);
	x1 = palignr<n>(x2, x1);
	x2 = palignr<n>(t, x2);
}

/*
	byte rotl [x2:x1:x0]
*/
template<int n>
static MIE_FORCE_INLINE void rotlByte(V128& x0, V128& x1, V128& x2)
{
	V128 t(x2);
	x2 = palignr<16 - n>(x2, x1);
	x1 = palignr<16 - n>(x1, x0);
	x0 = palignr<16 - n>(x0, t);
}

//MIE_FORCE_INLINE V128::V128(uint32_t x0, uint32_t x1, uint32_t x2, uint32_t x3)
//{
//	x_ = (pslldq<12>(V128(x0)) | pslldq<8>(V128(x1)) | pslldq<4>(V128(x2)) | V128(x3)).x_;
//}

/**
	index = bsf(val)
	return val != 0
*/
#ifdef _WIN32
	#define mie_bsf(pindex, val) _BitScanForward(reinterpret_cast<unsigned long *>(pindex), val)
#else
	#define mie_bsf(pindex, val) val ? *pindex = __builtin_ctz(val), 1 : 0
#endif
/**
	index = bsr(val)
	return val != 0
*/
#ifdef _WIN32
	#define mie_bsr(pindex, val) _BitScanReverse(reinterpret_cast<unsigned long *>(pindex), val)
#else
	#define mie_bsr(pindex, val) val ? (*pindex = __builtin_clz(val) ^ 0x1f), 1 : 0
#endif

#endif /* MIE_V128_H_ */

/*
  Local Variables:
  c-basic-offset: 4
  indent-tabs-mode: t
  tab-width: 4
  End:
*/
