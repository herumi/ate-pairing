#pragma once
/**
	@file
	@brief square root(a part of mie/include/gmp_util.hpp; this code should be unified with it)
	@author MITSUNARI Shigeo(@herumi)
	@license modified new BSD license
	http://opensource.org/licenses/BSD-3-Clause
*/
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <assert.h>
#ifdef _MSC_VER
	#pragma warning(push)
	#pragma warning(disable : 4616)
	#pragma warning(disable : 4800)
	#pragma warning(disable : 4244)
	#pragma warning(disable : 4127)
	#pragma warning(disable : 4512)
	#pragma warning(disable : 4146)
#endif
#include <gmpxx.h>
#include <stdint.h>
#ifdef _MSC_VER
	#pragma warning(pop)
#endif
#ifdef _MSC_VER
#if _MSC_VER == 1900
#ifdef _DEBUG
#pragma comment(lib, "14/mpird.lib")
#pragma comment(lib, "14/mpirxxd.lib")
#else
#pragma comment(lib, "14/mpir.lib")
#pragma comment(lib, "14/mpirxx.lib")
#endif
#elif _MSC_VER == 1800
#ifdef _DEBUG
#pragma comment(lib, "12/mpird.lib")
#pragma comment(lib, "12/mpirxxd.lib")
#else
#pragma comment(lib, "12/mpir.lib")
#pragma comment(lib, "12/mpirxx.lib")
#endif
#else
#ifdef _DEBUG
#pragma comment(lib, "mpird.lib")
#pragma comment(lib, "mpirxxd.lib")
#else
#pragma comment(lib, "mpir.lib")
#pragma comment(lib, "mpirxx.lib")
#endif
#endif
#endif

namespace bn {

/*
	Tonelli-Shanks
*/
class SquareRoot {
	bool isValid;
	mpz_class p;
	mpz_class g;
	int r;
	mpz_class q; // p - 1 = 2^r q
	mpz_class s; // s = g^q
public:
	SquareRoot() : isValid(false) {}
	void set(const mpz_class& p)
	{
		if (p <= 2) throw std::runtime_error("SquareRoot:bad p");
		isValid = isPrime(p);
		if (!isValid) return; // don't throw until get() is called
		this->p = p;
		// g is quadratic nonresidue
		g = 2;
		while (legendre(g, p) > 0) {
			g++;
		}
		// p - 1 = 2^r q, q is odd
		r = 0;
		q = p - 1;
		while ((q & 1) == 0) {
			r++;
			q /= 2;
		}
		powMod(s, g, q, p);
	}
	/*
		solve x^2 = a mod p
	*/
	bool get(mpz_class& x, const mpz_class& a) const
	{
		if (!isValid) throw std::runtime_error("SquareRoot:get:not prime");
		if (legendre(a, p) < 0) return false;
		if (r == 1) {
			powMod(x, a, (p + 1) / 4, p);
			return true;
		}
		mpz_class c = s, d;
		int e = r;
		powMod(d, a, q, p);
		powMod(x, a, (q + 1) / 2, p); // destroy a if &x == &a
		while (d != 1) {
			int i = 1;
			mpz_class dd = (d * d) % p;
			while (dd != 1) {
				dd = (dd * dd) % p;
				i++;
			}
			mpz_class b = 1;
			b <<= e - i - 1;
			powMod(b, c, b, p);
			x = (x * b) % p;
			c = (b * b) % p;
			d = (d * c) % p;
			e = i;
		}
		return true;
	}
	// z = x^y mod m (y >=0)
	static inline void powMod(mpz_class& z, const mpz_class& x, const mpz_class& y, const mpz_class& m)
	{
		mpz_powm(z.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t(), m.get_mpz_t());
	}
	/*
		assume p : odd prime
		return  1 if x^2 = a mod p for some x
		return -1 if x^2 != a mod p for any x
	*/
	static inline int legendre(const mpz_class& a, const mpz_class& p)
	{
		return mpz_legendre(a.get_mpz_t(), p.get_mpz_t());
	}
	static inline bool isPrime(const mpz_class& x)
	{
		return mpz_probab_prime_p(x.get_mpz_t(), 25) != 0;
	}
};

} // mie
