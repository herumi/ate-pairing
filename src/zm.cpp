#include "zm.h"
#include "xbyak/xbyak.h"
#include <cstdio>

using namespace mie;
using namespace Xbyak;
/**
	out[] = x[] + y[]
	@note the sizeof out >= n
	@return size of x[] + y[]
*/
static inline bool in_addN(Unit *out, const Unit *x, const Unit *y, size_t n)
{
	printf("in_addN\n");
	assert(n > 0);

	Unit c = 0;
	for (size_t i = 0; i < n; i++) {
		Unit xc = x[i] + c;
		if (xc < c) {
			// x[i] = Unit(-1) and c = 1
			out[i] = y[i];
		} else {
			xc += y[i];
			c = y[i] > xc ? 1 : 0;
			out[i] = xc;
		}
	}
	return c != 0;
}

/**
	out[] = x[] + y
*/
static inline bool in_add(Unit *out, const Unit *x, size_t n, Unit y)
{
	assert(n > 0);

	Unit xc = x[0] + y;
	Unit c = y > xc ? 1 : 0;
	out[0] = xc;
	for (size_t i = 1; i < n; i++) {
		Unit xc = x[i] + c;
		if (xc < c) {
			out[i] = 0;
		} else {
			out[i] = xc;
			c = 0;
		}
	}
	return c != 0;
}
/**
	out[] = x[] - y[]
*/
static inline bool in_subN(Unit *out, const Unit *x, const Unit *y, size_t n)
{
	assert(n > 0);
	Unit c = 0;
	for (size_t i = 0; i < n; i++) {
		Unit yc = y[i] + c;
		if (yc < c) {
			// y[i] = Unit(-1) and c = 1
			out[i] = x[i];
		} else {
			c = x[i] < yc ? 1 : 0;
			out[i] = x[i] - yc;
		}
	}
	return c != 0;
}

/**
	out[] = x[] - y
*/
static inline bool in_sub(Unit *out, const Unit *x, size_t n, Unit y)
{
	assert(n > 0);
	Unit c = x[0] < y ? 1 : 0;
	out[0] = x[0] - y;
	for (size_t i = 1; i < n; i++) {
		if (x[i] < c) {
			out[i] = Unit(-1);
		} else {
			out[i] = x[i] - c;
			c = 0;
		}
	}
	return c != 0;
}

/*
	[H:L] <= a * b
	@return L
*/
static inline Unit mulUnit(Unit *H, Unit a, Unit b)
{
#ifdef MIE_USE_UNIT32
#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
	uint64_t t = __emulu(a, b);
#else
	uint64_t t = uint64_t(a) * b;
#endif
	uint32_t L;
	split64(H, &L, t);
	return L;
#else
#if defined(_WIN64) && !defined(__INTEL_COMPILER)
	return _umul128(a, b, H);
#else
	fprintf(stderr, "not implemented mulUnit %p %d %d\n", H, (int)a, (int)b);
	exit(1);
	// use _mm_mul_epu32(a, b)
#endif
#endif
}

/*
	out[0..n + 1] = x[0..n] * y
	@note accept out == x
*/
static inline void in_mul(Unit *out, const Unit *x, size_t n, Unit y)
{
	printf("in_mul\n");

	assert(n > 0);
	Unit H = 0;
	for (size_t i = 0; i < n; i++) {
		Unit t = H;
		Unit L = mulUnit(&H, x[i], y);
		out[i] = t + L;
		if (out[i] < t) {
			H++;
		}
	}
	out[n] = H;
}

/*
	q = [H:L] / y
	r = [H:L] % y
	return q
*/
#ifdef MIE_USE_UNIT32
static inline Unit divUnit(Unit *r, Unit H, Unit L, Unit y)
{
	uint64_t t = make64(H, L);
	uint32_t q = uint32_t(t / y);
	*r = t % y;
	return q;
}
#else
static inline Unit divUnit(Unit *, Unit, Unit, Unit)
{
	fprintf(stderr, "not implemented divUnit\n");
	exit(1);
}
#endif

/*
	q = x[] / y
	@retval r = x[] % y
	@note accept q == x
*/
static inline Unit in_div(Unit *q, const Unit *x, size_t xn, Unit y)
{
	Unit r = 0;
	for (int i = (int)xn - 1; i >= 0; i--) {
		q[i] = divUnit(&r, r, x[i], y);
	}
	return r;
}

static inline Unit in_mod(const Unit *x, size_t xn, Unit y)
{
	Unit r = 0;
	for (int i = (int)xn - 1; i >= 0; i--) {
		divUnit(&r, r, x[i], y);
	}
	return r;
}

//int (*mie::local::PrimitiveFunction::compare)(const Unit *x, size_t xn, const Unit *y, size_t yn) = &in_compare;
bool (*mie::local::PrimitiveFunction::addN)(Unit *out, const Unit *x, const Unit *y, size_t n) = &in_addN;
bool (*mie::local::PrimitiveFunction::add1)(Unit *out, const Unit *x, size_t n, Unit y) = &in_add;
bool (*mie::local::PrimitiveFunction::subN)(Unit *out, const Unit *x, const Unit *y, size_t n) = &in_subN;
bool (*mie::local::PrimitiveFunction::sub1)(Unit *out, const Unit *x, size_t n, Unit y) = &in_sub;
void (*mie::local::PrimitiveFunction::mul1)(Unit *out, const Unit *x, size_t n, Unit y) = &in_mul;
Unit (*mie::local::PrimitiveFunction::div1)(Unit *q, const Unit *x, size_t n, Unit y) = &in_div;
Unit (*mie::local::PrimitiveFunction::mod1)(const Unit *x, size_t n, Unit y) = &in_mod;

class Code : public Xbyak::CodeGenerator {
	void genAddSub(bool isAdd)
	{
		using namespace Xbyak;
		inLocalLabel();
#ifdef XBYAK32
		mov(ecx, ptr [esp + 16]); // n
		cmp(ecx, 1);
		jnz(".main");
		mov(ecx, ptr [esp + 8]); // x
		mov(eax, ptr [ecx]);
		mov(ecx, ptr [esp + 12]); // y
		if (isAdd) {
			add(eax, ptr [ecx]);
		} else {
			sub(eax, ptr [ecx]);
		}
		mov(ecx, ptr [esp + 4]); // out
		mov(ptr [ecx], eax);
		mov(eax, 0);
		setc(al);
		ret();
	L(".main");
		push(ebx);
		push(esi);
		push(edi);
		const int P = 4 * 3;
		mov(edx, ptr [esp + P + 4]); // out
		mov(esi, ptr [esp + P + 8]); // x
		mov(edi, ptr [esp + P + 12]); // y
		test(ecx, 1);
		jnz(".odd");
		shr(ecx, 1); // clear CF because ecx is even
		jmp(".lp");
	L(".odd");
		shr(ecx, 1);
		mov(eax, ptr [esi]);
		lea(esi, ptr [esi + 4]);
		if (isAdd) {
			add(eax, ptr [edi]);
		} else {
			sub(eax, ptr [edi]);
		}
		lea(edi, ptr [edi + 4]);
		mov(ptr [edx], eax);
		lea(edx, ptr [edx + 4]);
		align(16);
	L(".lp");
		mov(eax, ptr [esi]);
		mov(ebx, ptr [esi + 4]);
		if (isAdd) {
			adc(eax, ptr [edi]);
			adc(ebx, ptr [edi + 4]);
		} else {
			sbb(eax, ptr [edi]);
			sbb(ebx, ptr [edi + 4]);
		}
		lea(esi, ptr [esi + 8]);
		lea(edi, ptr [edi + 8]);
		mov(ptr [edx], eax);
		mov(ptr [edx + 4], ebx);
		lea(edx, ptr [edx + 8]);
		dec(ecx);
		jnz(".lp");
	L(".exit");
		mov(eax, 0); // don't clear CF
		setc(al);
		pop(edi);
		pop(esi);
		pop(ebx);
		ret();
#else
		const Reg64& a = rax;
#ifdef XBYAK64_WIN
		const Reg64& out = rcx;
		const Reg64& x = rdx;
		const Reg64& y = r8;
		const Reg64& n = r9;
		const Reg64& t0 = r10;
		const Reg64& t1 = r11;
		const Reg64& t2 = rsi;
#else
		const Reg64& out = rdi;
		const Reg64& x = rsi;
		const Reg64& y = rdx;
		const Reg64& n = rcx;
		const Reg64& t0 = r8;
		const Reg64& t1 = r9;
		const Reg64& t2 = r10;
#endif
		cmp(n, 4);
		jge(".main", T_NEAR);
		cmp(n, 1);
		jne("@f");
		// n == 1
		mov(a, ptr [x]);
		if (isAdd) {
			add(a, ptr [y]);
		} else {
			sub(a, ptr [y]);
		}
		mov(ptr [out], a);
		mov(a, 0);
		setc(al);
		ret();
	L("@@");
		cmp(n, 2);
		jne("@f");
		// n == 2
		mov(a , ptr [x + 8 * 0]);
		mov(t0, ptr [x + 8 * 1]);
		if (isAdd) {
			add(a , ptr [y + 8 * 0]);
			adc(t0, ptr [y + 8 * 1]);
		} else {
			sub(a , ptr [y + 8 * 0]);
			sbb(t0, ptr [y + 8 * 1]);
		}
		mov(ptr [out + 8 * 0], a);
		mov(ptr [out + 8 * 1], t0);
		mov(a, 0);
		setc(al);
		ret();
	L("@@");
		// n == 3
		mov(a , ptr [x + 8 * 0]);
		mov(t0, ptr [x + 8 * 1]);
		mov(t1, ptr [x + 8 * 2]);
		if (isAdd) {
			add(a , ptr [y + 8 * 0]);
			adc(t0, ptr [y + 8 * 1]);
			adc(t1, ptr [y + 8 * 2]);
		} else {
			sub(a , ptr [y + 8 * 0]);
			sbb(t0, ptr [y + 8 * 1]);
			sbb(t1, ptr [y + 8 * 2]);
		}
		mov(ptr [out + 8 * 0], a);
		mov(ptr [out + 8 * 1], t0);
		mov(ptr [out + 8 * 2], t1);
		mov(a, 0);
		setc(al);
		ret();
	L(".main"); // n >= 4
#ifdef XBYAK64_WIN
		mov(ptr [rsp + 8 * 1], t2);
#endif
		mov(a, n);
		shr(n, 2);
		and(a, 3);
		jz(".lp");
		cmp(a, 1);
		jne("@f");
		// 4x + 1
		mov(a, ptr [x + 8 * 0]);
		if (isAdd) {
			add(a, ptr [y + 8 * 0]);
		} else {
			sub(a, ptr [y + 8 * 0]);
		}
		mov(ptr [out + 8 * 0], a);
		lea(x, ptr [x + 8]);
		lea(y, ptr [y + 8]);
		lea(out, ptr [out + 8]);
		jmp(".lp");
	L("@@");
		cmp(a, 2);
		jne("@f");
		// 4x + 2
		mov(a , ptr [x + 8 * 0]);
		mov(t0, ptr [x + 8 * 1]);
		if (isAdd) {
			add(a , ptr [y + 8 * 0]);
			adc(t0, ptr [y + 8 * 1]);
		} else {
			sub(a , ptr [y + 8 * 0]);
			sbb(t0, ptr [y + 8 * 1]);
		}
		mov(ptr [out + 8 * 0], a);
		mov(ptr [out + 8 * 1], t0);
		lea(x, ptr [x + 8 * 2]);
		lea(y, ptr [y + 8 * 2]);
		lea(out, ptr [out + 8 * 2]);
		jmp(".lp");
	L("@@");
		// 4x + 3
		mov(a , ptr [x + 8 * 0]);
		mov(t0, ptr [x + 8 * 1]);
		mov(t1, ptr [x + 8 * 2]);
		if (isAdd) {
			add(a , ptr [y + 8 * 0]);
			adc(t0, ptr [y + 8 * 1]);
			adc(t1, ptr [y + 8 * 2]);
		} else {
			sub(a , ptr [y + 8 * 0]);
			sbb(t0, ptr [y + 8 * 1]);
			sbb(t1, ptr [y + 8 * 2]);
		}
		mov(ptr [out + 8 * 0], a);
		mov(ptr [out + 8 * 1], t0);
		mov(ptr [out + 8 * 2], t1);
		lea(x, ptr [x + 8 * 3]);
		lea(y, ptr [y + 8 * 3]);
		lea(out, ptr [out + 8 * 3]);
		align(16);
	L(".lp");
		mov(a , ptr [x + 8 * 0]);
		mov(t0, ptr [x + 8 * 1]);
		mov(t1, ptr [x + 8 * 2]);
		mov(t2, ptr [x + 8 * 3]);
		if (isAdd) {
			adc(a , ptr [y + 8 * 0]);
			adc(t0, ptr [y + 8 * 1]);
			adc(t1, ptr [y + 8 * 2]);
			adc(t2, ptr [y + 8 * 3]);
		} else {
			sbb(a , ptr [y + 8 * 0]);
			sbb(t0, ptr [y + 8 * 1]);
			sbb(t1, ptr [y + 8 * 2]);
			sbb(t2, ptr [y + 8 * 3]);
		}
		mov(ptr [out + 8 * 0], a);
		mov(ptr [out + 8 * 1], t0);
		mov(ptr [out + 8 * 2], t1);
		mov(ptr [out + 8 * 3], t2);
		lea(x, ptr [x + 8 * 4]);
		lea(y, ptr [y + 8 * 4]);
		lea(out, ptr [out + 8 * 4]);
		dec(n);
		jnz(".lp");
	L(".exit");
		mov(a, 0);
		setc(al);
#ifdef XBYAK64_WIN
		mov(t2, ptr [rsp + 8 * 1]);
#endif
		ret();
#endif
		outLocalLabel();
	}
#ifdef XBYAK64
	// add1(Unit *out, const Unit *x, size_t n, Unit y);
	void genAddSub1(bool isAdd)
	{
		using namespace Xbyak;
		inLocalLabel();
		const Reg64& a = rax;
		const Reg64& c = rcx;
#ifdef XBYAK64_WIN
		mov(r10, c);
		mov(c, r8); // n
		const Reg64& out = r10;
		const Reg64& x = rdx;
		const Reg64& y = r9;
		const Reg64& t = r11;
#else
		mov(r10, c);
		mov(c, rdx); // n
		const Reg64& out = rdi;
		const Reg64& x = rsi;
		const Reg64& y = r10;
		const Reg64& t = r8;
#endif
		lea(out, ptr [out + c * 8]);
		lea(x, ptr [x + c * 8]);
		xor(a, a);
		neg(c);
		mov(t, ptr [x + c * 8]);
		if (isAdd) {
			add(t, y);
		} else {
			sub(t, y);
		}
		mov(ptr [out + c * 8], t);
		inc(c);
#if 0
		jrcxz(".exit");
	L(".lp");
		mov(t, ptr [x + c * 8]);
		if (isAdd) {
			adc(t, a);
		} else {
			sbb(t, a);
		}
		mov(ptr [out + c * 8], t);
		inc(c);
		jrcxz(".exit");
		jmp(".lp");
#else
	// faster on Core i3
		jz(".exit");
	L(".lp");
		mov(t, ptr [x + c * 8]);
		if (isAdd) {
			adc(t, a);
		} else {
			sbb(t, a);
		}
		mov(ptr [out + c * 8], t);
		inc(c);
		jnz(".lp");
#endif

	L(".exit");
		setc(al);
		ret();
		outLocalLabel();
	}
#endif
	void genMul()
	{
		using namespace Xbyak;
		inLocalLabel();

		// void in_mul(Unit *out, const Unit *x, size_t n, Unit y)
#ifdef XBYAK32
		const Reg32& a = eax;
		const Reg32& d = edx;
		const Reg32& t = ebp;

		const Reg32& out = edi;
		const Reg32& x = esi;
		const Reg32& n = ecx;
		const Reg32& y = ebx;
		push(ebx);
		push(esi);
		push(edi);
		push(ebp);
		const int P = 4 * 4;
		mov(out, ptr [esp + P + 4]);
		mov(x, ptr [esp + P + 8]);
		mov(n, ptr [esp + P + 12]);
		mov(y, ptr [esp + P + 16]);

#else

		const Reg64& a = rax;
		const Reg64& d = rdx;
		const Reg64& t = r11;
		mov(r10, rdx);

#ifdef XBYAK64_WIN

		const Reg64& out = rcx;
		const Reg64& x = r10; // rdx
		const Reg64& n = r8;
		const Reg64& y = r9;
#else
		const Reg64& out = rdi;
		const Reg64& x = rsi;
		const Reg64& n = r10; // rdx
		const Reg64& y = rcx;
#endif
#endif
		const int s = (int)sizeof(Unit);
		xor(d, d);
	L(".lp");
		mov(t, d);
		mov(a, ptr [x]);
		mul(y); // [d:a] = [x] * y
		add(t, a);
		adc(d, 0);
		mov(ptr [out], t);
		add(x, s);
		add(out, s);
		sub(n, 1);
		jnz(".lp");
		mov(ptr [out], d);

#ifdef XBYAK32
		pop(ebp);
		pop(edi);
		pop(esi);
		pop(ebx);
#endif
		ret();
		outLocalLabel();
	}
	void genDiv()
	{
		using namespace Xbyak;
		inLocalLabel();

		// Unit in_div(Unit *q, const Unit *x, size_t xn, Unit y)
#ifdef XBYAK32
		const Reg32& a = eax;
		const Reg32& d = edx;

		const Reg32& q = edi;
		const Reg32& x = esi;
		const Reg32& n = ecx;
		const Reg32& y = ebx;
		push(ebx);
		push(esi);
		push(edi);
		const int P = 4 * 3;
		mov(q, ptr [esp + P + 4]);
		mov(x, ptr [esp + P + 8]);
		mov(n, ptr [esp + P + 12]);
		mov(y, ptr [esp + P + 16]);

#else

		const Reg64& a = rax;
		const Reg64& d = rdx;
		mov(r10, rdx);

#ifdef XBYAK64_WIN

		const Reg64& q = rcx;
		const Reg64& x = r10; // rdx
		const Reg64& n = r8;
		const Reg64& y = r9;
#else
		const Reg64& q = rdi;
		const Reg64& x = rsi;
		const Reg64& n = r10; // rdx
		const Reg64& y = rcx;
#endif
#endif
		const int s = (int)sizeof(Unit);
		lea(x, ptr [x + n * s - s]); // x = &x[xn - 1]
		lea(q, ptr [q + n * s - s]); // q = &q[xn - 1]
		xor(d, d); // r = 0
	L(".lp");
		mov(a, ptr [x]);
		div(y); // [d:a] / y = a ... d ; q = a, r = d
		mov(ptr [q], a);
		sub(x, s);
		sub(q, s);
		sub(n, 1);
		jnz(".lp");
		mov(a, d);

#ifdef XBYAK32
		pop(edi);
		pop(esi);
		pop(ebx);
#endif
		ret();
		outLocalLabel();
	}
	void genMod()
	{
		using namespace Xbyak;
		inLocalLabel();

		// Unit mod1(const Unit *x, size_t n, Unit y);
#ifdef XBYAK32
		const Reg32& a = eax;
		const Reg32& d = edx;

		const Reg32& x = esi;
		const Reg32& n = ecx;
		const Reg32& y = ebx;
		push(ebx);
		push(esi);
		const int P = 4 * 2;
		mov(x, ptr [esp + P + 4]);
		mov(n, ptr [esp + P + 8]);
		mov(y, ptr [esp + P + 12]);

#else

		const Reg64& a = rax;
		const Reg64& d = rdx;
		mov(r10, rdx);

#ifdef XBYAK64_WIN

		const Reg64& x = rcx;
		const Reg64& n = r10; // rdx
		const Reg64& y = r8;
#else
		const Reg64& x = rdi;
		const Reg64& n = rsi;
		const Reg64& y = r10; // rdx
#endif
#endif
		const int s = (int)sizeof(Unit);
		lea(x, ptr [x + n * s - s]); // x = &x[xn - 1]
		xor(d, d); // r = 0
	L(".lp");
		mov(a, ptr [x]);
		div(y); // [d:a] / y = a ... d ; q = a, r = d
		sub(x, s);
		sub(n, 1);
		jnz(".lp");
		mov(a, d);

#ifdef XBYAK32
		pop(esi);
		pop(ebx);
#endif
		ret();
		outLocalLabel();
	}
public:
	Code()
	{
		mie::local::PrimitiveFunction::addN = (bool (*)(Unit *, const Unit *, const Unit *, size_t))getCurr();
		genAddSub(true);
		align(16);
#ifdef XBYAK64
		mie::local::PrimitiveFunction::add1 = (bool (*)(Unit *, const Unit *, size_t, Unit))getCurr();
		genAddSub1(true);
		align(16);
#endif
		mie::local::PrimitiveFunction::subN = (bool (*)(Unit *, const Unit *, const Unit *, size_t))getCurr();
		genAddSub(false);
		align(16);
		mie::local::PrimitiveFunction::mul1 = (void (*)(Unit *, const Unit *, size_t, Unit))getCurr();
		genMul();
		align(16);
		mie::local::PrimitiveFunction::div1 = (Unit (*)(Unit *, const Unit *, size_t, Unit))getCurr();
		genDiv();
		align(16);
		mie::local::PrimitiveFunction::mod1 = (Unit (*)(const Unit *, size_t, Unit))getCurr();
		genMod();
	}
};

void mie::zmInit()
{
	static bool isInit = false;
	if (isInit) return;
	isInit = true;
	try {
		printf("zmInit ");
#ifdef XBYAK32
		puts("32bit");
#else
		puts("64bit");
#endif
		static Code code;
		return;
	} catch (std::exception& e) {
		fprintf(stderr, "zmInit ERR:%s\n", e.what());
	}
	::exit(1);
}

/*
  Local Variables:
  c-basic-offset: 4
  indent-tabs-mode: t
  tab-width: 4
  End:
*/
