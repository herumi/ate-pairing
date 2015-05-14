/*
	a tiny sample of optimal ate pairing
*/
#ifdef __linux__
	#define MIE_ATE_USE_GMP
#endif
#include "bn_pairing.hpp"
#include "test_point.hpp"

static int errNum = 0;

template<class T, class S>
void verify(const char *msg, const T& a, const S& b)
{
	if (a == b) {
		printf("%s : ok\n", msg);
	} else {
		printf("%s : ng\n", msg);
		PUT(a);
		PUT(b);
		errNum++;
	}
}

void sample2(const bn::CurveParam& cp)
{
	using namespace bn;
	// init my library
	Param::init(cp);
	const Point& pt = selectPoint(cp);
	const Ec2 g2(
		Fp2(Fp(pt.g2.aa), Fp(pt.g2.ab)),
		Fp2(Fp(pt.g2.ba), Fp(pt.g2.bb))
	);
	const Ec1 g1(pt.g1.a, pt.g1.b);
	// verify g2 and g1 on curve
	verify("g1 is on EC", g1.isValid(), true);
	verify("g2 is on twist EC", g2.isValid(), true);
	puts("order of group");
	PUT(Param::r);
	{
		Ec1 t = g1 * Param::r;
		// Ec1::mul(t, g1, Param::r);
		verify("orgder of g1 == r", t.isZero(), true);
	}
	{
		Ec2 t = g2 * Param::r;
		// Ec2::mul(t, g2, Param::r);
		verify("order of g2 == r", t.isZero(), true);
	}
	const char *aStr = "123456789012345";
	const char *bStr = "998752342342342342424242421";
	const mie::Vuint a(aStr);
	const mie::Vuint b(bStr);

	// scalar-multiplication sample
	{
		const mie::Vuint c = a + b;
		Ec1 Pa = g1 * a;
		Ec1 Pb = g1 * b;
		Ec1 Pc = g1 * c;
		Ec1 out = Pa + Pb;

		verify("check g1 * c = g1 * a + g1 * b", Pc, out);
#ifdef MIE_ATE_USE_GMP
		{
			mpz_class aa(aStr);
			mpz_class bb(bStr);
			mpz_class cc = aa + bb;
			Ec1 Paa = g1 * aa;
			Ec1 Pbb = g1 * bb;
			Ec1 Pcc = g1 * cc;
			verify("gmp Paa == Pa", Paa, Pa);
			verify("gmp Pbb == Pb", Pbb, Pb);
			verify("gmp Pcc == Pc", Pcc, Pc);
		}
#endif
	}

	Fp12 e;
	// calc e : G2 x G1 -> G3 pairing
	opt_atePairing(e, g2, g1); // e = e(g2, g1)
	PUT(e);
	{
		Fp12 t = power(e, Param::r);
		verify("order of e == r", t, 1);
	}
	Ec2 g2a = g2 * a;
	// Ec2::mul(g2a, g2, a);
	Fp12 ea1;
	opt_atePairing(ea1, g2a, g1); // ea1 = e(g2a, g1)
	Fp12 ea2 = power(e, a); // ea2 = e^a
	verify("e(g2 * a, g1) = e(g2, g1)^a", ea1, ea2);

	Ec1 g1b = g1 * b;
	// Ec1::mul(g1b, g1, b);
	Fp12 eb1;
	opt_atePairing(eb1, g2, g1b); // eb1 = e(g2, g1b)
	Fp12 eb2 = power(e, b); // eb2 = e^b
	verify("e(g2a, g1 * b) = e(g2, g1)^b", eb1, eb2);

	Ec1 q1 = g1 * 12345;
	verify("q1 is on EC", q1.isValid(), true);
	Fp12 e1, e2;
	opt_atePairing(e1, g2, g1); // e1 = e(g2, g1)
	opt_atePairing(e2, g2, q1); // e2 = e(g2, q1)
	Ec1 q2 = g1 + q1;
	opt_atePairing(e, g2, q2); // e = e(g2, q2)
	verify("e = e1 * e2", e, e1 * e2);

	/*
		reduce one copy as the following
	*/
	Ec2::mul(g2a, g2, a); // g2a = g2 * a
	Ec1::mul(g1b, g1, b);
	verify("g2a == g2 * a", g2a, g2 * a);
	verify("g1b == g1 * b", g1b, g1 * b);
}


int main()
{
	bn::CurveParam cp = bn::CurveFp254BNb;
	sample2(cp);
	printf("errNum = %d\n", errNum);
}
