/*
	a tiny sample of optimal ate pairing
*/
#ifdef __linux__
	#define MIE_ATE_USE_GMP
#endif
#include "bn.h"

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

void sample1()
{
	using namespace bn;
	// init my library
	Param::init();
	// prepair a generator
	const Fp2 g2[3] = {
		Fp2(
			Fp("12723517038133731887338407189719511622662176727675373276651903807414909099441"),
			Fp("4168783608814932154536427934509895782246573715297911553964171371032945126671")
		),
		Fp2(
			Fp("13891744915211034074451795021214165905772212241412891944830863846330766296736"),
			Fp("7937318970632701341203597196594272556916396164729705624521405069090520231616")
		),
		Fp2(
			Fp("1"),
			Fp("0")
		)
	};
	const Fp g1[3] = {
		Fp("1674578968009266105367653690721407808692458796109485353026408377634195183292"),
		Fp("8299158460239932124995104248858950945965255982743525836869552923398581964065"),
		Fp("1")
	};
	// verify g2 and g1 on curve
	verify("g1 is on EC", ecop::isOnECJac3(g1), true);
	verify("g2 is on twist EC", ecop::isOnTwistECJac3(g2), true);
	puts("order of group");
	PUT(Param::r);
	{
		Fp t[3];
		ecop::ScalarMult(t, g1, Param::r);
		// (x, y, 0) means 0 at Jacobi coordinate
		verify("orgder of g1 == r", t[2], 0);
	}
	{
		Fp2 t[3];
		ecop::ScalarMult(t, g2, Param::r);
		verify("order of g2 == r", t[2], 0);
	}
	const char *aStr = "123456789012345";
	const char *bStr = "998752342342342342424242421";
	const mie::Vuint a(aStr);
	const mie::Vuint b(bStr);

	// scalar-multiplication sample
	{
		Fp Pa[3];
		Fp Pb[3];
		Fp Pc[3];
		Fp out[3];
		const mie::Vuint c = a + b;

		ecop::ScalarMult(Pa, g1, a); // Pa = g1 * a
		ecop::ScalarMult(Pb, g1, b); // Pb = g1 * b
		ecop::ScalarMult(Pc, g1, c); // Pc = g1 * (a + b)
		ecop::ECAdd(out, Pa, Pb); // g1 * a + g1 * b
		ecop::NormalizeJac(Pc, Pc);
		ecop::NormalizeJac(out, out);
		std::cout << std::hex;
		verify("check g1 * c = g1 * a + g1 * b", Pc[0] == out[0] && Pc[1] == out[1] && Pc[2] == out[2], true);

#ifdef MIE_ATE_USE_GMP
		{
			mpz_class aa(aStr);
			mpz_class bb(bStr);
			mpz_class cc = aa + bb;
			Fp Paa[3];
			Fp Pbb[3];
			Fp Pcc[3];
			ecop::ScalarMult(Paa, g1, aa); // Pa = g1 * a
			ecop::ScalarMult(Pbb, g1, bb); // Pb = g1 * b
			ecop::ScalarMult(Pcc, g1, cc); // Pc = g1 * (a + b)
			ecop::NormalizeJac(Pcc, Pcc);
			verify("gmp Paa == Pa", Paa[0] == Pa[0] && Paa[1] == Pa[1] && Paa[2] == Pa[2], true);
			verify("gmp Pbb == Pb", Pbb[0] == Pb[0] && Pbb[1] == Pb[1] && Pbb[2] == Pb[2], true);
			verify("gmp Pcc == Pc", Pcc[0] == Pc[0] && Pcc[1] == Pc[1] && Pcc[2] == Pc[2], true);
		}
#endif
	}

	Fp12 e;
	// calc e : G2 x G1 -> G3 pairing
	opt_atePairingJac<Fp>(e, g2, g1); // e = e(g2, g1)
	PUT(e);
	{
		Fp12 t = power(e, Param::r);
		verify("order of e == r", t, 1);
	}
	Fp2 g2a[3];
	ecop::ScalarMult(g2a, g2, a); // g2a = g2 * a
	Fp12 ea1;
	opt_atePairingJac<Fp>(ea1, g2a, g1); // ea1 = e(g2a, g1)
	Fp12 ea2 = power(e, a); // ea2 = e^a
	verify("e(g2 * a, g1) = e(g2, g1)^a", ea1, ea2);

	Fp g1b[3];
	ecop::ScalarMult(g1b, g1, b); // g1b = g1 * b
	Fp12 eb1;
	opt_atePairingJac<Fp>(eb1, g2, g1b); // eb1 = e(g2, g1b)
	Fp12 eb2 = power(e, b); // eb2 = e^b
	verify("e(g2a, g1 * b) = e(g2, g1)^b", eb1, eb2);

	const Fp q1[3] = {
		Fp("2"),
		Fp("16740896641879863340107777353588575149660814923656713498672603551465628253431"),
		Fp("1")
	};
	verify("q1 is on EC", ecop::isOnECJac3(q1), true);
	Fp12 e1, e2;
	opt_atePairingJac<Fp>(e1, g2, g1); // e1 = e(g2, g1)
	opt_atePairingJac<Fp>(e2, g2, q1); // e2 = e(g2, q1)
	Fp q2[3];
	ecop::ECAdd(q2, g1, q1); // q2 = g1 + q1
	opt_atePairingJac<Fp>(e, g2, q2); // e = e(g2, q2)
	verify("e = e1 * e2", e, e1 * e2);
}

void sample2()
{
	using namespace bn;
	// init my library
	Param::init();
	const Ec2 g2(
		Fp2(
			Fp("12723517038133731887338407189719511622662176727675373276651903807414909099441"),
			Fp("4168783608814932154536427934509895782246573715297911553964171371032945126671")
		),
		Fp2(
			Fp("13891744915211034074451795021214165905772212241412891944830863846330766296736"),
			Fp("7937318970632701341203597196594272556916396164729705624521405069090520231616")
		)
	);
	const Ec1 g1(
		Fp("1674578968009266105367653690721407808692458796109485353026408377634195183292"),
		Fp("8299158460239932124995104248858950945965255982743525836869552923398581964065")
	);
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

	const Ec1 q1(
		Fp("2"),
		Fp("16740896641879863340107777353588575149660814923656713498672603551465628253431")
	);
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
	puts("sample1");
	sample1();
	puts("sample2");
	sample2();
	printf("errNum = %d\n", errNum);
}
