/*
	a tiny sample of optimal ate pairing
*/
#include "bn.h"

int main()
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
	printf("g1 is on EC : %s\n", ecop::isOnECJac3(g1) ? "ok" : "ng");
	printf("g2 is on twist EC : %s\n", ecop::isOnTwistECJac3(g2) ? "ok" : "ng");
	puts("order of group");
	PUT(Param::r);
	{
		Fp t[3];
		ecop::ScalarMult(t, g1, Param::r);
		printf("order of g1 == r :%s\n", t[2] == 0 ? "ok" : "ng"); // (x, y, 0) means 0 at Jacobi coordinate
	}
	{
		Fp2 t[3];
		ecop::ScalarMult(t, g2, Param::r);
		printf("order of g2 == r :%s\n", t[2] == 0 ? "ok" : "ng"); // (x, y, 0) means 0 at Jacobi coordinate
	}
	const mie::Vuint a("123456789012345");
	const mie::Vuint b("998752342342342342424242421");

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
		puts("check g1 * c = g1 * a + g1 * b");
		PUT(Pc[0] - out[0]);
		PUT(Pc[1] - out[1]);
		PUT(Pc[2] - out[2]);
	}

	Fp12 e;
	// calc e : G2 x G1 -> G3 pairing
	opt_atePairingJac<Fp>(e, g2, g1); // e = e(g2, g1)
	PUT(e);
	{
		Fp12 t = power(e, Param::r);
		printf("order of e == r :%s\n", t == 1 ? "ok" : "ng");
	}
	Fp2 g2a[3];
	ecop::ScalarMult(g2a, g2, a); // g2a = g2 * a
	Fp12 ea1;
	opt_atePairingJac<Fp>(ea1, g2a, g1); // ea1 = e(g2a, g1)
//	PUT(ea1);
	Fp12 ea2 = power(e, a); // ea2 = e^a
//	PUT(ea2);
	printf("verify e(g2 * a, g1) = e(g2, g1)^a : %s\n", ea1 == ea2 ? "ok" : "ng");

	Fp g1b[3];
	ecop::ScalarMult(g1b, g1, b); // g1b = g1 * b
	Fp12 eb1;
	opt_atePairingJac<Fp>(eb1, g2, g1b); // eb1 = e(g2, g1b)
//	PUT(eb1);
	Fp12 eb2 = power(e, b); // eb2 = e^b
//	PUT(eb2);
	printf("verify e(g2a, g1 * b) = e(g2, g1)^b : %s\n", eb1 == eb2 ? "ok" : "ng");

	const Fp q1[3] = {
		Fp("2"),
		Fp("16740896641879863340107777353588575149660814923656713498672603551465628253431"),
		Fp("1")
	};
	printf("q1 is on EC : %s\n", ecop::isOnECJac3(q1) ? "ok" : "ng");
	Fp12 e1, e2;
	opt_atePairingJac<Fp>(e1, g2, g1); // e1 = e(g2, g1)
	opt_atePairingJac<Fp>(e2, g2, q1); // e2 = e(g2, q1)
	Fp q2[3];
	ecop::ECAdd(q2, g1, q1); // q2 = g1 + q1
	opt_atePairingJac<Fp>(e, g2, q2); // e = e(g2, q2)
	printf("verify e = e1 * e2 : %s\n", e == e1 * e2 ? "ok" : "ng");
}
