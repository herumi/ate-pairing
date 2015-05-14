import java.io.*;
import mcl.bn.*;

public class BnTest {
	static {
		System.loadLibrary("bn_if_wrap");
	}
	public static void assertEquals(String msg, Mpz lhs, Mpz rhs) {
		if (lhs.equals(rhs)) {
			System.out.println("OK : " + msg);
		} else {
			System.out.println("NG : " + msg + ", lhs = " + lhs + ", rhs = " + rhs);
		}
	}
	public static void main(String argv[]) {
		try {
//			Fp aa = new Fp("12723517038133731887338407189719511622662176727675373276651903807414909099441");
//			Fp ab = new Fp("4168783608814932154536427934509895782246573715297911553964171371032945126671");
//			Fp ba = new Fp("13891744915211034074451795021214165905772212241412891944830863846330766296736");
//			Fp bb = new Fp("7937318970632701341203597196594272556916396164729705624521405069090520231616");
			Bn.SystemInit();
			Ec1 g1 = new Ec1(new Fp(-1), new Fp(1));
			System.out.println("g1=" + g1);
//			Ec2 g2 = new Ec2(new Fp2(aa, ab), new Fp2(ba, bb));
		} catch (RuntimeException e) {
			System.out.println("unknown exception :" + e);
		}
	}
}
