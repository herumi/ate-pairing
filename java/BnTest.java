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
/*
	public static void assertEquals(String msg, int lhs, int rhs) {
		if (lhs == rhs) {
			System.out.println("OK : " + msg);
		} else {
			System.out.println("NG : " + msg + ", lhs = " + lhs + ", rhs = " + rhs);
		}
	}
*/
	public static void main(String argv[]) {
		try {
//			Mpz i1 = new Mpz(1);
//			Mpz i3 = new Mpz(3);
//			assertEquals("pub.add", i1, i3);

		} catch (RuntimeException e) {
			System.out.println("unknown exception :" + e);
		}
	}
}
