# Java class files(under construction)

## Build

* Install [swig](http://www.swig.org/) and Java.

### Windows

* set SWIG to the path to swig in make_wrap.bat
* set JAVA_DIR to the path to java in set-java-path.bat.
* Use the follogin commands:

    > cd java
    > make_wrap.bat

* bin/bn254_if_wrap.dll is a dll for java.

### Linux

* set JAVA_INC to the path to Java in Makefile

    > make test

* bin/libbn254_if_wrap.so is a shared library for java.

## API and Class

### Setup

* At first, call these functions.

    > System.loadLibrary("bn254_if_wrap");
	> BN254.SystemInit();

### class Mpz

* a wrapped class of mpz_class of GMP.

* Mpz(int), Mpz(String)
* set(int), set(String)
* add(Mpz x)
    * Set (this + x) to this.
* sub(Mpz x)
    * Set (this - x) to this.
* mul(Mpz x)
    * Set (this * x) to this.
* mod(Mpz x)
    * Set (this % x) to this.

### class Fp, Fp2, Fp12

* a wrapped class of bn::Fp, bn::Fp2, bn::Fp12.

** common method

* Fp(int), Fp(String)
* set(int), set(String)
* add(Fp x)
    * Set (this + x) to this.
* sub(Fp x)
    * Set (this - x) to this.
* mul(Fp x)
    * Set (this * x) to this.
* power(Mpz x)
    * Set (this ^ x) to this.

** Fp
** Fp2
** Fp12
* pairing(Ec2 ec2, Ec1 ec1)
    * Set opt_ate_pairing(ec2, ec1) to this.

* Ec1, Ec2

* a wrapped class of bn::Ec1 and bn::Ec2.


** Ec1(Fp x, Fp y)
    * set (x, y, 1) to this.
** Ec1(Fp x, Fp y, Fp z)
    * set (x, y, z) to this.


