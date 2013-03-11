
High-Speed Software Implementation of the Optimal Ate Pairing over Barreto-Naehrig Curves
=============

Abstruct
-------------

The library is able to compute the optimal ate pairing over a Barreto-Naerig curve defined over a 254-bit prime field Fp,
where p = 36z^4 + 36z^3 + 24z^2 + 6z + 1, z = -(2^62 + 2^55 + 1).

Requirements
-------------

* OS : 64-bit Windows/64-bit Linux/Mac OS X
* CPU : x64 Intel/AMD processor
* C++ compiler : Visual Studio 2008, gcc 4.4.1 or later

Build
-------------

### Windows

Open ate/ate.sln and compile test_bn with Release mode,
then you can get the binary in ate/x64/Release/test_bn.exe .

### Windows with cygwin

* Install mingw64-x86_64-gcc-g++ (run cygwin setup and search mingw64)
* PATH=/usr/x86_64-w64-mingw32/sys-root/mingw/bin/:$PATH
* make -j
* test/bn.exe

### Linux

Type:

    $ make

or, if you have OMake:

    $ omake

then you can get ate/test/bn .

>Remark:The other binaries except bn are for only test.

Benchmark
-------------

1.35M clock cycles(0.399msec) on Core i7 2620 3.4GHz Windows 7 with turbo boost.

Operation costs
-------------

We compare our library with Aranha et al.(http://eprint.iacr.org/2010/526).

* mu ; unreduced multiplication producing double-precision result(256-bit int x 256-bit int to 512-bit int).
* r : modular reduction of double-precision integers(512-bit int to 256-bit int in Fp).
* Fp:mul = mu + r
* Fp2:mul = 3mu + 2r
* Fp2:square = 2mu + 2r

Phase               | Aranha et al. | This work
--------------------|---------------|---------------
Miller Loop         | 6792mu + 3022r| 6785mu + 3022r
Final Exponentiation| 3753mu + 2006r| 3526mu + 1932r
Optimal Ate Pairing |10545mu + 5028r|10311mu + 4954r

Remark : Their Table 2 in p.17 does not contain the cost of (m, r) so
I add the costs of (282m + 6mu + 4r), (30m + 75mu + 50r) for ML and FE respectively.

How to use
-------------

see sample2() in https://github.com/herumi/ate-pairing/blob/master/test/sample.cpp

Class
-------------
* elliptic curve : y^2 = x^3 + 2
* Fp : an finite field of characteristic p = 16798108731015832284940804142231733909889187121439069848933715426072753864723
* Fp2 : Fp[u] / (u^2 + 1)
* Fp6 = Fp2[v] / (v^3 - Xi), Xi = -u - 1
* Fp12 : Fp6[w] / (w^2 - v)
* Ec1 : E(Fp)[n]
* Ec2 : inverse image of E'(Fp^2)[n] under twisting iso E' to E.
* opt_atePairing : Ec2 x Ec1 to Fp12

gmp
-------------
You can use mpz_class for scalar multiplication of points on the elliptic curves,
if MIE_ATE_USE_GMP is defined.

Xbyak
-------------

Xbyak is a x86/x86-64 JIT assembler for C++.
I made this library for developping pairing functions efficiently.

>http://homepage1.nifty.com/herumi/soft/xbyak_e.html

>https://github.com/heurmi/xbyak/

License
-------------

modified new BSD License

>http://opensource.org/licenses/BSD-3-Clause

If you have any questions or problems, just let me know, then I'm happy.

Reference
-------------

http://homepage1.nifty.com/herumi/crypt/ate-pairing.html

History
-------------
* 2013/Mar/08 add elliptic curve class
* 2012/Jan/30 rewrite ate pairing according to
  "Faster explicit formulas for computing pairings over ordinary curves.
  see http://www.patricklonga.bravehost.com/speed_pairing.html
* [2010/Sep/8 change xi from u + 12 to u
* 2010/Jul/15 use cyclotomic squaring for final_exp
* 2010/Jun/18 first release

AUTHOR
-------------

* MITSUNARI Shigeo(herumi@nifty.com)
* TERUYA Tadanori(tadanori.teruya@gmail.com)
