
High-Speed Software Implementation of the Optimal Ate Pairing over Barreto-Naehrig Curves
=============

Abstruct
-------------

The library is able to compute the optimal ate pairing over a Barreto-Naerig curve defined over a 254-bit prime field Fp,
where p = 36z^4 + 36z^3 + 24z^2 + 6z + 1, z = -(2^62 + 2^55 + 1).

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

Our current implementation does not support the algoirthm in Pereira et al.
[A family of mplementation-friendly BN elliptic curves](http://www.sciencedirect.com/science/article/pii/S0164121211000914).

Benchmark
-------------

1.17M clock cycles on Core i7 4700MQ (Haswell) 2.4GHz processor with TurboBoost technology disable
(see http://eprint.iacr.org/2013/362).

Clock cycle counts of operations on Core i7 2600 3.4GHz/Xeon X5650 2.6GHz/Core i7 4700MQ 2.4GHz.

operation   | i7 2600|Xeon X5650|Haswell|Haswell with mulx
------------|--------|----------|-------|-----------------
TurboBoost  |on      |on        |off    |off
            |        |          |       |
mu          | 50     |60        |42     |38
r           | 80     |98        |69     |65
Fp:mul      |124     |146       |98     |90
Fp2:mul     |360     |412       |       |
Fp2:square  |288     |335       |       |
            |        |          |       |
G1::double  |1150    |1300      |       |
G1::add     |2200    |2600      |       |
G2::double  |2500    |2900      |       |
G2::add     |5650    |6500      |       |
Fp12::square|4500    |5150      |       |
Fp12::mul   |6150    |7000      |       |
            |        |          |       |
Miller loop |0.83M   |0.97M     |0.82M  |0.71M
final_exp   |0.53M   |0.63M     |0.51M  |0.46M
            |        |          |       |
pairing     |1.36M   |1.60M     |1.33M  |1.17M

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

### mulx supported with Haswell

test/bn.exe use mulx if possible and test/bn.exe -mulx 0 does not use mulx.
You can verify the difference with/without mulx on Haswell.

### Linux

Type:

    $ git clone git://github.com/herumi/xbyak.git
    $ git clone git://github.com/herumi/ate-pairing.git
    $ cd ate-pairing
    $ make -j
    $ test/bn

>Remark:The other binaries except bn are for only test.

Parameter
-------------
* elliptic curve : y^2 = x^3 + 2
* Fp : an finite field of characteristic p = 16798108731015832284940804142231733909889187121439069848933715426072753864723
* Fp2 : Fp[u] / (u^2 + 1)
* Fp6 = Fp2[v] / (v^3 - Xi), Xi = -u - 1
* Fp12 : Fp6[w] / (w^2 - v)
* Ec1 : E(Fp)[n]
* Ec2 : inverse image of E'(Fp^2)[n] under twisting iso E' to E.
* opt_atePairing : Ec2 x Ec1 to Fp12

How to use
-------------

see sample2() in https://github.com/herumi/ate-pairing/blob/master/test/sample.cpp

gmp
-------------
You can use mpz_class for scalar multiplication of points on the elliptic curves,
if MIE_ATE_USE_GMP is defined.

    using namespace bn;
    Param::init();
    const Ec2 g2(...);
    const Ec1 g1(...);
    mpz_class a("123456789");
    mpz_class b("98765432");
    Ec1 g1a = g1 * a;
    Ec2 g2b = g2 * b;
    Fp12 e;
    // calc e : G2 x G1 -> G3 pairing
    opt_atePairing(e, g2b, g1a);


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

* http://eprint.iacr.org/2010/354
* [High-speed software implementation of the optimal ate pairing over Barreto-Naehrig curves](http://dl.acm.org/citation.cfm?id=1948969)
* http://homepage1.nifty.com/herumi/crypt/ate-pairing.html

History
-------------
* 2013/Jun/02 support mulx of Haswell
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
