
High-Speed Software Implementation of the Optimal Ate Pairing over BN Curves
=============

Abstruct
-------------

The library is able to compute the optimal ate pairing over a 254-bit prime field Fp.
p = 36t^4 + 36t^3 + 24t^2 + 6t + 1, t = -(2^62 + 2^55 + 1).

> old version:
>The library is able to compute the optimal ate pairing over a 254-bit prime field Fp.
>p = 36t^4 + 36t^3 + 24t^2 + 6t + 1, t = 2^62 - 2^54 + 2^44.

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

### Linux

Type

    cd ate && make

then you can get ate/test/bn .

>Remark:The other binaries except bn are for only test.

Benchmark
-------------

Core i7 2620 3.4GHz Windows 7 1.385M clock cycles(0.409msec).

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

* 2012/Jan/30 rewrite ate pairing according to "Faster explicit formulas for computing pairings over ordinary curves.  see http://www.patricklonga.bravehost.com/speed_pairing.html
* [2010/Sep/8 change xi from u + 12 to u
* 2010/Jul/15 use cyclotomic squaring for final_exp
* 2010/Jun/18 first release

AUTHOR
-------------

* MITSUNARI Shigeo(herumi at nifty dot com)
* TERUYA Tadanori

