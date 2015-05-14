#pragma once
/**
	@file
	@brief BN parameter
	@author herumi and t_teruya
	@note modified new BSD license
	http://opensource.org/licenses/BSD-3-Clause
*/
#define MIE_ATE_USE_GMP
#include "bn.h"

namespace mcl { namespace bn {

inline void SystemInit() throw(std::exception)
{
	::bn::Param::init();
}

class Fp2;
class Fp12;
class Ec1;
class Ec2;

class Mpz {
	mpz_class self_;
	friend class Fp12;
	friend class Ec1;
	friend class Ec2;
public:
	Mpz() {}
	Mpz(int x) throw(std::exception) : self_(x) {}
	Mpz(const std::string& str) throw(std::exception)
	{
		fromStr(str);
	}
	std::string toStr() const throw(std::exception)
	{
		return self_.get_str();
	}
	void fromStr(const std::string& str) throw(std::exception)
	{
		self_.set_str(str, 0);
	}
	std::string toString() const throw(std::exception) { return toStr(); }
	bool equals(const Mpz& rhs) const { return self_ == rhs.self_; }
	int compareTo(const Mpz& rhs) const { return mpz_cmp(self_.get_mpz_t(), rhs.self_.get_mpz_t()); }
	void add(const Mpz& rhs) throw(std::exception) { self_ += rhs.self_; }
	void sub(const Mpz& rhs) throw(std::exception) { self_ -= rhs.self_; }
	void mul(const Mpz& rhs) throw(std::exception) { self_ *= rhs.self_; }
	void mod(const Mpz& rhs) throw(std::exception) { self_ %= rhs.self_; }
};

class Fp {
	::bn::Fp self_;
	friend class Fp2;
	friend class Ec1;
public:
	Fp() {}
	Fp(int x) : self_(x) {}
	Fp(const std::string& str) throw(std::exception)
	{
		self_.set(str);
	}
	std::string toStr() const throw(std::exception)
	{
		return self_.toString();
	}
	void fromStr(const std::string& str) throw(std::exception)
	{
		self_.set(str);
	}
	std::string toString() const throw(std::exception) { return toStr(); }
	bool equals(const Fp& rhs) const { return self_ == rhs.self_; }
	void add(const Fp& rhs) throw(std::exception) { self_ += rhs.self_; }
	void sub(const Fp& rhs) throw(std::exception) { self_ -= rhs.self_; }
	void mul(const Fp& rhs) throw(std::exception) { self_ *= rhs.self_; }
};

class Fp2 {
	::bn::Fp2 self_;
	friend class Ec2;
public:
	Fp2() {}
	Fp2(int x) : self_(x) {}
	Fp2(const Fp& a, Fp& b) throw(std::exception)
		: self_(a.self_, b.self_)
	{
	}
	std::string toStr() const throw(std::exception)
	{
		return self_.toString();
	}
	void fromStr(const std::string& str) throw(std::exception)
	{
		self_.set(str);
	}
	std::string toString() const throw(std::exception) { return toStr(); }
	bool equals(const Fp2& rhs) const { return self_ == rhs.self_; }
	void add(const Fp2& rhs) throw(std::exception) { self_ += rhs.self_; }
	void sub(const Fp2& rhs) throw(std::exception) { self_ -= rhs.self_; }
	void mul(const Fp2& rhs) throw(std::exception) { self_ *= rhs.self_; }
};

class Fp12 {
	::bn::Fp12 self_;
public:
	Fp12() {}
	Fp12(int x) : self_(x) {}
	std::string toStr() const throw(std::exception)
	{
		std::ostringstream oss;
		oss << self_;
		return oss.str();
	}
	void fromStr(const std::string& str) throw(std::exception)
	{
		std::istringstream iss(str);
		iss >> self_;
	}
	std::string toString() const throw(std::exception) { return toStr(); }
	bool equals(const Fp12& rhs) const { return self_ == rhs.self_; }
	void add(const Fp12& rhs) throw(std::exception) { self_ += rhs.self_; }
	void sub(const Fp12& rhs) throw(std::exception) { self_ -= rhs.self_; }
	void mul(const Fp12& rhs) throw(std::exception) { self_ *= rhs.self_; }
	void pairing(const Ec2& ec2, const Ec1& ec1);
	void power(const Mpz& x)
	{
		mie::power(self_, x.self_);
	}
};

class Ec1 {
	::bn::Ec1 self_;
	friend class Fp12;
public:
	Ec1() {}
	Ec1(const Fp& x, const Fp& y) throw(std::exception)
	{
		set(x, y);
	}
	Ec1(const Fp& x, const Fp& y, const Fp& z) throw(std::exception)
	{
		set(x, y, z);
	}
	bool isValid() const { return self_.isValid(); }
	void set(const Fp& x, const Fp& y) throw(std::exception)
	{
		self_.set(x.self_, y.self_);
	}
	void set(const Fp& x, const Fp& y, const Fp& z) throw(std::exception)
	{
		self_.set(x.self_, y.self_, z.self_);
	}
	std::string toStr() const throw(std::exception)
	{
		std::ostringstream oss;
		oss << self_;
		return oss.str();
	}
	void fromStr(const std::string& str) throw(std::exception)
	{
		std::istringstream iss(str);
		iss >> self_;
	}
	std::string toString() const throw(std::exception) { return toStr(); }
	bool equals(const Ec1& rhs) const { return self_ == rhs.self_; }
	bool isZero() const { return self_.isZero(); }
	void clear() { self_.clear(); }
	void dbl() { ::bn::Ec1::dbl(self_, self_); }
	void neg() { ::bn::Ec1::neg(self_, self_); }
	void add(const Ec1& rhs) { ::bn::Ec1::add(self_, self_, rhs.self_); }
	void sub(const Ec1& rhs) { ::bn::Ec1::sub(self_, self_, rhs.self_); }
	void mul(const Mpz& rhs) { ::bn::Ec1::mul(self_, self_, rhs.self_); }
};

class Ec2 {
	::bn::Ec2 self_;
	friend class Fp12;
public:
	Ec2() {}
	Ec2(const Fp2& x, const Fp2& y) throw(std::exception)
	{
		set(x, y);
	}
	Ec2(const Fp2& x, const Fp2& y, const Fp2& z) throw(std::exception)
	{
		set(x, y, z);
	}
	bool isValid() const { return self_.isValid(); }
	void set(const Fp2& x, const Fp2& y) throw(std::exception)
	{
		self_.set(x.self_, y.self_);
	}
	void set(const Fp2& x, const Fp2& y, const Fp2& z) throw(std::exception)
	{
		self_.set(x.self_, y.self_, z.self_);
	}
	std::string toStr() const throw(std::exception)
	{
		std::ostringstream oss;
		oss << self_;
		return oss.str();
	}
	void fromStr(const std::string& str) throw(std::exception)
	{
		std::istringstream iss(str);
		iss >> self_;
	}
	std::string toString() const throw(std::exception) { return toStr(); }
	bool equals(const Ec2& rhs) const { return self_ == rhs.self_; }
	bool isZero() const { return self_.isZero(); }
	void clear() { self_.clear(); }
	void dbl() { ::bn::Ec2::dbl(self_, self_); }
	void neg() { ::bn::Ec2::neg(self_, self_); }
	void add(const Ec2& rhs) { ::bn::Ec2::add(self_, self_, rhs.self_); }
	void sub(const Ec2& rhs) { ::bn::Ec2::sub(self_, self_, rhs.self_); }
	void mul(const Mpz& rhs) { ::bn::Ec2::mul(self_, self_, rhs.self_); }
};

void Fp12::pairing(const Ec2& ec2, const Ec1& ec1)
{
	::bn::opt_atePairing(self_, ec2.self_, ec1.self_);
}

} } // mcl::bn

