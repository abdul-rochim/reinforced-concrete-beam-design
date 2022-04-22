//source: The Reinforced Concrete Design Handbook, SP-17M(14)

#ifndef MOMENT_HPP
#define MOMENT_HPP

#define _USE_MATH_DEFINES
#include <cmath>
//#include <utility>
#include <iostream>
#include <iomanip>
#include <fstream>

struct Moment {
	//constructors
	Moment() noexcept = default;
	~Moment() {} //destructor
	Moment(const double&, const double&, const double&, const double&, const double&) noexcept;
	Moment(const Moment&) noexcept ;	//copy constructor
	Moment(Moment&&) noexcept;		//move constructor
	Moment& operator = (const Moment&) noexcept;	//copy assignment
	Moment& operator = (Moment&&) noexcept;			//move assignment

	//functions
	const double initial_As_req();					//cross section area of longitudinal reinforcement
	//concrete compressive stress distribution
	//Code allows the use of an equivalent rectangular compressive stress distribution
	//of 0.85fc' with a depth of: a = betha1.c
	const double a(const double& As);				//depth of compressive concrete
	void swap(Moment&) noexcept;
	void print_As_required() { std::cout << initial_As_req() << '\n'; }

//private:
	//fields or members
	double Mu{};					//moment ultimate				
	double fy{};					//steel strength
	double fc{};					//concrete strength
	double b{};						//section width
	double d{};						//effective depth
	double betha1{};				//function of concrete compressive strength
	double Es = double(200000);		//young moduli of steel reinf.
	double Ec = double(25000);		//young moduli of concrete
	double eps_c = double(0.003);
};

inline Moment::Moment(const double& Mu, const double& fy, const double& fc, const double& b, const double& d) noexcept {
	this->Mu = Mu; //kN.m // this->Mu = Mu * 1000000.; //N.mm if multiply 1000000 -> do not use constructor delegation.
	this->fy = fy;
	this->fc = fc;
	this->b = b;
	this->d = d;

	if (fc >= double(17) && fc <= double(28))
		this->betha1 = double(0.85);
	else if (fc >= double(55))
		this->betha1 = double(0.65);
	else
		this->betha1 = double(0.85) - double(0.05) * (fc - double(28)) / double(7);

	this->Ec = double(4700) * std::sqrt(fc);
}

inline Moment::Moment(const Moment& m) noexcept : Moment(m.Mu, m.fy, m.fc, m.b, m.d)
{
/*	this->Mu = m.Mu; //N.mm
	this->fy = m.fy;
	this->fc = m.fc;
	this->b = m.b;
	this->d = m.d;
	this->betha1 = m.betha1;
	this->Ec = m.Ec;
*/
	//std::cout << "cpy constructor\n";
}

inline Moment::Moment(Moment&& moved) noexcept : //Moment(std::move(moved.fy), std::move(moved.fc), std::move(moved.As), std::move(moved.b), std::move(moved.d))
	Mu(std::exchange(moved.Mu, double(0))),
	fy(std::exchange(moved.fy, double(0))),
	fc(std::exchange(moved.fc, double(0))),
	b(std::exchange(moved.b, double(0))),
	d(std::exchange(moved.d, double(0))),
	betha1(std::exchange(moved.betha1, double(0))),
	Ec(std::exchange(moved.Ec, double(0)))
{
	//std::cout << "mv constructor\n";
}

inline Moment& Moment::operator = (const Moment& other) noexcept {
	//std::cout << "cpy assignment\n";
	if (this != &other) {
		Moment copy{ other };
		swap(copy);
	}
	return *this;
}

inline Moment& Moment::operator = (Moment&& other) noexcept {
	//std::cout << "mv assignment\n";
	if (this != &other) {
		Moment moved{ std::move(other) };
		swap(moved);
	}
	return *this;
}

//moment strength reduction factor
double phi_moment(Moment& m, const double& As_, double& eps_t_);

//moment nominal capacity
double phiMn(Moment& m, const double& dia_long, const double& d_hoop, const double& cover, double* As_provided, double* nb_bar);

#endif
