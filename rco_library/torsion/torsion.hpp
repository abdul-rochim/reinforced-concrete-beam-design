//source: The Reinforced Concrete Design Handbook, SP-17M(14)

#ifndef TORSION_HPP
#define TORSION_HPP

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>
#include <tuple>
#include <fstream>

std::tuple<bool, double, double, double> can_torsion_be_neglected(
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&);

struct Torsion {
	double Tth;
	double Tcr;
	double At_per_s_req;
	double Al_req;
};

Torsion torsion_design(
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&);

#endif
