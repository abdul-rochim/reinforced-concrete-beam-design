//source: The Reinforced Concrete Design Handbook, SP-17M(14)

#ifndef SHEAR_HPP
#define SHEAR_HPP

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>

struct Shear {
	double Vc, Vs_prov, Av_per_s_required, Av_per_s_prov;
};

Shear shear_design(
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	const double&, 
	double&,
	const double&);

#endif