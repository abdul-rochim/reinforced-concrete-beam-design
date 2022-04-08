//source: The Reinforced Concrete Design Handbook, SP-17M(14)

#ifndef BEAM_DESIGN_HPP
#define BEAM_DESIGN_HPP

#include <array>
#include <typeinfo>

#include "material.hpp"
#include "moment.hpp"
#include "shear.hpp"
#include "torsion.hpp"

constexpr double PI = double(3.14159265358979323846);
//const double pi = std::acos(-1);

std::array<double, 5> beam_design(
	std::tuple<
	std::string,
	double,
	double,
	double,
	double,
	double,
	double,
	double,
	double,
	double,
	double,
	double,
	double,
	double*,
	double*,
	double*,
	double*,
	double,
	double,
	double,
	double,
	double,
	double,
	double*,
	double,
	double,
	double,
	double,
	double,
	double,
	double>);

#endif