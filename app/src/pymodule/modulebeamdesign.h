#ifndef MODULEBEAMDESIGN_H
#define MODULEBEAMDESIGN_H

#include "moduleheader.h"
#include<list>

// module beam_design from beam_design.cpp
PYBIND11_EMBEDDED_MODULE(cpp_beam_design, handle) {
	handle.def("beamDesign", &beam_design);
}

/////////////////////////////// begin ///////////////////////////////
// module data from data.py
std::tuple<double, double, double,
	double, double, double, double, double,
	double, double, double, double, double,
	double, double, double, double, double,
	double, double, double, double, double,
	double, double, double, double, double
> _data_analysis_() {
	auto pyData = py::module::import("scripts.data"); //import from data.py
	auto func = pyData.attr("data_analysis"); //import from python file -> data.py

	return func().cast< std::tuple<double, double, double,
		double, double, double, double, double,
		double, double, double, double, double,
		double, double, double, double, double,
		double, double, double, double, double,
		double, double, double, double, double
	> >();
}

// module data from function data_analysis()
PYBIND11_EMBEDDED_MODULE(module_data, handle) {
	handle.def("beamDataAnalysis", &_data_analysis_); //C++ export to python
}
// end of module data
/////////////////////////////// end ///////////////////////////////


/////////////////////////////// begin data from excel ///////////////////////////////
// module data_excel from data.py
std::list<std::tuple<std::string, double, double, double,
	double, double, double, double, double, double,
	double, double, double, double, double, double,
	double, double, double, double, double, 
	double, double, double, double, double,
	double, double, double, double, double> > 
	_data_design_excel_() {
	auto pyData = py::module::import("scripts.data"); //import from data.py
	auto func = pyData.attr("data_design_excel"); //import from python file -> data.py

/*	std::array<std::tuple< std::string, double, double, double,
		double, double, double, double, double, double,
		double, double, double, double, double, double,
		double, double, double, double, double,
		double, double, double, double, double,
		double, double, double, double, double>, 5> result;
	*/

	return func().cast< std::list<std::tuple< std::string, double, double, double,
		double, double, double, double, double, double,
		double, double, double, double, double, double,
		double, double, double, double, double,
		double, double, double, double, double,
		double, double, double, double, double> > >();
}

// module data_excel from function data_design_excel()
PYBIND11_EMBEDDED_MODULE(module_data_excel, handle) {
	handle.def("beamDataDesignExcel", &_data_design_excel_); //C++ export to python
}
// end of module data_excel
/////////////////////////////// end of data from excel ///////////////////////////////


void concrete_beam_design() {
	std::cout << "[C++] call function beam_design_summary from Python" << '\n';

	auto pyModule = py::module::import("scripts.main");
	auto _beam_design_ = pyModule.attr("beam_design_summary");
	//std::array<double, 5> func_array = _beam_design_().cast<std::array<double, 5>>();
	
	std::list<std::array<double, 5>> func_array = _beam_design_().cast<std::list<std::array<double, 5>> >();

/*	std::cout << "print element from C++" << '\n';
	std::cout << std::fixed << std::setprecision(2);
	for (const auto& element : func_array) {
		unsigned int j{};
		for (const auto& a : element)
		{
			std::cout << "element[" << j << "] : " << a << '\n';
			++j;
		}
		std::cout << typeid(element).name() << "\n";
	}
	std::cout << typeid(func_array).name() << "\n";
	std::cout << "end print element from C++" << '\n'; */
}

/*
void check_deflection() {
	auto pyMod = py::module::import("scripts.beam_deflection");
	auto beam_deflection = pyMod.attr("beamDeflection");
	beam_deflection();
}
*/

#endif
