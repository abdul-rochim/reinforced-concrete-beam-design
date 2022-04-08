#ifndef EXAMPLE_H
#define EXAMPLE_H

#include "moduleheader.h"

template<typename T>
struct Rect {
	Rect(T nWidth, T nLength) : width(nWidth), length(nLength) {}
	T getWidth() const { return width; }
	T getLength() const { return length; }
	T getArea() const { return width * length; }


private:
	T width; T length;
};

PYBIND11_EMBEDDED_MODULE(cpp_rect, handle) {
	py::class_<Rect<double>>(handle, "Rect")
		.def(py::init<double, double>())
		.def("get_width", &Rect<double>::getWidth)
		.def("get_length", &Rect<double>::getLength)
		.def("get_area", &Rect<double>::getArea);
}


template<typename T>
T multi(T n) {
	return n * n;
}

PYBIND11_EMBEDDED_MODULE(cpp_math, handle) {
	handle.def("multiplication", &multi<float>);
}

inline std::tuple<double, double, double> s_provided_;  //declaration of tuple
//inline double s_provided_{};  //declaration of variable
void foo() {
	auto pyModule = py::module::import("scripts.example_");
	auto fy_420 = pyModule.attr("fy420");
	auto fc_30 = pyModule.attr("fc30");

	float fy420 = fy_420().cast<float>();
	float fc30 = fc_30().cast<float>();

	std::cout << "fy420 : " << fy420 << '\n';
	std::cout << "fc30 : " << fc30 << '\n';

	auto exampleModule = py::module::import("scripts.example_");
	auto func1 = exampleModule.attr("sayHello");
	func1();

	auto func2 = exampleModule.attr("add");
	float result = func2(1.2, 2.2).cast<float>();
	std::cout << "[C++] call function add from Python\n" << result << '\n';

	auto func3 = exampleModule.attr("calculate");
	float res = func3(9.2).cast<float>();
	std::cout << "[C++] call function calculate from Python\n" << res << '\n';

	auto func4 = exampleModule.attr("plots");
	func4();
	std::cout << "[C++] call function plots from Python\n";

	auto get_point1 = exampleModule.attr("point1");
	auto get_point2 = exampleModule.attr("point2");
	//double p1 = get_point1().cast<double>();
	//double p2 = get_point2().cast<double>();

	auto func_area = exampleModule.attr("get_area_rect");
	Rect re(get_point1().cast<double>(), get_point2().cast<double>());
	double area_rect = func_area(re).cast<double>();
	std::cout << "[C++] call function area Rect from Python\n" << area_rect << " mm2\n";

	std::cout << ".......................................\n";
	auto get_array = exampleModule.attr("array__");
	auto v = get_array();// get_array().cast<std::array<std::array, 3>,2>(); // can not do this cause get_array() is pybind11::object
	for (auto elements : v) { // el.cast<std::array, 3>(); // can not do this cause elements is pybind11::handle of type array
		for (auto value : elements) {
			std::cout << value.cast<int>() << ' '; // can do this cause value is element of array
			std::cout << typeid(value).name() << ' ';
		}
		std::cout << '\n' << typeid(elements).name() << '\n';
		std::cout << "\n---------------------------------------\n";
	}
	std::cout << typeid(v).name() << std::endl;
	std::cout << ".......................................\n";

	std::cout << "+++++++++++++++++++++++++++++++++++++++\n";
	auto get_tup_array = exampleModule.attr("tuple_in_array");
	auto v2 = get_tup_array();
	unsigned int i{};
	for (auto elements : v2) {
		for (auto a : elements) {
			std::cout << "i : " << i << ' ';
			if (i == 1) {
				s_provided_ = std::make_tuple(a.cast<double>(), 0.0, 0.0);
				//s_provided_ = a.cast<double>();
			}
			std::cout << a.cast<int>() << ' ';
			//std::cout << typeid(a).name() << ' ';
			++i;
		}
		//std::cout << '\n' << typeid(elements).name() << '\n';
		//std::cout << "\n---------------------------------------\n";
	}
	//std::cout << typeid(v).name() << std::endl;
	std::cout << "\n+++++++++++++++++++++++++++++++++++++++\n";

}

#endif