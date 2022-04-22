//Reinforced Concrete Optimization
//Version	: 1.0.0

#include "pymodule/modulebeamdesign.h"
#include "pymodule/example.h"

//#include <limits>

int main()
{
	//std::cout << "[C++] Program started\n\n";

	py::scoped_interpreter guard{};
//	py::exec(R"(
//		print("[Python] called from C++")
//	)");
	try {
		//from python module
		concrete_beam_design();

		//check beam deflection
		check_deflection();

		std::cout << "\nthe analysis complete!\n";
	}
	catch (py::error_already_set& e) {
		std::cout << e.what() << '\n';
	}
	//std::cout << "[C++] Program finished" << '\n';

	return EXIT_SUCCESS;
}