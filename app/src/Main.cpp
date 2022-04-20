//REINFORCED CONCRETE DESIGN - RECTANGULAR SECTION BEAM
//Version 1.0.0
//email : abdul.rochim.civeng@gmail.com

#include "pymodule/modulebeamdesign.h"
#include "pymodule/example.h"

//#include <limits>

int main()
{
	std::cout << std::fixed << std::setprecision(2);
	std::cout << "[C++] Program started" << '\n';
	
	py::scoped_interpreter guard{};
//	py::exec(R"(
//		print("[Python] Python says Hello")
//	)");
	try {
		//from python module
		concrete_beam_design();

	}
	catch (py::error_already_set& e) {
		std::cout << e.what() << '\n';
	}

	std::cout << "\n\n";
	//std::cin.get();
	return EXIT_SUCCESS;
	//return 0;
}
