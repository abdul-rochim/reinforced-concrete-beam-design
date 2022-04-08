
cmake_minimum_required(VERSION 3.10...3.22)

project("RCO_Libs" VERSION 1.0 DESCRIPTION "the optimization of reinforced concrete")

add_library(rco_lib_obj OBJECT)
#add_library(rco_lib_obj SHARED "material/material.cpp" "moment/moment.cpp")

target_sources(rco_lib_obj
	PRIVATE
		"material/material.cpp"
		"moment/moment.cpp"
		"shear/shear.cpp"
		"torsion/torsion.cpp"
		"reinforced_concrete/beam_design.cpp"

	PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}/material/material.hpp"
		"${CMAKE_CURRENT_SOURCE_DIR}/moment/moment.hpp"
		"${CMAKE_CURRENT_SOURCE_DIR}/shear/shear.hpp"
		"${CMAKE_CURRENT_SOURCE_DIR}/torsion/torsion.hpp"
		"${CMAKE_CURRENT_SOURCE_DIR}/reinforced_concrete/beam_design.hpp"

 )

target_include_directories(rco_lib_obj
	PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}/material"
		"${CMAKE_CURRENT_SOURCE_DIR}/moment"
		"${CMAKE_CURRENT_SOURCE_DIR}/shear"
		"${CMAKE_CURRENT_SOURCE_DIR}/torsion"
		"${CMAKE_CURRENT_SOURCE_DIR}/reinforced_concrete"

)

#choose Python Version used
###set(PYBIND11_PYTHON_VERSION 3.10 CACHE STRING "Embedding Python using Pybind11")

#copied from vcpkg installed
###find_package(pybind11 CONFIG REQUIRED) #if using vcpkg and placed at source directory

#copied from vcpkg installed
# Note: 8 target(s) were omitted.
###target_link_libraries(rco_lib_obj PUBLIC pybind11::lto pybind11::embed pybind11::module) #if location is in the library, to share must be in PUBLIC (or non PRIVATE)

target_compile_features( rco_lib_obj
	PUBLIC
		cxx_constexpr
)

#create ALIAS target
add_library( rco::lib ALIAS rco_lib_obj )

set_target_properties(rco_lib_obj PROPERTIES DEBUG_POSTFIX "d")
