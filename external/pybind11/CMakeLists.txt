# ./external/pybind11/CMakeLists.txt

add_library(external_pybind11 INTERFACE IMPORTED GLOBAL)

#choose Python Version used
set(PYBIND11_PYTHON_VERSION 3.10 CACHE STRING "Embedding Python using Pybind11")

#copied from vcpkg installed
###find_package(pybind11 CONFIG REQUIRED) #if using vcpkg and must be placed at source directory

#copied from vcpkg installed
# Note: 8 target(s) were omitted.
target_link_libraries(external_pybind11 INTERFACE pybind11::lto pybind11::embed pybind11::module) #must be INTERFACE due to add_library set"INTERFACE IMPORTED GLOBAL"

#create ALIAS target
add_library( external::pybind11 ALIAS external_pybind11 )