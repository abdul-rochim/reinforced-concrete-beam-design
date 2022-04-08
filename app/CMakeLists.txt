
cmake_minimum_required(VERSION 3.10...3.22)

project("RCO_App" VERSION 1.0 DESCRIPTION "the main program")

set(TARGET_NAME RCO)

#choose Python Version used
###set(PYBIND11_PYTHON_VERSION 3.10 CACHE STRING "Embedding Python using Pybind11")

#copied from vcpkg installed
###find_package(pybind11 CONFIG REQUIRED) #if using vcpkg

add_executable(${TARGET_NAME}
                "src/Main.cpp"
                "src/pymodule/moduleheader.h"
		        "src/pymodule/modulebeamdesign.h"
                "src/pymodule/example.h"
)

######
#pybind11_add_module()
######

#copied from vcpkg installed
# Note: 8 target(s) were omitted.
###target_link_libraries(${TARGET_NAME} PRIVATE pybind11::lto pybind11::embed pybind11::module)  #in this execution directory, External Library could set as PRIVATE

target_link_libraries(${TARGET_NAME}
    PUBLIC
        rco::lib
        external::pybind11
)

#set binary filenames
set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME_DEBUG ${TARGET_NAME}-debug )
set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME_RELWITHDEBINFO ${TARGET_NAME}-relwithdebinfo )
set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME_RELEASE ${TARGET_NAME}-release )
set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME_MINSIZEREL ${TARGET_NAME}-minsizerel )

#set warning levels
if(MSVC)
    target_compile_options(${TARGET_NAME} PUBLIC /W4 /permissive)
else()
    target_compile_options(${TARGET_NAME} PUBLIC -Wall -Wextra -pedantic)
endif()
