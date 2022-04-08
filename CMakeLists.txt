##source :
##https://www.youtube.com/watch?v=hmMcl2TVViw
##https://github.com/Bagira80/More-Modern-CMake

cmake_minimum_required(VERSION 3.10...3.22)

project("Reinforced_Concrete_Optimization" VERSION 1.0 DESCRIPTION "the stuff of Thesis")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#set C++ language standard
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)


#include(GNUInstallDirs)
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
#    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
#    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
#    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )

if(MSVC)
    #enable code analysis
    set_property(GLOBAL PROPERTY VS_GLOBAL_EnableCppCoreCheck true)
    set_property(GLOBAL PROPERTY VS_GLOBAL_CodeAnalysisRuleSet CppCoreCheckRules.ruleset)
    set_property(GLOBAL PROPERTY VS_GLOBAL_RunCodeAnalysis true)

    #choose static standard library
    set(CompilerFlags
        CMAKE_CXX_FLAGS
        CMAKE_CXX_FLAGS_DEBUG
        CMAKE_CXX_FLAGS_RELEASE
        CMAKE_C_FLAGS
        CMAKE_C_FLAGS_DEBUG
        CMAKE_C_FLAGS_RELEASE
    )
    foreach(CompilerFlag ${CompilerFlags})
        string(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
    endforeach()

    #remove /W3 from defaults
    string(REGEX REPLACE "/W3" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
    string(REGEX REPLACE "-W3" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
endif()

add_subdirectory(pybind11)

add_subdirectory(external)

add_subdirectory(rco_library)

add_subdirectory(app)

#enable_testing()
#add_subdirectory(tests)
