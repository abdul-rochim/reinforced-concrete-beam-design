#ifndef MODULEHEADER_H
#define MODULEHEADER_H

#include "beam_design.hpp"

#ifdef _MSC_VER
#include <CodeAnalysis/Warnings.h> // ALL_CODE_ANALYSIS_WARNINGS
#pragma warning (push)
#pragma warning (disable : ALL_CODE_ANALYSIS_WARNINGS)
#endif
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#ifdef _MSC_VER
#pragma warning (pop)
#endif

namespace py = pybind11;

#endif