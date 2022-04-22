# module main.py

import cpp_beam_design
from module_data import beamDataAnalysis
from module_data_excel import beamDataDesignExcel

# function with global variable
def beam_design_summary():
    # print("[Python] call function beam_design from C++\n")

    # calc = cpp_beam_design.beamDesign(data_analysis)
    # calc = cpp_beam_design.beamDesign(data_analysis())
    # calc = cpp_beam_design.beamDesign(beamDataAnalysis())
    # calc = cpp_beam_design.beamDesign(beamDataDesignExcel()[0])
    # calc = cpp_beam_design.beamDesign(beamDataDesignExcel()[1])
    
    calc = []
    for i in range(len(beamDataDesignExcel())):
        calc.append(cpp_beam_design.beamDesign(beamDataDesignExcel()[i]))
    """
    print("print element from Python")
    for x in range(len(calc)):
        print(calc[x])
    print("end print element from Python\n")
    """
    return calc
