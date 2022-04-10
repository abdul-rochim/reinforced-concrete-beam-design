#module data.py

##source :
#python read excel rows into tuple
#python read and write excel
#https://stackoverflow.com/questions/61675812/from-excel-to-list-of-tuples
#https://stackoverflow.com/questions/33655127/how-to-read-certain-columns-from-excel-using-pandas-python
#https://qavalidation.com/2021/03/python-openpyxl-library-read-excel-details-or-multiple-rows-as-tuple.html/
#https://automatetheboringstuff.com/2e/chapter13/
#https://zetcode.com/python/openpyxl/

import pandas as pd
#import openpyxl

#data given from structure analysis
Mu_pos = 270.0          #kN.m
Mu_neg = 250.0          #kN.m
fy = 420.0              #MPa
fc = 35.0               #MPa
dia_As_pos = 19.0       #mm
dia_As_neg = 19.0       #mm
dia_hoop = 13.0         #mm
bw = 450.0              #mm
h = 760.0               #mm
cover = 40.0            #mm
dt_pos = h - cover - dia_hoop - dia_As_pos / 2.0          #mm
dt_neg = h - cover - dia_hoop - dia_As_neg / 2.0          #mm

#initial definition of As and number of longitudinal bars
As_req_pos = 0.0        #mm2
As_req_neg = 0.0        #mm2
n_long_bar_pos = 0.0    #pieces
n_long_bar_neg = 0.0    #pieces

Vu = 211.0              #kN
fyv = 420.0             #MPa
dia_stirrup = 13.0      #mm
n_stirrup = 2.0         #leg
phi_shear = 0.75
lambda_ = 1.0

#initial definition of stirrup spacing
s_prov = 500.0              #mm

Tu = 22.50                  #kN.m
fyt = fyv                   #MPa
dia_trans = dia_stirrup     #mm
phi_torsion = 0.75
tetha = 45.0                #degree

#data_analysis = (Mu_pos,Mu_neg,fy,fc,dia_As_pos,dia_As_neg,dia_hoop,bw,h,cover,dt_pos,dt_neg,As_req_pos,As_req_neg,n_long_bar_pos,n_long_bar_neg,Vu,fyv,dia_stirrup,n_stirrup,phi_shear,lambda_,s_prov,Tu,fyt,dia_trans,phi_torsion,tetha)

def data_analysis() :
    return (Mu_pos,Mu_neg,fy,fc,dia_As_pos,dia_As_neg,dia_hoop,bw,h,cover,dt_pos,dt_neg,As_req_pos,As_req_neg,n_long_bar_pos,n_long_bar_neg,Vu,fyv,dia_stirrup,n_stirrup,phi_shear,lambda_,s_prov,Tu,fyt,dia_trans,phi_torsion,tetha)


def data_design_excel() :
    """
    #using pandas
    #https://stackoverflow.com/questions/61675812/from-excel-to-list-of-tuples
    df = pd.read_excel(r'D:\My Software\Pybind11\Thesis\Reinforced_Concrete_Optimization\input_output\input\data_excel.xlsx', sheet_name='data')    #, header=None)
    list_of_tuples = list(df.to_records(index=False))
    i = 0
    for element in list_of_tuples :
        print('element : [', i, '] =' , element)
        i = i + 1
    return list_of_tuples[0]
    """
    """
    #using openpyxl
    #https://stackoverflow.com/questions/61675812/from-excel-to-list-of-tuples
    wb = openpyxl.load_workbook(r'D:\My Software\Pybind11\Thesis\Reinforced_Concrete_Optimization\input_output\input\data_excel.xlsx')
    ws = wb.active
    cells = ws['A2:AB2']
    l = []
    for c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27, c28 in cells:
        l.append((c1.value, c2.value, c3.value, c4.value, c5.value, c6.value, c7.value, c8.value, c9.value, c10.value, c11.value, c12.value, c13.value, c14.value, c15.value, c16.value, c17.value, c18.value, c19.value, c20.value, c21.value, c22.value, c23.value, c24.value, c25.value, c26.value, c27.value, c28.value))
    print(l)
    return l #Also Display The Excel Formula -> Can not use this method.
    """
    #using pandas
    #https://stackoverflow.com/questions/33655127/how-to-read-certain-columns-from-excel-using-pandas-python
    file_location = "~\Reinforced_Concrete_Optimization\input_output\input\data_excel.xlsx"
    df = pd.read_excel(file_location, sheet_name='data', usecols = "A:AE")    #, header=None)
    list_of_tuples = list(df.to_records(index=False))
    """
    i = 0
    for element in list_of_tuples :
        print('element : [', i, '] =' , element)
        i = i + 1
    """
    return list_of_tuples #[0]

