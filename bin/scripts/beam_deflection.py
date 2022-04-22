# module beam_deflection.py
from module_data import beamDataAnalysis
from module_data_excel import beamDataDesignExcel

# function with global variable
def beamDeflection():
	# (Mu_pos, Mu_neg, fy, fc, dia_As_pos, dia_As_neg, dia_hoop, bw, h, cover, dt_pos, dt_neg,
	# As_req_pos, As_req_neg, n_long_bar_pos, n_long_bar_neg, Vu, fyv, dia_stirrup, n_stirrup,
	# phi_shear, lambda_, s_prov, Tu, fyt, dia_trans, phi_torsion, tetha) = beamDataAnalysis()

	(beam_name, Mu_pos, Mu_neg, fy, fc, dia_As_pos, dia_As_neg, dia_hoop, bw, h, cover, dt_pos, dt_neg, As_req_pos, As_req_neg, n_long_bar_pos, n_long_bar_neg, Vu, fyv, dia_stirrup, n_stirrup, phi_shear, lambda_, s_prov, Pu, Tu, fyt, dia_trans, phi_torsion, tetha, dia_As_sideFace) = beamDataDesignExcel()[0]
	# (beam_name, *x, dia_As_sideFace) = beamDataDesignExcel()[0]
	# print("beam_name", beam_name)
	# print("diameter of side face bar", dia_As_sideFace)
