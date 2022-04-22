#include "torsion.hpp"
/*
Calculating torsion capacity due to the responses of the structure analysis.
*/

//Determine the concrete section resisting torsion.
//calculation the threshold torsion value to check
//can torsional effects be neglected ?
std::tuple<bool, double, double, double> can_torsion_be_neglected(
	const double& Tu,
	const double& fc,
	const double& b,
	const double& h,
	const double& phi_torsion,
	const double& lambda
	) {
	//check threshold torsion Tth
	//Tth = lambda * sqrt(fc) * Acp^2 / Pcp
	double Acp = b * h;		//mm2 is the area enclosed by outside perimeter of concrete
	double Pcp = 2 * (b + h);	//mm is the perimeter of concrete gross area
	double Tth = lambda * 0.083 * std::pow(fc, 0.5) * std::pow(Acp, 2.) / Pcp / 1000000.;
	double phi_Tth = phi_torsion * Tth;

	std::ofstream out_file;
	std::string filepath = "../input_output/output/data_output.txt";
	out_file.open(filepath, std::ios::app);
	if (!out_file) {
		std::cerr << "Acces denied! file not created.\n";
		exit(1);
	}
	out_file << std::fixed << std::setprecision(2);

	bool no_torsion_reinf;
	//check if torsion can be ignored; does Tu < phi_Tth?
	if (Tu < phi_Tth) {
		out_file << "Tu : " << Tu << " kN.m < phi_Tth : " << phi_Tth << " kN.m\n";
		out_file << "torsion reinforcement is not required\n";
		no_torsion_reinf = true;
	}
	else {
		out_file << "Tu : " << Tu << " kN.m > phi_Tth : " << phi_Tth << " kN.m\n";
		out_file << "torsion effect cannot be neglected, reinforcement and detailing requirements \nfor torsion must be considered\n";
		no_torsion_reinf = false;
	}
	return std::make_tuple(no_torsion_reinf, Acp, Pcp, Tth);
}

//torsional reinforcement ->reference SP-17M(14)
Torsion torsion_design(
	const double& Tu,
	const double& Vu,
	const double& fyt,
	const double& fy,
	const double& fc,
	const double& b,
	const double& h,
	const double& d,
	const double& phi_torsion,
	const double& lambda,
	const double& dia_transverse,
	const double& cover,
	const double& tetha
	) {
	std::ofstream out_file;
	std::string filepath = "../input_output/output/data_output.txt";
	out_file.open(filepath, std::ios::app);
	if (!out_file) {
		std::cerr << "Acces denied! file not created.\n";
		exit(1);
	}
	else {
		//std::cout << "The file to write the analysis result is available in the directory.\n";
	}

	out_file << "\n-----------------------------------------------------------------------------------\n";
	out_file << "starting to compute torsion reinforcement of concrete beam!\n";
	out_file << "-----------------------------------------------------------------------------------\n";
	out_file << std::fixed << std::setprecision(2);

	//calculate the threshold value:
	auto [can_torsion_be_ignored, Acp, Pcp, Tth] = can_torsion_be_neglected(Tu, fc, b, h, phi_torsion, lambda); //C++17 structured binding -> en.cppreference.com

	if (!can_torsion_be_ignored) {
		//Torsion reinforcement
		//calculate cracking torsion:
		double Tcr = 0.33 * lambda * std::pow(fc, 0.5) * std::pow(Acp, 2.) / Pcp / 1000000.;
		double Tu_used;

		out_file << '\n';
		//check if cross section will crack under the torsional moment.
		if (Tu < Tcr) {
			out_file << "Tu : " << Tu << " kN.m < Tcr : " << Tcr << " kN.m\n";
			out_file << "reducing Tu to Tcr is not required\n";
			Tu_used = Tu;
		}
		else {
			out_file << "Tu : " << Tu << " kN.m > Tcr : " << Tcr << " kN.m\n";
			out_file << "Tu can be reduced to Tcr\n";
			Tu_used = Tcr;
		}

		//check if cross section is adequate to resist the torsional moment
		//Aoh = xo.yo
		//Ph = 2 (xo + yo)
		double Aoh = (b - 2. * cover - dia_transverse) * (h - 2. * cover - dia_transverse);
		double Ph = 2 * ((b - 2. * cover - dia_transverse) + (h - 2. * cover - dia_transverse));

		double Vc = 0.17 * lambda * std::pow(fc, 0.5) * b * d / 1000.; //kN
		double Vu_per_bd = Vu * 1000. / (b * d); // N/mm2
		double Vc_per_bd = Vc * 1000. / (b * d); // N/mm2
		double Tu_Ph_per_Aoh2 = Tu_used * 1000000. * Ph / (1.7 * std::pow(Aoh, 2.)); // N/mm2

		double lhs = std::pow((std::pow(Vu_per_bd, 2.) + std::pow(Tu_Ph_per_Aoh2, 2.)), 0.5);
		double rhs = phi_torsion * (Vc_per_bd + 0.66 * std::pow(fc, 0.5));
		out_file << "\ncheck if cross section is adequate to resist the torsional moment\n";
		try {
			if (lhs <= rhs) {
				out_file << "torsional shear flow : " << lhs << " MPa < resistance : " << rhs << " MPa\n";
				out_file << "Ok, section is adequate to resist torsion\n";
			}
			else {
				std::cout << "torsional shear flow : " << lhs << " MPa > resistance : " << rhs << " MPa\n";
				std::cout << "section is not adequate to resist torsion\n";
				throw std::invalid_argument("WARNING! need to enlarge the cross section!");
			}
		}
		catch (std::invalid_argument& e) {
			std::cerr << e.what() << '\n';
			throw;
			//return -1;
		}

		out_file << "\ncalculate required transverse and longitudinal torsion reinforcement :\n";
		//transverse
		double Ao = 0.85 * Aoh; //is the gross area enclosed by torsional shear flow path
		double At_per_s_req = Tu_used * 1000000. / phi_torsion / (2. * Ao * fyt) * std::tan(double(tetha) * M_PI / 180.); //mm2/mm
		out_file << "At/s required : " << At_per_s_req << " mm2/mm\n";

		//longitudinal
		double Al_req = Tu_used * 1000000. / phi_torsion / (2. * Ao * fy) / std::tan(double(tetha) * M_PI / 180.) * Ph;	//mm2
		out_file << "Al required : " << Al_req << " mm2\n";

		out_file << "\ncalculation of the torsional requirement complete!\n";
		out_file << "-----------------------------------------------------------------------------------\n";
		//out_file.close();
		return{ Tth, Tcr , At_per_s_req, Al_req };
	}
	else {
		out_file << "due to Tu < phi_Tth\n";
		out_file << "hence torsion reinforcement is not required.\n";
		out_file << "\ncalculation of the torsional requirement complete!\n";
		out_file << "-----------------------------------------------------------------------------------\n";
		//out_file.close();
		return{ Tth, 4. * Tth, 0.0, 0.0 };
	}
}
