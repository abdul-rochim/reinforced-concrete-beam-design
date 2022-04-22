#include "shear.hpp"
/*
Calculating shear capacity due to the responses of the structure analysis.
data of material such as concrete and steel reinforced strength and the dimensions are given.
*/

Shear shear_design(
	const double& Vu,
	const double& fyv,
	const double& fc,
	const double& b,
	const double& d,
	const double& phi_shear,
	const double& lambda_,
	const double& dia_stirrup,
	const double& n_stirrup,
	double& s_provided,
	const double& Pu
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
	out_file << "starting the iteration to compute shear design of concrete beam!\n";
	out_file << "-----------------------------------------------------------------------------------\n";
	out_file << std::fixed << std::setprecision(2);

	//shear strength
	using std::pow;
	double Vc = 0.17 * lambda_ * pow(fc, 0.5) * b * d / 1000.; //kN
	double phi_Vc = phi_shear * Vc;

	//Check if phi_Vc >= Vu
	if (phi_Vc < Vu) {
		out_file << "phi Vc : " << phi_Vc << " kN < Vu : " << Vu << " kN\n";
		out_file << "therefore, shear reinforcement is required!\n";
	}
	else {
		out_file << "phi Vc : " << phi_Vc << " kN > Vu : " << Vu << " kN\n";
		out_file << "therefore, shear reinforcement is not required!\n";
		out_file << "but, recommend to use the minimum shear reinforcement\n";
	}

	//prior to calculating shear reinforcement, check if the cross-sectional dimensions satisfy
	double phi_Vn_maximum = phi_shear * (Vc + 0.66 * pow(fc, 0.5) * b * d / 1000.); //kN
	try {
		if (Vu < phi_Vn_maximum) {
			out_file << "phi_Vn_maximum : " << phi_Vn_maximum << " kN > Vu : " << Vu << " kN\n";
			out_file << "Ok therefore, cross-sectional dimensions satisfy!\n";
		}
		else {
			std::cout << "phi_Vn_maximum : " << phi_Vn_maximum << " kN < Vu : " << Vu << " kN\n";
			std::cout << "recommend to enlarge cross-sectional dimensions or to revise material properties\n";
			throw std::invalid_argument("WARNING! Vu > phi_Vn_maximum!");
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what() << '\n';
		throw;
		//return -1;
	}

	// check number of stirrup
	try {
		if (n_stirrup >= 2.0) {
			//do nothing
		}
		else {
			std::cout << "\nWrong input of number of stirrup, the input is only 1 leg.\n";
			throw std::invalid_argument("WARNING! minimum number of stirrup is 2 legs.");
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what() << '\n';
		throw;
		//return -1;
	}

	//shear reinforcement
	//tranverse reinforcement satisfying Eq.(22.5.10.1) is required at each section where Vu > phi_Vc
	double Vs_req = Vu / phi_shear - Vc;

	//where Vs = Av fyt d / s
	//calculate Av per s required = Vs / (fyt * d)
	double Av_per_s_req = Vs_req * 1000.0 / (fyv * d);
	double s_req;

	if (Av_per_s_req < 0.0) {
		Av_per_s_req = 0.0;
		s_req = 600.0;	//assumption at the maximum stirrup spacing requirement
	}
	else {
		//calculate maximum allowable stirrup spacing:
		s_req = (n_stirrup * 0.25 * M_PI * dia_stirrup * dia_stirrup) / Av_per_s_req;
	}
	out_file << "initial s required(maximum allowable stirrup spacing): " << s_req << " mm\n";

	if (s_provided < 40. || s_provided > 600.)
		s_provided = 600.; //assumption initial s

	out_file << "initial s provided : " << s_provided << " mm\n\n";

	//first does the beam transverse reinforcement value need to exceed
	//the threshold value of Vs <= 0.33 * sqrt(fc') bd?
	//calculate Vs_maximum -> will be evaluated inside of the iteration
	double Vs_maximum = 0.33 * pow(fc, 0.5) * b * d / 1000.; //kN

	//In the region where Vu <= phi_Vc /2, shear reinforcement is not required.
	//Check is Vu <= phi_Vc /2 ?
	if (Vu <= phi_Vc / 2.) {
		out_file << "Vu : " << Vu << " kN < phi_Vc/2 : " << phi_Vc / 2. << " kN\n";
		out_file << "therefore, shear reinforcement is not required!\n";
		out_file << "but, recommend to use minimum shear reinforcement\n\n";
	}
	else {
		out_file << "Vu : " << Vu << " kN > phi_Vc/2 : " << phi_Vc / 2. << " kN\n";
		out_file << "check, is Av_provided greater than minimum shear reinforcement?\n\n";
	}

	//specified shear reinforcement must be at least:
	//Av_min per s = 0.062 * sqrt(fc') bw/ fyt	//requirement 1
	//Av_min per s = 0.35 * bw/ fyt				//requirement 2
	double Av_per_s_minimum_1 = 0.062 * pow(fc, 0.5) * b / fyv;
	double Av_per_s_minimum_2 = 0.35 * b / fyv;
	double Av_per_s_minimum = std::max(Av_per_s_minimum_1, Av_per_s_minimum_2);

	//area of shear bars
	double Av_prov = n_stirrup * 0.25 * M_PI * dia_stirrup * dia_stirrup;
	
	//declare variables of Av per s prov, Vs prov, phi_Vn prov
	double Av_per_s_provided;
	double Vs_provided;
	double phi_Vn_provided;

	//calculate s maximum
	//the maximum stirrup spacing is the lesser of d/2 and 600 mm, and Av_prov/(Av/s minimum)
	double s_max = std::min(Av_prov / Av_per_s_minimum, std::min(d / 2., 600.));

	//iteration for computing s provided dan check to the minimum shear bars
	unsigned int j = 1;
	do {
		out_file << "iteration_" << j << '\n';

		Av_per_s_provided = Av_prov / s_provided;
		out_file << "Av_per_s_provided : " << Av_per_s_provided << " mm2/m\n";

		//Vs = Av * fyt * d/s
		Vs_provided = Av_prov * fyv * d / s_provided /1000.; //kN
		out_file << "Vs_provided : " << Vs_provided << " kN\n";

		//first does the beam transverse reinforcement value need to exceed
		//the threshold value of Vs <= 0.33 * sqrt(fc') bd?
		//Vs provided should be lesser than Vs maximum
		try {
			if (Vs_provided <= Vs_maximum)
				out_file << "Ok, Vs_provided : " << Vs_provided << " kN < Vs maximum : " << Vs_maximum << " kN\n";
			else {
				std::cout << "Vs_provided : " << Vs_provided << " kN > Vs maximum : " << Vs_maximum << " kN\n";
				std::cout << "enlarge cross-sectional dimensions\n";
				throw std::invalid_argument("WARNING! Vs_provided > Vs_maximum");
			}
		}
		catch (std::invalid_argument& e) {
			std::cerr << e.what() << '\n';
			throw;
			//return -1;
		}

		//check Av per s minimum
		if (Av_per_s_provided >= Av_per_s_minimum) {
			out_file << "Av_per_s prov : " << Av_per_s_provided << " mm2/mm > Av_per_s_minimum: " << Av_per_s_minimum << " mm2/mm\n";
			out_file << "Ok, stirrup spacing satisfies the minimum shear reinforcement\n";
		}
		else {
			out_file << "Av_per_s prov : " << Av_per_s_provided << " mm2/mm < Av_per_s_minimum: " << Av_per_s_minimum << " mm2/mm\n";
			out_file << "stirrup spacing is lesser than the minimum shear reinforcement\n";
		}
		out_file << "s required : " << s_req << " mm\n";
		if (s_provided > s_req)
			out_file << "s provided iteration : " << s_provided << " mm > s required\n";
		else
			out_file << "s provided iteration : " << s_provided << " mm < s required\n";
		
		//because the required shear strength is below the threshold value,
		//the maximum stirrup spacing is the lesser of d/2 and 600 mm, and Av_prov/(Av/s minimum)
		out_file << "\ns maximum requirement :\n";
		out_file << "s provided must be lesser than d/2 : " << d / 2. << " mm and " << 600.0 << " mm\n";
		out_file << "s provided must be lesser than Av prov/ (Av/s min) : " << Av_prov / Av_per_s_minimum << " mm\n";
		//display s maximum required
		out_file << "s maximum required : " << s_max << " mm\n";

		if(s_provided > s_max )
			out_file << "s provided iteration : " << s_provided << " mm > s maximum\n";
		else
			out_file << "s provided iteration : " << s_provided << " mm < s maximum\n";

		//compute nominal shear capacity
		phi_Vn_provided = phi_shear * (Vc + Vs_provided);
		//check is Vu < phi Vn ?
		if (phi_Vn_provided >= Vu) {
			out_file << "\nphi Vn provided iteration : " << phi_Vn_provided << " kN > Vu : " << Vu << " kN\n";
			if (s_provided > s_max)
				out_file << "however s provided > s max; therefore re-calculate again!\n";
		}
		else {
			out_file << "\ns provided > s required\n";
			out_file << "phi Vn provided iteration : " << phi_Vn_provided << " kN < Vu : " << Vu << " kN, therefore re-calculate again!\n";
		}
		out_file << '\n';

		s_provided = s_provided - 10.;
		++j;
	} while (s_provided + 10. > s_max || phi_Vn_provided < Vu);	//s should be lesser than s_max, and phi Vn provided should be greater than Vu.

	out_file << "calculation of the nominal shear capacity complete!\n";
	out_file << "-----------------------------------------------------------------------------------\n";
	s_provided = s_provided + 10.;
	//out_file.close();
	return { Vc, Vs_provided, Av_per_s_req, Av_per_s_provided };	//return rvalue (using constructor) //return Shear{ Vc, Vs_provided, Av_per_s_req, Av_per_s_provided };
}
