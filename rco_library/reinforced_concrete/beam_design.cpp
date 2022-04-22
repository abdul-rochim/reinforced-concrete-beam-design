#include "beam_design.hpp"

/*
summary of the reinforced concrete beam design.
*/

std::array<double, 5> beam_design(
	const std::tuple<
	std::string,
	double,
	double,
	double,
	double,
	double,
	double,
	double,
	double,
	double,
	double,
	double,
	double,
	double*,
	double*,
	double*,
	double*,
	double,
	double,
	double,
	double,
	double,
	double,
	double*,
	double,
	double,
	double,
	double,
	double,
	double,
	double >& tup
	)
{
	//C++17 structured binding: <en.cppreference.com>
	auto [beam_name, Mu_pos, Mu_neg, fy, fc, dia_As_pos, dia_As_neg, dia_hoop, bw, h, cover, dt_pos, dt_neg, As_req_pos, As_req_neg, n_long_bar_pos, n_long_bar_neg,
		Vu, fyv, dia_stirrup, n_stirrup, phi_shear, lambda, s_prov, Pu,
		Tu, fyt, dia_trans, phi_torsion, tetha, dia_As_sideFace] = tup;

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
	std::cout << "the analysis of beam: " << beam_name << " is being computed...\n\n";

	out_file << "-----------------------------------------------------------------------------------\n";
	out_file << "                                      Beam : " << beam_name << '\n';
	out_file << "-----------------------------------------------------------------------------------\n";

	out_file << "Calculation of nominal moment capacity - positive\n";
	out_file.close();
	Moment m_pos(Mu_pos, fy, fc, bw, dt_pos);
	//binding As_req and n_long_bar
	//double& As_req_pos_used = *As_req_pos;
	//double& n_long_bar_pos_used = *n_long_bar_pos;
	//std::cout << "As_req_pos_initial : " << As_req_pos_used << " mm2 and n_bar_pos_initial : " << n_long_bar_pos_used << " pieces\n";
	auto momentDesign_positive = phiMn(m_pos, dia_As_pos, dia_hoop, cover, As_req_pos, n_long_bar_pos);
	//return double phiMn, *As_provided and number of bars provided - > *nb_bar
	//std::cout << "As_req_pos_used : " << As_req_pos_used << " mm2 and n_bar_pos_used : " << n_long_bar_pos_used << " pieces\n";

	out_file.open(filepath, std::ios::app);
	out_file << "\nCalculation of nominal moment capacity - negative\n";
	out_file.close();
	Moment m_neg(Mu_neg, fy, fc, bw, dt_neg);
	//binding As_req and n_long_bar
	//double& As_req_neg_used = *As_req_neg;
	//double& n_long_bar_neg_used = *n_long_bar_neg;
	//std::cout << "As_req_neg_initial : " << As_req_neg_used << " mm2 and n_bar_neg_initial : " << n_long_bar_neg_used << " pieces\n";
	auto momentDesign_negative = phiMn(m_neg, dia_As_neg, dia_hoop, cover, As_req_neg, n_long_bar_neg);
	//return double phiMn, *As_provided and number of bars provided -> *nb_bar
	//std::cout << "As_req_neg_used : " << As_req_neg_used << " mm2 and n_bar_neg_used : " << n_long_bar_neg_used << " pieces\n";
	
	//binding reference s_prov
	//double& s_prov_used = *s_prov;
	//std::cout << "s provided_initial : " << s_prov_used << " mm\n";
	auto shearDesign = shear_design(Vu, fyv, fc, bw, dt_pos, phi_shear, lambda, dia_stirrup, n_stirrup, *s_prov, Pu);
	//return object of class Shear -> Vc, Vs_prov, Av_per_s_required, Av_per_s_prov
	//std::cout << "s provided_used : " << s_prov_used << " mm\n";

	auto torsionalDesign = torsion_design(Tu, Vu, fyt, fy, fc, bw, h, dt_pos, phi_torsion, lambda, dia_trans, cover, tetha);
	//return object of class Torsion -> Tth, Tcr, At_per_s_req, Al_req

	out_file.open(filepath, std::ios::app);
	out_file << std::fixed << std::setprecision(2);
	//out_file << "\n-----------------------------------------------------------------------------------\n";
	//out_file << "starting to compute summary of beam design\n\n";

	if (torsionalDesign.At_per_s_req == 0.0 || torsionalDesign.Al_req == 0.0) {
		out_file << "\n-----------------------------------------------------------------------------------\n";
		out_file << "                  the summary of reinforced concrete beam design\n";
		out_file << "-----------------------------------------------------------------------------------\n";
		out_file << "reinforcement detailing:\n";

		out_file << std::fixed << std::setprecision(0);
		out_file << "detailing of longitudinal bars :\n";
		out_file << "\tpositive reinforcement : " << *n_long_bar_pos << " D" << dia_As_pos << '\n';
		out_file << "\tnegative reinforcement : " << *n_long_bar_neg << " D" << dia_As_neg << '\n';

		out_file << '\n';
		out_file << "detailing of transverse bars :\n";
		out_file << "\tnumber of leg-space  : " << n_stirrup << " D" << dia_trans << '-' << *s_prov << '\n';
	}
	else { //Calculation and Recapitulation
		out_file << "Compute summary of transverse reinforcement.\n";

		//the required area for shear and torsion transverse reinforcement are additive:
		out_file << "Av/s required : " << shearDesign.Av_per_s_required << " mm2/mm\n";
		out_file << "At/s required : " << torsionalDesign.At_per_s_req << " mm2/mm\n";
		double Avt_per_s_req = shearDesign.Av_per_s_required + 2.0 * torsionalDesign.At_per_s_req; //Avt per s required = Av/s + 2 At/s
		out_file << "Avt/s required : " << Avt_per_s_req << " mm2/mm\n";

		// Avt/s minimum
		double Avt_per_s_min = std::max(
			0.062 * std::pow(fc, 0.5) * bw / fyt,
			0.35 * bw / fyt
		);
		out_file << "Avt/s minimum : " << Avt_per_s_min << " mm2/mm\n\n";
		
		// initial declaration Avt/s provided
		double Avt_per_s_prov_alt_1{};
		double Avt_per_s_prov_alt_2{};

		out_file << "maximum spacing of transverse shear reinforcement must not exceed,\nthe lesser of d/2 and 600 mm.\n";
		out_file << "maximum spacing of transverse torsional reinforcement must not exceed,\nthe lesser of ph/8 and 300 mm.\n";
		double max_shear_spacing = std::min(600.0, std::min(dt_pos, dt_neg) / 2.0);
		out_file << "maximum spacing of transverse shear reinforcement [requirement] : " << max_shear_spacing << " mm\n";
		
		double ph = 2.0 * ((bw - 2.0 * cover - dia_trans) + (h - 2.0 * cover - dia_trans));
		double max_torsion_spacing = std::min(ph/8.0, 300.0);
		out_file << "maximum spacing of transverse torsional reinforcement [requirement] : " << max_torsion_spacing << " mm\n";

		double max_s = std::min(max_torsion_spacing, max_shear_spacing);
		out_file << "maximum spacing of transverse reinforcement [requirement] used : " << max_s << " mm\n\n";

		//maximum space from the threshold At/s required.
		////=[ max spacing = 1.0 * 0.25 * PI * dia_trans * dia_trans / torsionalDesign.At_per_s_req ]=////
		double max_spacing_for_torsion{ 600.0 }; //initial declaration // assumtion -> maximum spacing = 600.0 mm
		
		double s_prov_alt_1{ }; // { *s_prov }; 
		double s_prov_alt_2{ }; // { *s_prov };
		double n_stirrup_alt_1{ n_stirrup };
		double n_stirrup_alt_2{ n_stirrup };

		//s torsion requirement
		out_file << "calculate maximum spacing torsion required :\n";
		double At_per_s_prov_dummy{}; // dummy At/s to calculate maximum spacing due to torsional
		unsigned int i{ 1 };
		do {
			out_file << "iter i: " << i << '\n';

			if (max_spacing_for_torsion > max_s)
				out_file << "s torsion max iter : " << max_spacing_for_torsion << " mm > max spacing requirement : " << max_s << " mm\n";

			At_per_s_prov_dummy = 1.0 * 0.25 * PI * dia_trans * dia_trans / max_spacing_for_torsion; // 1 leg torsion
			out_file << "At_per_s_prov_dummy : " << At_per_s_prov_dummy << " mm2/mm\n";

			max_spacing_for_torsion = max_spacing_for_torsion - 10.0;
			++i;
		} while (At_per_s_prov_dummy < torsionalDesign.At_per_s_req || max_spacing_for_torsion + 10.0 > max_s);
		max_spacing_for_torsion = max_spacing_for_torsion + 10.0;
		
		out_file << "s torsion required : " << 1.0 * 0.25 * PI * dia_trans * dia_trans/ torsionalDesign.At_per_s_req << " mm\n";
		out_file << "s torsion maximum : " << max_spacing_for_torsion << " mm\n";
		
		out_file << '\n';
		out_file << "Alternative 1 [for 2 legs stirrup] :\n";
		out_file << "number of leg before iter alternative 1 : " << n_stirrup_alt_1 << " legs\n";

		s_prov_alt_1 = max_spacing_for_torsion; //initial assignment for spacing provided
		out_file << "spacing of stirrup before iter alternative 1 (assume = max_spacing_for_torsion) : " << s_prov_alt_1 << " mm\n";

		//alternative 1
		double n_stirrup_alt_1_dummy{ n_stirrup };
		
		if (n_stirrup_alt_1_dummy > 2.0){
			// do nothing
			out_file << "there's no alternative 1 for 2 legs stirrup.\n";
		}
		else {
			if (n_stirrup_alt_1 == 2.0) {
				unsigned int j{ 1 };
				do {
					out_file << "iter j: " << j << '\n';
					Avt_per_s_prov_alt_1 = n_stirrup_alt_1 * 0.25 * PI * dia_trans * dia_trans / s_prov_alt_1;
					out_file << "Avt/s prov : " << Avt_per_s_prov_alt_1 << " mm2/mm\n";

					s_prov_alt_1 = s_prov_alt_1 - 10.0;
					++j;
				} while (Avt_per_s_prov_alt_1 <= std::max(Avt_per_s_req, Avt_per_s_min));
			}
			s_prov_alt_1 = s_prov_alt_1 + 10.0;
		}
		
		out_file << "\nAvt/s prov alternative 1 : " << Avt_per_s_prov_alt_1 << " mm2/mm\n";
		out_file << "number of leg after iter alternative 1 : " << n_stirrup_alt_1 << " legs\n";
		out_file << "spacing of stirrup after iter alternative 1 : " << s_prov_alt_1 << " mm\n\n";

		out_file << "Alternative 2 [for more than 2 legs stirrup] :\n";
		out_file << "number of leg before iter alternative 2 : " << n_stirrup_alt_2 << " legs\n";
		out_file << "spacing of stirrup before iter alternative 2 : " << s_prov_alt_2 << " mm\n";

		//alternative 2
		if (n_stirrup_alt_2 > 2.0) {
			n_stirrup_alt_2 = n_stirrup;
		}
		else {
			n_stirrup_alt_2 = n_stirrup_alt_2 + 1.0;
		}

		if (n_stirrup_alt_2 > 2.0) { // alternative is provided stirrup reinf. more than 2 legs
			unsigned int k{ 1 };
			do {
				//maximum spacing of transverse torsional reinforcement must not exceed the lesser of ph/8 and 300 mm
				s_prov_alt_2 = 300.0; //assumption 300 mm
				out_file << "iter k: " << k << '\n';

				unsigned int j{ 1 };
				do {
					out_file << "\titer j: " << j << '\n';
					Avt_per_s_prov_alt_2 = n_stirrup_alt_2 * 0.25 * PI * dia_trans * dia_trans / s_prov_alt_2;
					out_file << "\tAvt/s prov : " << Avt_per_s_prov_alt_2 << " mm2/mm\n";

					s_prov_alt_2 = s_prov_alt_2 - 10.0;
					++j;
				} while (Avt_per_s_prov_alt_2 <= std::max(Avt_per_s_req, Avt_per_s_min));
				s_prov_alt_2 = s_prov_alt_2 + 10.0;

				out_file << "n_stirrup_alt_2 : " << n_stirrup_alt_2 << " legs\n";
				++n_stirrup_alt_2;
				++k;
			} while (s_prov_alt_2 < 90.0); // 100mm
			--n_stirrup_alt_2;
		}

		out_file << "\nAvt/s prov alternative 2 : " << Avt_per_s_prov_alt_2 << " mm2/mm\n";
		out_file << "number of leg after iter alternative 2 : " << n_stirrup_alt_2 << " legs\n";
		out_file << "spacing of stirrup after iter alternative 2 : " << s_prov_alt_2 << " mm\n";
		if (s_prov_alt_2 > max_spacing_for_torsion) {
			out_file << "spacing stirrup : " << s_prov_alt_2 << " mm > maximum spacing required : " << max_spacing_for_torsion << " mm\n";
			s_prov_alt_2 = max_spacing_for_torsion;
			out_file << "spacing of stirrup after iter alternative 2 used : " << s_prov_alt_2 << " mm\n";
		}

		out_file << "\nCompute summary of logitudinal reinforcement.\n";
		
		//double Avt_per_s_req = shearDesign.Av_per_s_required + 2 * torsionalDesign.At_per_s_req; //Avt per s required = Av/s + 2 At/s
		//calculate At/s (aproximation)
		double Av_per_s_aprox_prov = std::min(Avt_per_s_prov_alt_1, Avt_per_s_prov_alt_2) * shearDesign.Av_per_s_required / Avt_per_s_req; // 1.0 * Av/s
		double At_per_s_aprox_prov = std::min(Avt_per_s_prov_alt_1, Avt_per_s_prov_alt_2) * torsionalDesign.At_per_s_req / Avt_per_s_req; // 2.0 * At/s

		out_file << "calculation of At/s prov [aproximation]\n";
		out_file << "minimum(Avt_per_s_prov_alt_1, Avt_per_s_prov_alt_2) : " << std::min(Avt_per_s_prov_alt_1, Avt_per_s_prov_alt_2) << " mm2/mm\n";
		out_file << "Av_per_s_aprox_prov = ratio Av/s required * Avt/s prov = " << Av_per_s_aprox_prov << " mm2/mm\n";
		out_file << "At_per_s_aprox_prov = ratio At/s required * Avt/s prov = " << At_per_s_aprox_prov << " mm2/mm\n";

		out_file << "\nThe torsional longitudinal reinforcement, Al_min must be the lesser of:\n";
		double Acp = bw * h; //mm2
		double Al_min_req_1 = 0.42 * std::pow(fc, 0.5) * Acp / fy - At_per_s_aprox_prov * ph * fyt / fy;
		double Al_min_req_2 = 0.42 * std::pow(fc, 0.5) * Acp / fy - (0.175 * bw /fyt) * ph * fyt / fy;

		out_file << "Al_min_req_1 = 0.42 * sqrt(fc, 0.5) * Acp / fy - At_per_s_aprox_prov * ph * fyt / fy = " << Al_min_req_1 << " mm2\n";
		out_file << "Al_min_req_2 = 0.42 * sqrt(fc, 0.5) * Acp / fy - (0.175 * bw /fyt) * ph * fyt / fy = " << Al_min_req_2 << " mm2\n";
		double Al_min_req = std::min(Al_min_req_1, Al_min_req_2);
		out_file << "Al_min required : " << Al_min_req << " mm2\n";

		//check Al required
		out_file << "\ncheck Al required calculation :\n";
		out_file << "Al required : " << torsionalDesign.Al_req << " mm2\n";
		double Al_req_calc_used{};
		if (torsionalDesign.Al_req > Al_min_req) {
			out_file << "Al required calc : " << torsionalDesign.Al_req << " mm2 > Al_min : " << Al_min_req << " mm2\n";
			Al_req_calc_used = torsionalDesign.Al_req;
			out_file << "Al required calc used : " << Al_req_calc_used << " mm2\n";
		}
		else {
			out_file << "Al required calc : " << torsionalDesign.Al_req << " mm2 < Al_min : " << Al_min_req << " mm2\n";
			Al_req_calc_used = Al_min_req;
			out_file << "Al required calc used : " << Al_req_calc_used << " mm2\n";
		}
		out_file << "the longitudinal reinforcement must be added to the flexural reinforcement.\n";

		out_file << "\nTorsion longitudinal reinforcement, Al, must be distributed around the cross section\n";
		out_file << "and the portion of Al that needs to be placed where As is needed is added to As.\n";
		double sideFace = Al_req_calc_used / 2.0; //mm2
		double n_bar_sideFace = sideFace / (0.25 * PI * dia_As_sideFace * dia_As_sideFace);
		out_file << "number of bar side face required : " << n_bar_sideFace << " bars\n";
		out_file << "number of bar side face required (round) : " << std::round(n_bar_sideFace) << " bars\n";
		if (static_cast<int>(std::round(n_bar_sideFace)) % 2 == 1) {
			n_bar_sideFace = std::round(n_bar_sideFace) + 1.0;
		}
		else {
			n_bar_sideFace = std::round(n_bar_sideFace);
		}
		out_file << "number of bar side face used : " << n_bar_sideFace << " bars\n";

		out_file << "calculate delta Al to be placed where As is needed:\n";
		double delta_Al = Al_req_calc_used - (n_bar_sideFace * 0.25 * PI * dia_As_sideFace * dia_As_sideFace);
		out_file << "delta Al : " << delta_Al << " mm2\n";
		if (delta_Al < 0.0)
			delta_Al = 0.0;
		out_file << "delta Al used : " << delta_Al << " mm2\n";

		out_file << "\nthe result of As flexural requirement calculation:\n";
		out_file << "As flexural required positive : " << *As_req_pos << " mm2\n";
		out_file << "As flexural required negative : " << *As_req_neg << " mm2\n";

		//As flexural added delta Al
		out_file << "calculate As flexural added by delta Al:\n";
		out_file << "As_flex_plus_Al_torsion = As flexural + delta Al/2\n";
		double As_flex_plus_Al_torsion_pos = *As_req_pos + delta_Al / 2.0;
		out_file << "As_flex_plus_Al_torsion positive : " << As_flex_plus_Al_torsion_pos << " mm2\n";

		double As_flex_plus_Al_torsion_neg = *As_req_neg + delta_Al / 2.0;
		out_file << "As_flex_plus_Al_torsion positive : " << As_flex_plus_Al_torsion_neg << " mm2\n";
		
		out_file << "\ncalculate number of longitudial bars:\n";
		double n_flex_plus_torsion_pos = As_flex_plus_Al_torsion_pos / (0.25 * PI * dia_As_pos * dia_As_pos);
		out_file << "number of flexural plus torsion positive : " << n_flex_plus_torsion_pos << " bars\n";
		n_flex_plus_torsion_pos = std::ceil (n_flex_plus_torsion_pos);
		out_file << "number of flexural plus torsion positive (roundup) : " << n_flex_plus_torsion_pos << " bars\n";

		double n_flex_plus_torsion_neg = As_flex_plus_Al_torsion_neg / (0.25 * PI * dia_As_neg * dia_As_neg);
		out_file << "\nnumber of flexural plus torsion negative : " << n_flex_plus_torsion_neg << " bars\n";
		n_flex_plus_torsion_neg = std::ceil(n_flex_plus_torsion_neg);
		out_file << "number of flexural plus torsion negative (roundup) : " << n_flex_plus_torsion_neg << " bars\n";

		out_file << '\n';
		//check n_bar provided for longitudinal reinforcement and check s_provided for transverse reinforcement
		//longitudinal reinforcement -> because As required(moment) is checked for As_min, therefore n_long_bar(moment) always > n_flex_plus_torsion
		
		//transverse reinforcement
		if (s_prov_alt_1 > *s_prov) {
			out_file << "s provided for transverse[shear + torsion] alternative 1: " << s_prov_alt_1 << " mm < s provided [judgement] : " << *s_prov << " mm\n";
			out_file << "therefore, s provided alternative 1 : " << *s_prov << " mm\n";
			s_prov_alt_1 = *s_prov;
		}
		if (s_prov_alt_2 > *s_prov) {
			out_file << "s provided for transverse[shear + torsion] alternative 2 : " << s_prov_alt_2 << " mm < s provided [judgement] : " << *s_prov << " mm\n";
			out_file << "therefore, s provided alternative 2 : " << *s_prov << " mm\n";
			s_prov_alt_2 = *s_prov;
		}

		out_file << "\n-----------------------------------------------------------------------------------\n";
		out_file << "                  the summary of reinforced concrete beam design\n";
		out_file << "-----------------------------------------------------------------------------------\n";
		out_file << "reinforcement detailing:\n";

		out_file << std::fixed << std::setprecision(0);
		out_file << "detailing of longitudinal bars :\n";
		out_file << "\tAs required positive reinforcement: " << As_flex_plus_Al_torsion_pos << " mm2\n";
		out_file << "\tpositive moment reinforcement     : " << *n_long_bar_pos << " D" << dia_As_pos << '\n';
		if (delta_Al == 0.0) {
			out_file << "\ttorsion reinforcement             : " << "< nothing >" << "\n\t[ no torsion bar added to As flexural, revise(reduce) bar diameter of side face ]\n";
		}
		else if (n_flex_plus_torsion_pos == *n_long_bar_pos) {
			out_file << "\ttorsion reinforcement             : " << "< small value >\n\t[ torsion bars added in longitudinal bars ]" << '\n';
		}
		else {
			out_file << "\ttorsion reinforcement             : " << n_flex_plus_torsion_pos - *n_long_bar_pos << " D" << dia_As_pos << '\n';
		}
		out_file << "\t-----------------------------------------\n";
		out_file << "\ttherefore, positive reinforcement : " << n_flex_plus_torsion_pos << " D" << dia_As_pos << '\n';

		out_file << '\n';
		out_file << "\tAs required negative reinforcement: " << As_flex_plus_Al_torsion_neg << " mm2\n";
		out_file << "\tnegative moment reinforcement     : " << *n_long_bar_neg << " D" << dia_As_neg << '\n';
		if (delta_Al == 0.0) {
			out_file << "\ttorsion reinforcement             : " << "< nothing >" << "\n\t[ no torsion bar added to As flexural, revise(reduce) bar diameter of side face ]\n";
		}
		else if (n_flex_plus_torsion_neg == *n_long_bar_neg) {
			out_file << "\ttorsion reinforcement             : " << "< small value >\n\t[ torsion bars added in longitudinal bars ]" << '\n';
		}
		else {
			out_file << "\ttorsion reinforcement             : " << n_flex_plus_torsion_neg - *n_long_bar_neg << " D" << dia_As_neg << '\n';
		}
		out_file << "\t-----------------------------------------\n";
		out_file << "\ttherefore, negative reinforcement : " << n_flex_plus_torsion_neg << " D" << dia_As_neg << '\n';
		out_file << '\n';

		out_file << "Bar(s) at each side face :\n";
		out_file << "\tleft side face reinforcement  : " << n_bar_sideFace / 2.0 << " D" << dia_As_sideFace << '\n';
		out_file << "\tright side face reinforcement : " << n_bar_sideFace / 2.0 << " D" << dia_As_sideFace << '\n';

		out_file << '\n';
		out_file << "detailing of transverse bars :\n";
		out_file << "Alternatif 1 :\n";
		if (n_stirrup_alt_1_dummy > 2.0) {
			out_file << "\tnumber of leg-space  : " << "< nothing >" << '\n';
		}
		else {
			out_file << "\tnumber of leg-space  : " << n_stirrup_alt_1 << " D" << dia_trans << '-' << s_prov_alt_1 << '\n';
		}
		
		out_file << "Alternatif 2 :\n";
		if(s_prov_alt_1 == s_prov_alt_2)
			out_file << "\tnumber of leg-space  : " << "< nothing >" << '\n';
		else
			out_file << "\tnumber of leg-space  : " << n_stirrup_alt_2 << " D" << dia_trans << '-' << s_prov_alt_2 << '\n';
	}

	out_file << "\n                                    end of summary\n";
	out_file << "-----------------------------------------------------------------------------------\n";
	out_file << "                             End of calculation beam : " << beam_name << '\n';
	out_file << "-----------------------------------------------------------------------------------\n\n";
	
	//out_file.close();
	return std::array{ 0., 0., 0., 0., 0. };
}

//https://www.cplusplus.com/doc/tutorial/files/
//std::ios::trunc -> If the file is opened for output operations and it already existed, its previous content is deleted and replaced by the new one.
void trunc_file() {
	std::ofstream out_file;
	std::string filepath = "../input_output/output/data_output.txt";
	out_file.open(filepath, std::ios::trunc);
	if (!out_file) {
		std::cerr << "File doesn't exist!\n";
		exit(1);
	}
	out_file << "***********************************************************************************\n";
	out_file << "*              REINFORCED CONCRETE DESIGN - RECTANGULAR SECTION BEAM              *\n";
	out_file << "*                                  Version 1.0.0                                  *\n";
	out_file << "***********************************************************************************\n\n";
	//out_file.close();
	std::cout << "***********************************************************************************\n";
	std::cout << "*              REINFORCED CONCRETE DESIGN - RECTANGULAR SECTION BEAM              *\n";
	std::cout << "*                                  Version 1.0.0                                  *\n";
	std::cout << "***********************************************************************************\n\n";
}
