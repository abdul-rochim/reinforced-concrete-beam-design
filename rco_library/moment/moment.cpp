#include "moment.hpp"
/*
Calculating moment capacity due to the responses of the structure analysis.

data of material such as concrete and steel reinforced strength and the dimensions are given,
compute the betha1 for supporting the calculation of the compression concrete block(a).
then compute the nominal capacity of moment(Mn).
for comparison the ultimate moment(Mu), we have to compute the reduction factor (phi).
and phi.Mn must be greater than Mu.

this moment.hpp file is using class with the new standard C++ or Modern C++
*/

inline const double Moment::initial_As_req() const {
	//the design assumption is that beams will be tension-controlled, phi = 0.9
	//this assumption will be checked later.
	using std::pow;
	double phi = double(0.9);
	double Rn = Mu * double(1000000.) / (phi * b * d * d); //convert the unit of Mu from kN.m to N.mm
	double As_required = b * d * double(0.85) * fc / fy * (1 - pow(1 - double(2) * Rn / (double(0.85) * fc), double(0.5)));
	
	return As_required;
}

inline const double Moment::a(const double& As) const {
	return As * fy / (double(0.85) * fc * b);
}

void Moment::swap(Moment& other) noexcept {
	using std::swap;
	swap(this->Mu, other.Mu);
	swap(this->fy, other.fy);
	swap(this->fc, other.fc);
	swap(this->b, other.b);
	swap(this->d, other.d);
	swap(this->betha1, other.betha1);
	swap(this->Ec, other.Ec);
}

inline double phi_moment(const Moment& m, const double& As_, double& eps_t_) {
	double eps_ty{};
	if (m.fy == 420)
		eps_ty = double(0.002);
	else
		eps_ty = m.fy / m.Es;

	double a = m.a(As_);
	double eps_t = (m.d - a / m.betha1) * m.eps_c / (a / m.betha1);

	//return epsilon t
	eps_t_ = eps_t;

	double phi{};
	if (eps_t <= eps_ty)
		phi = double(0.65);
	else if (eps_t >= eps_ty + double(0.003))
		phi = double(0.90);
	else
		phi = double(0.65) + double(0.25) * (eps_t - eps_ty) / double(0.003);
	
	return phi;
}

double phiMn(
	const Moment& m, 
	const double& dia_long, 
	const double& d_hoop, 
	const double& cover, 
	double* As_required, 
	double* nb_bar
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

	out_file << "-----------------------------------------------------------------------------------\n";
	out_file << "Starting the iteration to compute nominal moment capacity of concrete beam!\n";
	out_file << "-----------------------------------------------------------------------------------\n";
	out_file << std::fixed << std::setprecision(2);

	out_file << "initial As longitudinal required : " << m.initial_As_req() << " mm2\n";

	//compute number of bar reinforcement required
	double n_ = m.initial_As_req() / (double(0.25) * M_PI * dia_long * dia_long);
	out_file << "initial number of longitudinal bars required : " << n_ << " bars\n";
	if(n_ > 0.)
		out_file << "initial number of longitudinal bars required(roundup) : " << std::ceil(n_) << " bars\n";
	
	out_file << "initial number of bars provided : " << *nb_bar << " bars\n";

	if (*nb_bar >= n_) {
		out_file << "if initial number of bars provided : " << *nb_bar << " bars > initial number of bars required : " << n_ << " bars\n";
		out_file << "assumption, initial number of bars required is using number of bars provided : " << *nb_bar << " bars\n";
		n_ = *nb_bar;
	}
	else if (n_ <= 1) {
		out_file << "if initial number of bars required : " << n_ << " bars <= 1\n";
		out_file << "initial number of bars required used : " << 2.0 << " bars\n";
		n_ = 2.0; //minimum number of bar reinforcement placed at both corner
	}
	out_file << '\n';

	//declare variables of As_req_iteration, As_prov, phiMn
	double As_req_iteration;
	double As_prov;
	double phiMn;

	//minimum clear spacing between the horizontal
	double clear_spacing_long = std::max(double(25.), dia_long);

	//declare of number of bars per layer
	double n_dummy;
	double n1{};
	double n2{};
	double n3{};
	//double n4{};

	//declare clear spacing provided
	double clear_spacing_provided;

	//compute maximum bar spacing at tension face must not exceed the lesser of:
	double fs_max_spacing = 2. / 3 * m.fy; //MPa
	double Cc_ = cover + d_hoop; //Cc is the cover to the longitudinal bar, not to the tie.
	double max_bar_spacing1 = 380. * (280. / fs_max_spacing) - 2.5 * Cc_; //SP-17M(14)
	double max_bar_spacing2 = 300. * (280. / fs_max_spacing); //SP-17M(14)
	double max_bar_spacing = std::min(max_bar_spacing1, max_bar_spacing2);

	//declare As provided
	double As1;
	double As2;
	double As3;
	//double As4;

	//declare effective depth(d)
	double bar_spacing_layer = std::max(double(40.), clear_spacing_long+dia_long); //maximum(40, clear spacing+db)
	double d1 = m.d;
	double d2{};
	double d3{};
	//double d4{};

	//declare epsilon t, a depth of compressive concrete, phi reduction factor
	double eps_t{};
	double a_;
	double phi_;

	//declare variables of As_minimum
	double As_min1;
	double As_min2;
	double As_minimum;

	//declare variables of rho balance, rho maximum
	double rho_balance;
	double rho_max_per_rho_balance;
	double rho_maximum;
	double As_maximum;

	//iteration for computing phi Mn dan check to the minimum longitudinal bars
	unsigned int i = 1;
	do {
		out_file << "iteration_" << i << '\n';
		out_file << "n longitudinal required from the iteration : " << n_ << " bars\n";
		out_file << "n longitudinal provided from the iteration : " << std::ceil(n_) << " bars\n";
		As_req_iteration = n_ * double(0.25) * M_PI * dia_long * dia_long;

		//assigment of bars area requirement
		*As_required = As_req_iteration;

		As_prov = std::ceil(n_) * double(0.25) * M_PI * dia_long * dia_long; //round up from As required or As required iteration
		out_file << "As longitudinal required from the iteration : " << As_req_iteration << " mm2\n";
		out_file << "As longitudinal provided from the iteration : " << As_prov << " mm2\n";

		//beam's web required, bw_req = 2 (cover + d_stirrup + db) + (n-1) * db + (n-1) * clear spacing
		//assume n_dummy = n - 1
		n_dummy = (m.b - double(2) * (cover + d_hoop + dia_long)) / (dia_long + clear_spacing_long);
		//re-calculate n_dummy;
		n_dummy = n_dummy + 1.0;
		out_file << "maximum number of bars per layer from the calculation : " << n_dummy << " bars per layer\n";
		out_file << "maximum number of bars per layer recommended : " << std::trunc(n_dummy) << " bars per layer\n";

		//compute the distribution of bar layer
		if (std::ceil(n_) <= std::trunc(n_dummy)) {
			n1 = std::ceil(n_);
		}
/*
		///////////////// for n1,n2,n3,n4 /////////////////
		else if (std::ceil(n_) > 3.0 * std::trunc(n_dummy)) {
			if ((std::ceil(n_) - 3.0 * std::trunc(n_dummy)) == 1.0) {
				n4 = 2.0;
				n3 = std::trunc((std::ceil(n_) - n4) / 3.0);
				n2 = std::trunc((std::ceil(n_) - n4 - n3) / 2.0);
				n1 = std::ceil(n_) - n4 - n3 - n2;
			}
			else {
				n1 = n2 = n3 = std::trunc(n_dummy);
				n4 = std::ceil(n_) - n1 - n2 - n3;
			}
		}
		///////////////////////////////////////////////////
*/
		else if (std::ceil(n_) > 2.0 * std::trunc(n_dummy)) {
			//asumption number leg = 1 leg is converted 2 legs
		/*	if ((std::ceil(n_) - 2.0 * std::trunc(n_dummy)) == 1.0) {
				n3 = 2.0;
				n2 = std::trunc((std::ceil(n_) - n3) / 2.0);
				n1 = std::ceil(n_) - n3 - n2;
			}
			else {
				n1 = n2 = std::trunc(n_dummy);
				n3 = std::ceil(n_) - n1 - n2;
			}*/
			//asumption number leg = 1 leg is NOT converted to be 2 legs
			n1 = n2 = std::trunc(n_dummy);
			n3 = std::ceil(n_) - n1 - n2;
		}
		else {
			//asumption number leg = 1 leg is converted 2 legs
		/*	if ((std::ceil(n_) - std::trunc(n_dummy)) == 1.0) {
				n2 = 2.0;
				n1 = std::ceil(n_) - n2;
			}
			else {
				n1 = std::trunc(n_dummy);
				n2 = std::ceil(n_) - n1;
			}*/
			//asumption number leg = 1 leg is NOT converted to be 2 legs
			n1 = std::trunc(n_dummy);
			n2 = std::ceil(n_) - n1;
		}
		out_file << "number of bars, layer-1 : " << n1 << " bars" <<
			"\nnumber of bars, layer-2 : " << n2 << " bars" <<
			"\nnumber of bars, layer-3 : " << n3 << " bars\n";// <<
//			"\nnumber of bars, layer-4 : " << n4 << " bars\n";
		
		//display minimum and maximum clear spacing bar 
		out_file << "minimum clear spacing bar : " << clear_spacing_long << " mm\n";
		out_file << "maximum bar spacing : " << max_bar_spacing << " mm\n";

		//compute clear spacing bar provided
		//clear_spacing_provided = (m.b - double(2) * (cover + d_hoop + dia_long) - (std::max(n1, std::max(n2, std::max(n3, n4))) - 1.0) * dia_long) / (std::max(n1, std::max(n2, std::max(n3, n4))) - 1.0);
		clear_spacing_provided = (m.b - double(2) * (cover + d_hoop + dia_long) - (std::max(n1, std::max(n2, n3)) - 1.0) * dia_long) / (std::max(n1, std::max(n2, n3)) - 1.0);
		out_file << "clear spacing bar provided : " << clear_spacing_provided << " mm\n";
		
		//check maximum bar spacing!
		//this limit is intended to control flexural cracking width.
		try {
			if (clear_spacing_provided < max_bar_spacing)
				out_file << "clear spacing bar provided : " << clear_spacing_provided << " mm < max bar spacing : " << max_bar_spacing << " mm\n";
			else {
				std::cout << "clear spacing bar provided : " << clear_spacing_provided << " mm > max bar spacing : " << max_bar_spacing << " mm\n";
				std::cout << "need to add longitudinal bar(s)\n";
				throw std::invalid_argument("WARNING! clear spacing bar provided is greater than maximum bar spacing!");
			}
		}
		catch (std::invalid_argument& e) {
			std::cerr << e.what() << '\n';
			throw;
			//return -1;
		}
		//check minimum bar spacing!
		try {
			if (clear_spacing_provided > clear_spacing_long)
				out_file << "clear spacing bar provided : " << clear_spacing_provided << " mm > min bar spacing : " << clear_spacing_long << " mm\n";
			else {
				std::cout << "clear spacing bar provided : " << clear_spacing_provided << " mm < min bar spacing : " << clear_spacing_long << " mm\n";
				std::cout << "enlarge the cross-sectional dimensions\n";
				throw std::invalid_argument("WARNING! clear spacing bar provided is lesser than minimum bar spacing!");
			}		
		}
		catch (std::invalid_argument& e) {
			std::cerr << e.what() << '\n';
			throw;
			//return -1;
		}

		//number of bars assignment
		*nb_bar = n1 + n2 + n3;// +n4;

		As1 = n1 * double(0.25) * M_PI * dia_long * dia_long;
		As2 = n2 * double(0.25) * M_PI * dia_long * dia_long;
		As3 = n3 * double(0.25) * M_PI * dia_long * dia_long;
//		As4 = n4 * double(0.25) * M_PI * dia_long * dia_long;

		//std::cout << "As layer-1 provided : " << As1 << '\n';
		//std::cout << "As layer-2 provided : " << As2 << '\n';
		//std::cout << "As layer-3 provided : " << As3 << '\n';
		//std::cout << "As layer-4 provided : " << As4 << '\n';
		
		/*
		if (n2 > 0.0 && n3 == 0.0 && n4 == 0.0)
			d2 = m.d - bar_spacing_layer;
		else if (n3 > 0.0 && n4 == 0.0) {
			d2 = m.d - bar_spacing_layer;
			d3 = m.d - 2.0 * bar_spacing_layer;
		}
		else if (n4 > 0.0) {
			d2 = m.d - bar_spacing_layer;
			d3 = m.d - 2.0 * bar_spacing_layer;
			d4 = m.d - 3.0 * bar_spacing_layer;
		}
		*/
		if (n2 > 0.0 && n3 == 0.0)
			d2 = m.d - bar_spacing_layer;
		else if (n3 > 0.0) {
			d2 = m.d - bar_spacing_layer;
			d3 = m.d - 2.0 * bar_spacing_layer;
		}

		//calculation or assignment of epsilon t, a, phi
		a_ = m.a(As_prov);
		phi_ = phi_moment(m, As_prov, eps_t);
		out_file << std::fixed << std::setprecision(4);
		out_file << "\nepsilon t calc' : " << eps_t << '\n'; //epsilon t for computing rho maximum per rho balance
		out_file << std::fixed << std::setprecision(2);

		//compute nominal moment capacity
//		phiMn = phi_ * (As1 * m.fy * d1 + As2 * m.fy * d2 + As3 * m.fy * d3 + As4 * m.fy * d4 - (As1 + As2 + As3 + As4) * m.fy * a_ / 2.0) / 1000000.0; //kN.m
		phiMn = phi_ * (As1 * m.fy * d1 + As2 * m.fy * d2 + As3 * m.fy * d3 - (As1 + As2 + As3) * m.fy * a_ / 2.0) / 1000000.0; //kN.m

		//Checking the reinforcement ratio
		//minimum reinforcement ratio
		out_file << "minimum and maximum reinforcement ratio:\n";
		As_min1 = 0.25 * std::pow(m.fc, 0.5) / m.fy * m.b * m.d;
		As_min2 = 1.4 / m.fy * m.b * m.d;
		As_minimum = std::max(As_min1, As_min2);
		out_file << "As minimum requirement 1 [0.25 sqrt(fc')/fy *bd] : " << As_min1 << " mm2\n";
		out_file << "As minimum requirement 2 [1.4/fy *bd] : " << As_min2 << " mm2\n";
		out_file << "As minimum requirement [used] : " << As_minimum << " mm2\n";

		//maximum reinforcement ratio
		//calculate epsilon t for defining rho max per rho balance
		//if eps_t >= 0.005 use eps_t = 0.005, if eps_t < 0.005 use eps_t = 0.004
		if (eps_t >= 0.005) {
			eps_t = 0.005;
		}
		else {
			//epsilon t not less than 0.004 due to ensure the ductility level and to show
			//the visually visible signs before collapse.
			eps_t = 0.004;
		}
		//assignment of rho balance, and rho maximum
		rho_balance = 0.85 * m.betha1 * m.fc / m.fy * (600 / (600 + m.fy));	//rho balance
		rho_max_per_rho_balance = (0.003 + m.fy / m.Es) / (0.003 + eps_t);	//rho maximum per rho balance
		rho_maximum = rho_max_per_rho_balance * rho_balance;	//rho maximum
		out_file << std::fixed << std::setprecision(4);
		out_file << "epsilon t for computing rho maximum: " << eps_t << '\n';
		out_file << "rho maximum per rho balance : " << rho_max_per_rho_balance << '\n';
		out_file << "rho maximum : " << rho_maximum << '\n';
		As_maximum = rho_maximum * m.b * m.d;
		out_file << std::fixed << std::setprecision(2);
		out_file << "As maximum requirement : " << As_maximum << " mm2\n";

		//error handling
		try {
			if (As_prov >= As_minimum && As_prov <= As_maximum) {
				out_file << "Ok, As minimum : " << As_minimum << " mm2 < " << "As longitudinal reinforcement provided : " << As_prov << " mm2 < " << "As maximum : " << As_maximum << " mm2\n";
			}
			else if (As_prov < As_minimum && phiMn >= 4.0 / 3.0 * m.Mu) {
				out_file << "As longitudinal reinforcement provided : " << As_prov << " < As minimum : " << As_minimum << " and phiMn : " << phiMn << " kN.m > (4/3) Mu : " << 4. / 3 * m.Mu << " kN.m\n";
				out_file << "Ok, the element structure is big and massive\n";
			}
			else if (As_prov < As_minimum)
			{
				out_file << "As_provided : " << As_prov << " mm2 < As_minimum : " << As_minimum << " mm2\n";
				out_file << "reinforcement ratio provided is lesser then minimum reinforcement ratio\n";
			}
			else //if(As_prov > As_maximum)
			{
				std::cout << "As_provided : " << As_prov << " mm2 > As_maximum : " << As_maximum << " mm2\n";
				std::cout << "reinforcement ratio provided is greater then maximum reinforcement ratio\n";
				std::cout << "The beam section or material properties needs to be revised, \nenlarge the cross-sectional dimensions or revise the material properties!\n";
				throw std::invalid_argument("WARNING! As provided > As maximum");
			}
		}
		catch (std::invalid_argument& e){
			std::cerr << e.what() << '\n';
			throw;
			//return -1;
		}

		//check assumption of phi_moment
		out_file << "\ncheck the assumption of phi!\n";
		if (phi_ == 0.9)
			out_file << "the design assumption of using phi is 0.9 in tension-controlled is verified, is phiMn > Mu ?\n";
		else {
			out_file << "the phi calculation is " << phi_ << " , is phiMn > Mu ?\n";
		}
		
		//check is phi Mn > Mu and As provided?
		if (phiMn >= 4./3 * m.Mu && As_prov < As_minimum) {
			out_file << "whereas As_provided: " << As_prov << " mm2 << As_minimum: " << As_minimum << " mm2\n";
			out_file << "but it's Ok, since phiMn: " << phiMn << " kN.m > 4/3 Mu: " << 4./3 * m.Mu << " kN.m ";
			out_file << "[massive structure]\n";
		}
		else if (phiMn >= m.Mu && As_prov >= As_minimum) {
			out_file << "Ok, phiMn: " << phiMn << " kN.m > Mu: " << m.Mu << " kN.m ";
			out_file << "and As_provided: " << As_prov << " mm2 > As_minimum: " << As_minimum << " mm2\n";
		}
		else if (phiMn < m.Mu && As_prov > As_minimum) {
			out_file << "As_provided: " << As_prov << " mm2 > As_minimum: " << As_minimum << " mm2\n";
			out_file << "but phiMn: " << phiMn << " kN.m < Mu: " << m.Mu << " kN.m, re-calculate again!\n";
		}
		else {
			out_file << "phiMn: " << phiMn << " kN.m > Mu: " << m.Mu << " kN.m\n";
			out_file << "but As_provided: " << As_prov << " mm2 < As_minimum: " << As_minimum << " mm2, re-calculate again!\n";
		}

		out_file << '\n';
		++n_;
		++i;
	} while ((phiMn < m.Mu && As_prov < As_minimum) || phiMn < m.Mu || (phiMn < 4./3. * m.Mu && As_prov < As_minimum));
	
	out_file << "calculation of the nominal moment capacity complete!\n";
	out_file << "-----------------------------------------------------------------------------------\n";
	//out_file.close();
	return phiMn;
}
