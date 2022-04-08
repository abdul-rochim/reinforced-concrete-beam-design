#ifndef MATERIAL_HPP
#define MATERIAL_HPP

//Singleton
//modernescpp.com -> Rainer Grimm
struct SteelBar {
	static SteelBar& getInstance();

	constexpr double fy400() const { return 400.f; }
	constexpr double fy420() const { return 420.f; }
	constexpr double fy450() const { return 450.f; }
	constexpr double fy500() const { return 500.f; }
	constexpr double fy550() const { return 550.f; }
	constexpr double fy600() const { return 600.f; }

private:
	SteelBar() = default;
	~SteelBar() = default;
	SteelBar(const SteelBar&) = delete;
	SteelBar& operator = (const SteelBar&) = delete;
	//SteelBar(SteelBar&&) = delete;
	//SteelBar& operator = (SteelBar&&) = delete;
};

struct Concrete {
	static Concrete& getInstance();

	constexpr double fc20() const { return 20.f; }
	constexpr double fc25() const { return 25.f; }
	constexpr double fc30() const { return 30.f; }
	constexpr double fc35() const { return 35.f; }
	constexpr double fc40() const { return 40.f; }
	constexpr double fc45() const { return 45.f; }
	constexpr double fc50() const { return 50.f; }
	constexpr double fc55() const { return 55.f; }
	constexpr double fc60() const { return 60.f; }

private:
	Concrete() = default;
	~Concrete() = default;
	Concrete(const Concrete&) = delete;
	Concrete& operator = (const Concrete&) = delete;
	//Concrete(Concrete&&) = delete;
	//Concrete& operator = (Concrete&&) = delete;
};

#endif