#include "material.hpp"

SteelBar& SteelBar::getInstance() {
	static SteelBar instance;
	//volatile int dummy{};
	return instance;
}

Concrete& Concrete::getInstance() {
	static Concrete instance;
	//volatile int dummy{};
	return instance;
}