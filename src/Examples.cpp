#include <iostream> // cout
#include <iomanip> // setprecision

#include "Examples.h"
#include "GF.h"
#include "Constants.h"
#include "Experiment.h"
#include "Check.h"

// Calculates P(nubar_mu->nubar_e) for DUNE parameters
void Example1()
{
	int alpha = 1;					// muon type neutrino
	int beta = 0;					// electron type neutrino
	double delta = 3 * M_PI / 2;	// the CP phase
	int exp_num = 0;				// 0 - DUNE, 1 - NOvA, 2 - T2K(HK)
	double E = 2.5;					// neutrino energy in GeV
	int order = 0;					// Zeroth order in the perturbative expansion
	E *= -1;						// Negative energy for antineutrinos

	double Yrho, L;
	set_experimental_parameters(exp_num, &Yrho, &L); // get matter density and baseline for the experiment

	double a = Yrho * E * Y_to_a;   // Y_to_a converts units
	double LE = L / E;
	double P = GF::Palphabeta(alpha, beta, a, LE, delta, order);
	std::cout << std::setprecision(10);
	std::cout << "P(nubar_mu->nubar_e) at zeroth order: ";
	std::cout << P << std::endl;

	order = 1; // First order in the perturbative expansion
	P = GF::Palphabeta(alpha, beta, a, LE, delta, order);
	std::cout << "P(nubar_mu->nubar_e) at first order:  ";
	std::cout << P << std::endl;

	order = 2; // Second order in the perturbative expansion
	P = Check::Palphabeta(alpha, beta, a, LE, delta, order);
	std::cout << "P(nubar_mu->nubar_e) at second order: ";
	std::cout << P << std::endl;
}

int main()
{
	Example1();
}
