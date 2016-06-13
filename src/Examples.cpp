#include <iostream> // cout
#include <iomanip> // setprecision

#include "Examples.h"
#include "GF.h"
#include "Constants.h"
#include "Experiment.h"
#include "Check.h"
#include "Exact.h"

// Calculates P(nubar_mu->nubar_e) for DUNE parameters
void Example1()
{
	std::cout << "==== Example 1 ====" << std::endl;

	int alpha = 1;					// muon type neutrino
	int beta = 0;					// electron type neutrino
	double delta = 3 * M_PI / 2;	// The CP phase
	int exp_num = 0;				// 0 - DUNE, 1 - NOvA, 2 - T2K/T2HK
	double E = 2.5;					// neutrino energy in GeV
	int order = 0;					// Zeroth order in the perturbative expansion
	E *= -1;						// Negative energy for antineutrinos

	double Yrho, L;
	set_experimental_parameters(exp_num, &Yrho, &L); // Get matter density and baseline for the experiment

	double a = Yrho * E * Y_to_a;   // Y_to_a converts units
	double LE = L / E;
	double P = GF::Palphabeta(alpha, beta, a, LE, delta, order); // Calculate the transition probability using the general form described in the tables

	std::cout << std::setprecision(10); // Increases the number of decimal places to output
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

void Example2()
{
	std::cout << "==== Example 2 ====" << std::endl;

	int alpha = 0;		// electron type
	int beta = 0;		// electron type
	double delta = 0;	// irrelevant
	double E = 2.5e-3;	// 2.5 MeV in GeV
	E *= -1;			// anti-nus
	double L = 52.5;	// 52.5 km

	double LE = L / E;
	double Yrho = 2.6 * 0.5;
	double a = Yrho * E * Y_to_a;
	double P_vac = Exact::Palphabeta(alpha, beta, 0, LE, delta);
	double P_mat = Exact::Palphabeta(alpha, beta, a, LE, delta);

	std::cout << P_mat - P_vac << std::endl;

}

int main()
{
	Example1();
	Example2();
}
