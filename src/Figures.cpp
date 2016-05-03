#include <iostream>
#include <cmath>
#include <fstream>

#include "Figures.h"
#include "Constants.h"
#include "Exact.h"
#include "Check.h"
#include "Hat.h"
#include "GF.h"
#include "Experiment.h"
#include "Progress.h"
#include "Tilde.h"

namespace Figures
{
// See figure 1
void Eigenvalues()
{
	double a;
	std::ofstream data[2];
	data[0].open("data/Eigenvalues_NO.txt");
	data[1].open("data/Eigenvalues_IO.txt");

	for (int i = 0; i < 2; i++)
	{
		set_ordering(i == 0);
		for (double Y = -40; Y <= 40; Y += 0.01)
		{
			a = Y * Y_to_a;
			data[i] << Y << " " << Check::l10(a) << " " << Check::l20(a) << " " << Check::l30(a) << std::endl;
		}
		data[i].close();
	}
}
// See figure 1
void Phases()
{
	std::ofstream data[2];
	data[0].open("data/Phases_NO.txt");
	data[1].open("data/Phases_IO.txt");

	double a;

	for (int i = 0; i < 2; i++)
	{
		set_ordering(i == 0);
		for (double Yrho = -40; Yrho <= 40; Yrho += 0.02)
		{
			a = Yrho * Y_to_a;
			data[i] << Yrho << " " << Hat::phi(a) << " " << Check::psi(a) << std::endl;
		} // Yrho
		data[i].close();
	} // i, ordering
}
// See figure 2
void Expansion_Parameter()
{
	std::ofstream data[2];
	data[0].open("data/Expansion_Parameter_NO.txt");
	data[1].open("data/Expansion_Parameter_IO.txt");

	for (int i = 0; i < 2; i++)
	{
		set_ordering(i == 0);
		for (double Y = -40; Y <= 40; Y += 0.05)
		{
			double a = Y * Y_to_a;

			data[i] << Y << " ";
			data[i] << std::abs(c12 * s12 * eps) << " ";
			data[i] << std::abs(epsp(a)) << " ";
			data[i] << std::endl;
		}
		data[i].close();
	}

}
// See figure 3
void Pmu2e_Precision()
{
	std::ofstream data("data/Pmu2e_Precision.txt");

	int exp_num = 0; // DUNE
	double delta = 3 * M_PI / 2;
	int alpha = 1; // nu_mu
	int beta = 0; // nu_e
	double Emin = 0.1;
	double Emax = 10;
	double Einc = 1.001;

	double Yrho, L, LE, a;
	double P_exact, P_GF0, P_GF1, P_Check2;
	set_experimental_parameters(exp_num, &Yrho, &L);

	data << L << " ";
	data << delta << " ";
	data << std::endl;

	Progress_Bar Pbar;
	for (double E = Emin; E <= Emax; E *= Einc)
	{
		LE = L / E;
		a = Yrho * E * Y_to_a;

		P_exact = Exact::Palphabeta(alpha, beta, a, LE, delta);
		P_GF0 = GF::Palphabeta(alpha, beta, a, LE, delta, 0);
		P_GF1 = GF::Palphabeta(alpha, beta, a, LE, delta, 1);
		P_Check2 = Check::Palphabeta(alpha, beta, a, LE, delta, 2);

		data << E << " ";
		data << P_exact << " ";
		data << std::abs(1 - P_GF0 / P_exact) << " ";
		data << std::abs(1 - P_GF1 / P_exact) << " ";
		data << std::abs(1 - P_Check2 / P_exact) << " ";
		data << std::endl;

		Pbar.update(Emin, Emax, E, false);
	} // E
	data.close();
}
void Eigenvalues_Bases()
{
	double a;
	std::ofstream data[2];
	data[0].open("data/Eigenvalues_Bases_NO.txt");
	data[1].open("data/Eigenvalues_Bases_IO.txt");

	for (int i = 0; i < 2; i++)
	{
		set_ordering(i == 0);
		for (double Y = -40; Y <= 40; Y += 0.01)
		{
			a = Y * Y_to_a;

			data[i] << Y << " " << Tilde::la(a) << " " << Tilde::lb(a) << " " << Tilde::lc(a) << " ";
			data[i] << Hat::lm(a) << " " << Hat::l0(a) << " " << Hat::lp(a) << " ";
			data[i] << Check::l10(a) << " " << Check::l20(a) << " " << Check::l30(a) << " ";

			data[i] << std::endl;
		}
		data[i].close();
	}
}
} // namespace Figures
int main()
{
	// Paper figures
	Figures::Eigenvalues();
	Figures::Expansion_Parameter();
	Figures::Phases();
	Figures::Pmu2e_Precision();

	// Additional figures for talks
	Figures::Eigenvalues_Bases();

	return 0;
}
