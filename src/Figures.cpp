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
#include "Mixed.h"
#include "Hyperbolas.h"

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
			data[i] << Y;
			data[i] << " " << Exact::M1sq(a) << " " << Exact::M2sq(a) << " " << Exact::M3sq(a);
//			data[i] << " " << Check::l10(a) << " " << Check::l20(a) << " " << Check::l30(a);
			data[i] << " " << Hyperbolas::l1(a) << " " << Hyperbolas::l2(a) << " " << Hyperbolas::l3(a);
			data[i] << std::endl;
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
	set_ordering(true);

	int exp_num = 0; // DUNE
	double delta = 3 * M_PI / 2;
	int alpha = 1; // nu_mu
	int beta = 0; // nu_e
	double Emin = 0.1;
	double Emax = 10;
	double Einc = 1.001;

	double Yrho, L, LE, a;
	double P_exact, P_GF0, P_GF1, P_Check2, P_Mixed, P_Mixed2, P_Hyp;
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
		P_Mixed = Mixed::Palphabeta(alpha, beta, a, LE, delta, 0);
		P_Mixed2 = Mixed::Palphabeta(alpha, beta, a, LE, delta, 2);
		P_Hyp = Hyperbolas::Palphabeta(alpha, beta, a, LE, delta);

		data << E << " ";
		data << P_exact << " ";
		data << std::abs(1 - P_GF0 / P_exact) << " ";
		data << std::abs(1 - P_GF1 / P_exact) << " ";
		data << std::abs(1 - P_Check2 / P_exact) << " ";
		data << std::abs(1 - P_Mixed / P_exact) << " ";
		data << std::abs(1 - P_Mixed2 / P_exact) << " ";
		data << std::abs(1 - P_Hyp / P_exact) << " ";
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
void Eigenvalue_Precision()
{
	double a;
	std::ofstream data[2];
	data[0].open("data/Eigenvalue_Precision_NO.txt");
	data[1].open("data/Eigenvalue_Precision_IO.txt");

	for (int i = 0; i < 1; i++) // only NO for now
	{
		set_ordering(i == 0);
		for (double Y = -40; Y <= 40; Y += 0.01)
		{
			a = Y * Y_to_a;
			data[i] << Y;
			data[i] << " " << std::min(std::abs(1 - Hyperbolas::l1(a) / Exact::M1sq(a)), std::abs(Hyperbolas::l1(a) - Exact::M1sq(a)) / Dmsqee);
			data[i] << " " << std::min(std::abs(1 - Hyperbolas::l2(a) / Exact::M2sq(a)), std::abs(Hyperbolas::l2(a) - Exact::M2sq(a)) / Dmsqee);
			data[i] << " " << std::min(std::abs(1 - Hyperbolas::l3(a) / Exact::M3sq(a)), std::abs(Hyperbolas::l3(a) - Exact::M3sq(a)) / Dmsqee);
			data[i] << std::endl;
		}
		data[i].close();
	}
}
void Reno50_Matter()
{
	std::ofstream data("data/Reno50_Matter.txt");

	int alpha = 0;		// electron type
	int beta = 0;		// electron type
	double delta = 0;	// irrelevant
	double Emin = 1.8e-3;	// 1 MeV in GeV
	double Emax = 1.2e-2;	// 10 MeV in GeV
	double L = 52.5;	// 52.5 km

	double LE;
	double Yrho = 2.6 * 0.5;
	double a;
	double P_vac, P_mat;

	for (double E = Emin; E < Emax; E += (Emax - Emin) / 2000)
	{
		data << E << " ";
		LE = L / (-E);				// -E for anti-nus
		a = Yrho * (-E) * Y_to_a;	// -E for anti-nus
		P_vac = Exact::Palphabeta(alpha, beta, 0, LE, delta);
		P_mat = Exact::Palphabeta(alpha, beta, a, LE, delta);
		data << P_mat - P_vac << std::endl;
	}
}

} // namespace Figures
int main()
{
/*
	// Paper figures
	Figures::Eigenvalues();
	Figures::Expansion_Parameter();
	Figures::Phases();
	Figures::Pmu2e_Precision();

	// Additional figures for talks
	Figures::Eigenvalues_Bases();

	// Additional figures
	Figures::Eigenvalue_Precision();
*/
	Figures::Reno50_Matter();
	return 0;
}
