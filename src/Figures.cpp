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
#include "LWX.h"
#include "First_Peak.h"

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
			data[i] << " " << Check::l10(a) << " " << Check::l20(a) << " " << Check::l30(a);
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
	double P_exact, P_GF0, P_GF1, P_Check2, P_Mixed, P_Mixed2;
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

		data << E << " ";
		data << P_exact << " ";
		data << std::abs(1 - P_GF0 / P_exact) << " ";
		data << std::abs(1 - P_GF1 / P_exact) << " ";
		data << std::abs(1 - P_Check2 / P_exact) << " ";
		data << std::abs(1 - P_Mixed / P_exact) << " ";
		data << std::abs(1 - P_Mixed2 / P_exact) << " ";
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

	double l10, l20, l30;
	double l12, l22, l32;
	double M1sq, M2sq, M3sq;
	int order;
	for (int i = 0; i < 1; i++) // only NO for now
	{
		set_ordering(i == 0);
		for (double Y = -40; Y <= 40; Y += 0.011)
		{
			a = Y * Y_to_a;

			order = 0;
			l10 = Check::l1(a, order);
			l20 = Check::l2(a, order);
			l30 = Check::l3(a, order);

			order = 2;
			l12 = Check::l1(a, order);
			l22 = Check::l2(a, order);
			l32 = Check::l3(a, order);

			M1sq = Exact::M1sq(a);
			M2sq = Exact::M2sq(a);
			M3sq = Exact::M3sq(a);

			data[i] << Y;
			data[i] << " " << std::abs(1 - l10 / M1sq);
			data[i] << " " << std::abs(1 - l20 / M2sq);
			data[i] << " " << std::abs(1 - l30 / M3sq);
			data[i] << " " << std::abs(1 - l12 / M1sq);
			data[i] << " " << std::abs(1 - l22 / M2sq);
			data[i] << " " << std::abs(1 - l32 / M3sq);
			data[i] << std::endl;
		}
		data[i].close();
	}
}
void Reno50_Matter()
{
	std::ofstream data[2];
	data[0].open("data/Reno50_Matter_NO.txt");
	data[1].open("data/Reno50_Matter_IO.txt");

	int alpha = 0;		// electron type
	int beta = 0;		// electron type
	double delta = 0;	// irrelevant
	double Emin = 1.8e-3;	// 1 MeV in GeV
	double Emax = 1.2e-2;	// 10 MeV in GeV
	double L = 52.5;	// 52.5 km

	double LE;
	double Yrho = 2.6 * 0.5;
	double a;
	double P_vac, P_mat, P0_vac, P0_mat, PLWX_vac, PLWX_mat;

	for (int i = 0; i < 2; i++)
	{
		set_ordering(i == 0);
		for (double E = Emin; E < Emax; E += (Emax - Emin) / 2000)
		{
			LE = L / E;
			a = Yrho * E * Y_to_a;

			P_vac = Exact::Palphabeta(alpha, beta, 0, -LE, delta);	// -LE for anti-nus
			P_mat = Exact::Palphabeta(alpha, beta, -a, -LE, delta);	// -a, -LE for anti-nus

			P0_vac = GF::Palphabeta(alpha, beta, 0, -LE, delta, 0);	// -LE for anti-nus
			P0_mat = GF::Palphabeta(alpha, beta, -a, -LE, delta, 0);	// -a, -LE for anti-nus

			PLWX_vac = LWX::Pee(0, LE);	// formulas are already for anti-nus
			PLWX_mat = LWX::Pee(a, LE);	// formulas are already for anti-nus

			data[i] << E << " ";
			data[i] << P_mat - P_vac << " ";
			data[i] << P0_mat - P0_vac << " ";
			data[i] << PLWX_mat - PLWX_vac << " ";
			data[i] << std::endl;
		} // E
		data[i].close();
	} // i, ordering
}
void FP_Line()
{
	std::ofstream data[2];
	data[0].open("data/FP_Line_NO.txt");
	data[1].open("data/FP_Line_IO.txt");

	int alpha = 1;
	int beta = 0;
	double Yrho = 1.4;
	Yrho /= 2;
	bool max = true;
	double m, b;

	for (int i = 0; i < 2; i++)
	{
		set_ordering(i == 0);
		for (double delta = 0; delta <= 2 * M_PI; delta += 2 * M_PI / 40)
		{
			FP_line(alpha, beta, Yrho, delta, max, &m, &b);
			data[i] << delta << " " << m << " " << b << std::endl;
		} // delta
		data[i].close();
	} // i, ordering
}
void FP_Vac()
{
	std::ofstream data;
	data.open("data/FP_Vac.txt");

	int alpha = 1;
	int beta = 0;
	double Yrho = 0;
	double E = 1;
	bool max = true;

	double delta, Le, Ln, L0; // Le = exact L, Ln = nth order L, L0 = zeroth order L
	int N = 1e2;
	for (int i = 0; i < N; i++)
	{
		delta = 2 * M_PI * i / N;
		data << delta << " ";
		for (int j = 0; j < 3; j++)
		{
			Le = First_Extremum_Exact(alpha, beta, E, Yrho, delta, max);
			
			Ln = LE_SP1(delta, Yrho * E * Y_to_a, j);
			L0 = LE_SP1(delta, Yrho * E * Y_to_a, 0);
			data << (Le / E - Ln) / L0 << " ";
		} // j, order
		data << std::endl;
	} // i, delta
	data.close();
}
void FP_Mat()
{
	std::ofstream data;
	data.open("data/FP_Mat.txt");

	int alpha = 1;
	int beta = 0;
	double Yrho = 1.4;
	double E = 2.5;
	bool max = true;

	double delta, Le, Ln, L0; // Le = exact L, Ln = nth order L, L0 = zeroth order L
	int N = 1e2;
	for (int i = 0; i < N; i++)
	{
		delta = 2 * M_PI * i / N;
		data << delta << " ";
		for (int j = 0; j < 3; j++)
		{
			Le = First_Extremum_Exact(alpha, beta, E, Yrho, delta, max);
			
			Ln = LE_SP2(delta, Yrho * E * Y_to_a, j);
			L0 = LE_SP2(delta, Yrho * E * Y_to_a, 0);
			data << (Le / E - Ln) / L0 << " ";
		} // j, order
		data << std::endl;
	} // i, delta
	data.close();
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
	Figures::Reno50_Matter();
	Figures::Eigenvalue_Precision();
	Figures::FP_Line();
*/
	Figures::FP_Vac();
	Figures::FP_Mat();
	return 0;
}
