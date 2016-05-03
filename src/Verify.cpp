#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

#include "Verify.h"
#include "Exact.h"
#include "Tilde.h"
#include "Constants.h"
#include "GF.h"
#include "Progress.h"
#include "Hat.h"
#include "Check.h"

namespace Verify
{
bool is_accurate(double P1, double P2, double Precision)
{
	return std::abs(1 - P1 / P2) < Precision or std::abs(P1 - P2) < Precision;
}
void Precision_Scan(double (*P_To_Verify)(int, int, double, double, double, int), int order, double Precision)
{
	double YerhoE_max = 40;
	double a_max = Y_to_a * YerhoE_max;
	double LE_max = 2500;
	int N_Steps = 10;
	std::ofstream out_file("data/Precision_Scan.txt", std::fstream::app);

	std::string flavor_names[3] = {"e ", "mu", "tau"};
	int N_Total, N_Passed;
	double P_Tilde, P_Exact;

	for (int alpha = 0; alpha <= 1; alpha++)
	{
		for (int beta = 0; beta <= alpha; beta++)
		{
			N_Total = 0;
			N_Passed = 0;
			std::cerr << flavor_names[alpha] << " -> " << flavor_names[beta] << ", ord = " << order << ", prec = " << Precision << " ..." << std::endl;
			Progress_Bar Pbar;
			for (double a = -a_max; a <= a_max; a += 2 * a_max / N_Steps)
			{
				for (double LE = -LE_max; LE <= LE_max; LE += 2 * LE_max / N_Steps)
				{
					for (double delta = 0; delta <= 2 * M_PI; delta += 2 * M_PI / N_Steps)
					{
						for (int ord = 0; ord <= 1; ord++)
						{
							set_ordering(ord == 0);
							P_Tilde = P_To_Verify(alpha, beta, a, LE, delta, order);
							P_Exact = Exact::Palphabeta(alpha, beta, a, LE, delta);
							N_Total++;
							if (is_accurate(P_Tilde, P_Exact, Precision))
								N_Passed++;
						} // ord
					} // delta
				} // LE
				Pbar.update(-a_max, a_max, a, true);
			} // a
			out_file << flavor_names[alpha] << " -> " << flavor_names[beta] << ": ";
			out_file << std::setprecision(3) << std::setw(5) << (double)N_Passed / N_Total;
			out_file << " passed, ord = " << order << ", prec = " << Precision << std::endl;
		} // beta
	} // alpha
	out_file << std::setfill('=') << std::setw(60) << "=" << std::setfill(' ') << std::endl;
	out_file.close();
}
} // namespace Verify

int main()
{
	set_Progress_Bar_visibility(true);

	std::ofstream out_file("data/Precision_Scan.txt", std::fstream::out);
	out_file << "Tilde basis:" << std::endl;
	out_file.close();
	Verify::Precision_Scan(&Tilde::Palphabeta, 0, 1e-2);
	Verify::Precision_Scan(&Tilde::Palphabeta, 1, 1e-2);

	out_file.open("data/Precision_Scan.txt", std::fstream::app);
	out_file << "Hat basis (GF):" << std::endl;
	out_file.close();
	Verify::Precision_Scan(&Hat::Palphabeta, 0, 1e-2);
	Verify::Precision_Scan(&Hat::Palphabeta, 0, 1e-3);
	Verify::Precision_Scan(&Hat::Palphabeta, 1, 1e-2);
	Verify::Precision_Scan(&Hat::Palphabeta, 1, 1e-3);

	out_file.open("data/Precision_Scan.txt", std::fstream::app);
	out_file << "Check basis (GF):" << std::endl;
	out_file.close();
	Verify::Precision_Scan(&GF::Palphabeta, 0, 1e-2);
	Verify::Precision_Scan(&Check::Palphabeta, 0, 1e-2);
	Verify::Precision_Scan(&GF::Palphabeta, 0, 1e-3);
	Verify::Precision_Scan(&GF::Palphabeta, 1, 1e-3);

	out_file.close();
	return 0;
}
