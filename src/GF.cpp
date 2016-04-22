/*
This file contains the transition probabilities calculated to zeroth or first order using the general form presented in tables 1 and 2.
*/
#include <cassert>
#include <complex>
#include <cmath>

#include "GF.h"
#include "Constants.h"
#include "Check.h"
#include "Hat.h"
#include "Warning.h"

namespace GF
{
double ReUm4(double a, double delta, int alpha, int beta, int w, int x, int y, int z)
{
	assert (alpha >= 0 && alpha <= 2); // alpha valid for 0, 1, 2
	assert (beta >= 0 && beta <= 2); // beta valid for 0, 1, 2
	assert (w >= 0 && w <= 2); // w valid for 0, 1, 2
	assert (x >= 0 && x <= 2); // x valid for 0, 1, 2
	assert (y >= 0 && y <= 2); // y valid for 0, 1, 2
	assert (z >= 0 && z <= 2); // z valid for 0, 1, 2

	Matrix<std::complex<double> > U = Check::UMNSm(a, delta);

	return std::real(U(alpha, w) * std::conj(U(beta, x) * U(alpha, y)) * U(beta, z));
}
double F1(int alpha, int Beta, double a, double delta)
{
	double spsi = sin(Check::psi(a));
	return -spsi * (ReUm4(a, delta, alpha, Beta, 2, 0, 1, 1) + ReUm4(a, delta, alpha, Beta, 0, 2, 1, 1));
}
double F2(int alpha, int Beta, double a, double delta)
{
	double cpsi = cos(Check::psi(a));
	return cpsi * (ReUm4(a, delta, alpha, Beta, 2, 1, 0, 0) + ReUm4(a, delta, alpha, Beta, 1, 2, 0, 0));
}
double G1(int alpha, int Beta, double a, double delta)
{
	double spsi = sin(Check::psi(a));
	Matrix<std::complex<double> > U = Check::UMNSm(a, delta);
	double deltaalphabeta = alpha == Beta ? 1 : 0;
	return -spsi * std::real((U(alpha, 0) * std::conj(U(Beta, 2)) + U(alpha, 2) * std::conj(U(Beta, 0)))
				* (2. * std::conj(U(alpha, 2)) * U(Beta, 2) - deltaalphabeta));
}
double G2(int alpha, int Beta, double a, double delta)
{
	double cpsi = cos(Check::psi(a));
	Matrix<std::complex<double> > U = Check::UMNSm(a, delta);
	double deltaalphabeta = alpha == Beta ? 1 : 0;
	return cpsi * std::real((U(alpha, 1) * std::conj(U(Beta, 2)) + U(alpha, 2) * std::conj(U(Beta, 1)))
				* (2. * std::conj(U(alpha, 2)) * U(Beta, 2) - deltaalphabeta));
}
double C21(int alpha, int Beta, double a, double delta, int order)
{
	assert (order >= 0 && order <= 1); // order valid for 0, 1
	double C210 = -ReUm4(a, delta, alpha, Beta, 1, 1, 0, 0);
	double C1 = epsp(a) * Dmsqee;
	double F1 = GF::F1(alpha, Beta, a, delta);
	double F2 = GF::F2(alpha, Beta, a, delta);
	double Dl31 = Check::Dl31(a, order);
	double Dl32 = Check::Dl32(a, order);
	double C211 = C1 * (F1 / Dl31 + F2 / Dl32);
	return C210 + order * C211;
}
double C31(int alpha, int Beta, double a, double delta, int order)
{
	assert (order >= 0 && order <= 1); // order valid for 0, 1
	double C310 = -ReUm4(a, delta, alpha, Beta, 2, 2, 0, 0);
	double C1 = epsp(a) * Dmsqee;
	double F1 = GF::F1(alpha, Beta, a, delta);
	double F2 = GF::F2(alpha, Beta, a, delta);
	double G1 = GF::G1(alpha, Beta, a, delta);
	double Dl31 = Check::Dl31(a, order);
	double Dl32 = Check::Dl32(a, order);
	double C311 = C1 * ((F1 + G1) / Dl31 - F2 / Dl32);
	return C310 + order * C311;
}
double C32(int alpha, int Beta, double a, double delta, int order)
{
	assert (order >= 0 && order <= 1); // order valid for 0, 1
	double C320 = -ReUm4(a, delta, alpha, Beta, 2, 2, 1, 1);
	double C1 = epsp(a) * Dmsqee;
	double F1 = GF::F1(alpha, Beta, a, delta);
	double F2 = GF::F2(alpha, Beta, a, delta);
	double G2 = GF::G2(alpha, Beta, a, delta);
	double Dl31 = Check::Dl31(a, order);
	double Dl32 = Check::Dl32(a, order);
	double C321 = C1 * (-F1 / Dl31 + (F2 + G2) / Dl32);
	return C320 + order * C321;
}
double D(int alpha, int Beta, double a, double delta, int order)
{
	assert (order >= 0 && order <= 1); // order valid for 0, 1
	if (alpha == Beta)
		return 0;
	double sdelta = sin(delta);
	double phi = Hat::phi(a);
	double cphi = cos(phi);
	double sphi = sin(phi);
	double cphisq = pow(cphi, 2);
	double sphisq = pow(sphi, 2);
	double psi = Check::psi(a);
	double cpsi = cos(psi);
	double spsi = sin(psi);
	double cpsisq = pow(cpsi, 2);
	double spsisq = pow(spsi, 2);
	double D0 = sdelta * Check::Jrm(a);
	double D1 = sdelta * c23 * s23 * cphi * (epsp(a) * Dmsqee * pow(spsi, 2) * (sphisq - cphisq * cpsisq) / Check::Dl31(a)
				- epsp(a) * Dmsqee * pow(cpsi, 2) * (sphisq - cphisq * spsisq) / Check::Dl32(a));
	double D = D0 + order * D1;
	if ((alpha + Beta) % 2 == 0)
		D *= -1;
	if (alpha > Beta)
		D *= -1;
	return D;
}
double Palphabeta(int alpha, int beta, double a, double LE, double delta, int order)
{
	assert (alpha >= 0 && alpha <= 2); // alpha valid for 0, 1, 2
	assert (beta >= 0 && beta <= 2); // Beta valid for 0, 1, 2
	assert (order >= 0 && order <= 1); // order valid for 0, 1

	double L4E = km_per_GeV_to_per_eV2 * LE / 4; // in eV^-2 now

	double Dl21 = Check::Dl21(a, order);
	double Dl31 = Check::Dl31(a, order);
	double Dl32 = Check::Dl32(a, order);

	double Delta21 = Dl21 * L4E;
	double Delta31 = Dl31 * L4E;
	double Delta32 = Dl32 * L4E;

	double deltaalphabeta = alpha == beta ? 1 : 0;

	double C21 = GF::C21(alpha, beta, a, delta, order);
	double C31 = GF::C31(alpha, beta, a, delta, order);
	double C32 = GF::C32(alpha, beta, a, delta, order);
	double D = GF::D(alpha, beta, a, delta, order);

	double P = deltaalphabeta + 4 * C21 * pow(sin(Delta21), 2) + 4 * C31 * pow(sin(Delta31), 2) + 4 * C32 * pow(sin(Delta32), 2)
			+ 8 * D * sin(Delta21) * sin(Delta31) * sin(Delta32);
	ProbabilityWarning W(&P);
	return P;
}
} // GF
