/*
From Zaglauer and Schwarzer: Z. Phys. C40 (1988) 273.
Notes:
	WLOG, take m_1=0, as can be seen by eqs 12, 13
Corrections:
	Eq 24 for S should be ^1.5 in the denominator, not ^(1/3)
	Eq 50 for deltam should be s23c23 not s23c23^2 in the numerator
*/

#include <cmath>
#include <complex>
#include <cassert>

#include "Exact.h"
#include "Constants.h"
#include "Matrix.h"

namespace Exact
{
// nominal definition, true for normal ordering
double M1sqn(double D);
double M2sqn(double D);
double M3sqn(double D);
// actual definitions. For NO, the same as the nominal. For IO, 1=2n, 2=3n, 3=1n
double M1sq(double D);
double M2sq(double D);
double M3sq(double D);

double A(double D);
double B(double D);
double C(double D);
double S(double D);

double alpha();
double beta();

double E(double D);
double F(double D);

double s12msq(double D);
double c12msq(double D);
double s13msq(double D);
double c13msq(double D);
double s23msq(double D, double delta);
double c23msq(double D, double delta);
double cdeltam(double D, double delta);
double sdeltam(double D, double delta);

// See eq. 18
double M1sqn(double D)
{
	double x = sqrt(pow(A(D), 2) - 3 * B(D));
	return (A(D) / 3.) - x * S(D) / 3 - sqrt(3) * x * sqrt(1 - pow(S(D), 2)) / 3;
}
// See eq. 19
double M2sqn(double D)
{
	double x = sqrt(pow(A(D), 2) - 3 * B(D));
	return (A(D) / 3.) - x * S(D) / 3 + sqrt(3) * x * sqrt(1 - pow(S(D), 2)) / 3;
}
// See eq. 20
double M3sqn(double D)
{
	double x = sqrt(pow(A(D), 2) - 3 * B(D));
	return (A(D) / 3.) + 2 * x * S(D) / 3;
}
double M1sq(double D)
{
	return normal ? M1sqn(D) : M2sqn(D);
}
double M2sq(double D)
{
	return normal ? M2sqn(D) : M3sqn(D);
}
double M3sq(double D)
{
	return normal ? M3sqn(D) : M1sqn(D);
}
double DM21sq(double D)
{
	return M2sq(D) - M1sq(D);
}
double DM31sq(double D)
{
	return M3sq(D) - M1sq(D);
}
double DM32sq(double D)
{
	return M3sq(D) - M2sq(D);
}
// See eq. 21
double A(double D)
{
	return Dmsq21 + Dmsq31 + D;
}
// See eq. 22
double B(double D)
{
	return Dmsq31 * Dmsq21 + D * (Dmsq31 * c13sq + Dmsq21 * (c13sq * c12sq + s13sq));
}
// See eq. 23
double C(double D)
{
	return D * Dmsq31 * Dmsq21 * c13sq * c12sq;
}
// See eq. 24
double S(double D)
{
	double x = 2 * pow(A(D), 3) - 9 * A(D) * B(D) + 27 * C(D);
	double y = 2 * pow(pow(A(D), 2) - 3 * B(D), 1.5);
	return cos((1. / 3) * acos(x / y));
}
// See eq. 43
double alpha()
{
	return Dmsq31 * c13sq + Dmsq21 * (c13sq * c12sq + s13sq);
}
// See eq. 44
double beta()
{
	return Dmsq31 * c13sq * Dmsq21 * c12sq;
}
// See eq. 45
double E(double D)
{
	return (Dmsq31 * (M3sq(D) - Dmsq21) - Dmsq21 * (M3sq(D) - Dmsq31) * s12sq) * c13 * s13;
}
// See eq. 46
double F(double D)
{
	return (M3sq(D) - Dmsq31) * Dmsq21 * c12 * s12 * c13;
}
// angles in matter
// See eq. 47
double s12msq(double D)
{
	double x = (pow(M2sq(D), 2) - alpha() * M2sq(D) + beta()) * DM31sq(D);
	double s12msq_n = -x;
	double s12msq_d = DM32sq(D) * (pow(M1sq(D), 2) - alpha() * M1sq(D) + beta()) - x;
	return s12msq_n / s12msq_d;
}
double c12msq(double D)
{
	return 1 - s12msq(D);
}
// See eq. 48
double s13msq(double D)
{
	return (pow(M3sq(D), 2) - alpha() * M3sq(D) + beta()) / (DM31sq(D) * DM32sq(D));
}
double c13msq(double D)
{
	return 1 - s13msq(D);
}
// See eq. 49
double s23msq(double D, double delta)
{
	return (pow(E(D), 2) * s23sq + pow(F(D), 2) * c23sq + 2 * E(D) * F(D) * c23 * s23 * cos(delta)) / (pow(E(D), 2) + pow(F(D), 2));
}
double c23msq(double D, double delta)
{
	return 1 - s23msq(D, delta);
}
double cdeltam(double D, double delta)
{
	// tmp
	double Et = E(D);
	double Ft = F(D);
	// denominator
	double densq1 = pow(Et, 2) * s23sq + pow(Ft, 2) * c23sq + 2 * Et * Ft * c23 * s23 * cos(delta);
	double densq2 = pow(Et, 2) * c23sq + pow(Ft, 2) * s23sq - 2 * Et * Ft * c23 * s23 * cos(delta);
	// numerator c23sq -> c23 typo corrected
	return ((pow(Et, 2) - pow(Ft, 2)) * cos(delta) * s23 * c23 + Et * Ft * (c23sq - s23sq)) / sqrt(densq1 * densq2);
}
double sdeltam(double D, double delta)
{
	// tmp
	double Et = E(D);
	double Ft = F(D);
	// denominator
	double densq1 = pow(Et, 2) * s23sq + pow(Ft, 2) * c23sq + 2 * Et * Ft * c23 * s23 * cos(delta);
	double densq2 = pow(Et, 2) * c23sq + pow(Ft, 2) * s23sq - 2 * Et * Ft * c23 * s23 * cos(delta);
	// numerator c23sq -> c23 typo corrected
	return ((pow(Et, 2) + pow(Ft, 2)) * sin(delta) * s23 * c23) / sqrt(densq1 * densq2);
}
Matrix<std::complex<double> > UMNSm(double a, double delta)
{
	double t12m = atan2(sqrt(s12msq(a)), sqrt(c12msq(a)));
	double t13m = atan2(sqrt(s13msq(a)), sqrt(c13msq(a)));
	double t23m = atan2(sqrt(s23msq(a, delta)), sqrt(c23msq(a, delta)));
	double deltam = atan2(sdeltam(a, delta), cdeltam(a, delta));

	Matrix<std::complex<double> > U12m = U12(t12m);
	Matrix<std::complex<double> > U13m = U13(t13m);
	Matrix<std::complex<double> > U23m = U23(deltam, t23m);

	return U23m * U13m * U12m;
}
std::complex<double> Talphabetai(int alpha, int beta, int i, double a, double delta)
{
	assert (alpha >= 0 && alpha <= 2); // alpha valid for 0, 1, 2
	assert (beta >= 0 && beta <= 2); // beta valid for 0, 1, 2
	assert (i >= 0 && i <= 2); // i valid for 0, 1, 2

	Matrix<std::complex<double> > U = UMNSm(a, delta);
	return std::conj(U(alpha, i)) * U(beta, i);
}
double Palphabeta(int alpha, int beta, double a, double LE, double delta)
{
	assert (alpha >= 0 && alpha <= 2); // alpha valid for 0, 1, 2
	assert (beta >= 0 && beta <= 2); // beta valid for 0, 1, 2

	double L2E = km_per_GeV_to_per_eV2 * LE / 2; // in eV^-2 now

	// See eq. 4.0.1 or eq. 3.12 in arXiv:1505.01826
	std::complex<double> A(0, 0);
	A += Talphabetai(alpha, beta, 0, a, delta) * std::complex<double>(cos(M1sq(a) * L2E), -sin(M1sq(a) * L2E));
	A += Talphabetai(alpha, beta, 1, a, delta) * std::complex<double>(cos(M2sq(a) * L2E), -sin(M2sq(a) * L2E));
	A += Talphabetai(alpha, beta, 2, a, delta) * std::complex<double>(cos(M3sq(a) * L2E), -sin(M3sq(a) * L2E));

	return std::norm(A);
}
} // namespace Exact
