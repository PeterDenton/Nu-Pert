/* Use the exact solution but with check eigenvalues */

#include <cmath>
#include <complex>
#include <cassert>

#include "Mixed.h"
#include "Constants.h"
#include "Matrix.h"
#include "Check.h"

namespace Mixed
{
double alpha();
double beta();

double E(double D, int order);
double F(double D, int order);

double s12msq(double D, int order);
double c12msq(double D, int order);
double s13msq(double D, int order);
double c13msq(double D, int order);
double s23msq(double D, double delta, int order);
double c23msq(double D, double delta, int order);
double cdeltam(double D, double delta, int order);
double sdeltam(double D, double delta, int order);

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
double E(double D, int order)
{
	return (Dmsq31 * (Check::l3(D, order) - Dmsq21) - Dmsq21 * (Check::l3(D, order) - Dmsq31) * s12sq) * c13 * s13;
}
// See eq. 46
double F(double D, int order)
{
	return (Check::l3(D, order) - Dmsq31) * Dmsq21 * c12 * s12 * c13;
}
// angles in matter
// See eq. 47
double s12msq(double D, int order)
{
	double x = (pow(Check::l2(D, order), 2) - alpha() * Check::l2(D, order) + beta()) * Check::Dl31(D, order);
	double s12msq_n = -x;
	double s12msq_d = Check::Dl32(D, order) * (pow(Check::l1(D, order), 2) - alpha() * Check::l1(D, order) + beta()) - x;
	return s12msq_n / s12msq_d;
}
double c12msq(double D, int order)
{
	return 1 - s12msq(D, order);
}
// See eq. 48
double s13msq(double D, int order)
{
	return (pow(Check::l3(D, order), 2) - alpha() * Check::l3(D, order) + beta()) / (Check::Dl31(D, order) * Check::Dl32(D, order));
}
double c13msq(double D, int order)
{
	return 1 - s13msq(D, order);
}
// See eq. 49
double s23msq(double D, double delta, int order)
{
	return (pow(E(D, order), 2) * s23sq + pow(F(D, order), 2) * c23sq + 2 * E(D, order) * F(D, order) * c23 * s23 * cos(delta)) /
			(pow(E(D, order), 2) + pow(F(D, order), 2));
}
double c23msq(double D, double delta, int order)
{
	return 1 - s23msq(D, delta, order);
}
double cdeltam(double D, double delta, int order)
{
	// tmp
	double Et = E(D, order);
	double Ft = F(D, order);
	// denominator
	double densq1 = pow(Et, 2) * s23sq + pow(Ft, 2) * c23sq + 2 * Et * Ft * c23 * s23 * cos(delta);
	double densq2 = pow(Et, 2) * c23sq + pow(Ft, 2) * s23sq - 2 * Et * Ft * c23 * s23 * cos(delta);
	// numerator c23sq -> c23 typo corrected
	return ((pow(Et, 2) - pow(Ft, 2)) * cos(delta) * s23 * c23 + Et * Ft * (c23sq - s23sq)) / sqrt(densq1 * densq2);
}
double sdeltam(double D, double delta, int order)
{
	// tmp
	double Et = E(D, order);
	double Ft = F(D, order);
	// denominator
	double densq1 = pow(Et, 2) * s23sq + pow(Ft, 2) * c23sq + 2 * Et * Ft * c23 * s23 * cos(delta);
	double densq2 = pow(Et, 2) * c23sq + pow(Ft, 2) * s23sq - 2 * Et * Ft * c23 * s23 * cos(delta);
	// numerator c23sq -> c23 typo corrected
	return ((pow(Et, 2) + pow(Ft, 2)) * sin(delta) * s23 * c23) / sqrt(densq1 * densq2);
}
Matrix<std::complex<double> > UMNSm(double a, double delta, int order)
{
	double t12m = atan2(sqrt(s12msq(a, order)), sqrt(c12msq(a, order)));
	double t13m = atan2(sqrt(s13msq(a, order)), sqrt(c13msq(a, order)));
	double t23m = atan2(sqrt(s23msq(a, delta, order)), sqrt(c23msq(a, delta, order)));
	double deltam = atan2(sdeltam(a, delta, order), cdeltam(a, delta, order));

	Matrix<std::complex<double> > U12m = U12(t12m);
	Matrix<std::complex<double> > U13m = U13(t13m);
	Matrix<std::complex<double> > U23m = U23(deltam, t23m);

	return U23m * U13m * U12m;
}
std::complex<double> Talphabetai(int alpha, int beta, int i, double a, double delta, int order)
{
	assert (alpha >= 0 && alpha <= 2); // alpha valid for 0, 1, 2
	assert (beta >= 0 && beta <= 2); // beta valid for 0, 1, 2
	assert (i >= 0 && i <= 2); // i valid for 0, 1, 2

	Matrix<std::complex<double> > U = UMNSm(a, delta, order);
	return std::conj(U(alpha, i)) * U(beta, i);
}
double Palphabeta(int alpha, int beta, double a, double LE, double delta, int order)
{
	assert (alpha >= 0 && alpha <= 2); // alpha valid for 0, 1, 2
	assert (beta >= 0 && beta <= 2); // beta valid for 0, 1, 2
	assert (order >= 0 && order <= 2); // order valid for 0, (1), 2

	double L2E = km_per_GeV_to_per_eV2 * LE / 2; // in eV^-2 now

	// See eq. 4.0.1 or eq. 3.12 in arXiv:1505.01826
	std::complex<double> A(0, 0);
	A += Talphabetai(alpha, beta, 0, a, delta, order) * std::complex<double>(cos(Check::l1(a, order) * L2E), -sin(Check::l1(a, order) * L2E));
	A += Talphabetai(alpha, beta, 1, a, delta, order) * std::complex<double>(cos(Check::l2(a, order) * L2E), -sin(Check::l2(a, order) * L2E));
	A += Talphabetai(alpha, beta, 2, a, delta, order) * std::complex<double>(cos(Check::l3(a, order) * L2E), -sin(Check::l3(a, order) * L2E));

	return std::norm(A);
}
} // namespace Mixed
