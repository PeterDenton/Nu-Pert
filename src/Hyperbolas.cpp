/* Use the exact solution but with hyperbola inspired eigenvalues */

#include <cmath>
#include <complex>
#include <cassert>

#include "Mixed.h"
#include "Constants.h"
#include "Matrix.h"
#include "Hyperbolas.h"

namespace Hyperbolas
{
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

/* the eigenvalues */
double la(double a)
{
	return a + (s13sq + eps * s12sq) * Dmsqee;
}
double lb = eps * c12sq * Dmsqee;
double lc = (c13sq + eps * s12sq) * Dmsqee;

double X12(double a)
{
	return sin(2 * t12) * eps * Dmsqee;
}
double X23 = sin(2 * t13) * Dmsqee;

double l1(double a)
{
	double ax = cos(2 * t12) * eps * Dmsqee;
	return 0.5 * ((a - ax) - sqrt(pow(a - ax, 2) + pow(X12(a), 2))) + lb;
}
double l2(double a)
{
	double ax = cos(2 * t12) * eps * Dmsqee;
	double ax12 = lb - la(0);
	double ax13 = lc - la(0);
	double atr = (ax13 - ax12) / 2.; // ansatz that simple average is right, perhaps weighted average by the X's is better

	if (a < atr)
		return 0.5 * ((a - ax) + sqrt(pow(a - ax, 2) + pow(X12(a), 2))) + lb;
	else
		return 0.5 * ((la(a) + lc) - sqrt(pow(la(a) - lc, 2) + pow(X23, 2)));

}
double l3(double a)
{
	return 0.5 * ((la(a) + lc) + sqrt(pow(la(a) - lc, 2) + pow(X23, 2)));
}
double Dl21(double a)
{
	return l2(a) - l1(a);
}
double Dl31(double a)
{
	return l3(a) - l1(a);
}
double Dl32(double a)
{
	return l3(a) - l2(a);
}
/* end eigenvalues */

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
	return (Dmsq31 * (Hyperbolas::l3(D) - Dmsq21) - Dmsq21 * (Hyperbolas::l3(D) - Dmsq31) * s12sq) * c13 * s13;
}
// See eq. 46
double F(double D)
{
	return (Hyperbolas::l3(D) - Dmsq31) * Dmsq21 * c12 * s12 * c13;
}
// angles in matter
// See eq. 47
double s12msq(double D)
{
	double x = (pow(Hyperbolas::l2(D), 2) - alpha() * Hyperbolas::l2(D) + beta()) * Hyperbolas::Dl31(D);
	double s12msq_n = -x;
	double s12msq_d = Hyperbolas::Dl32(D) * (pow(Hyperbolas::l1(D), 2) - alpha() * Hyperbolas::l1(D) + beta()) - x;
	return s12msq_n / s12msq_d;
}
double c12msq(double D)
{
	return 1 - s12msq(D);
}
// See eq. 48
double s13msq(double D)
{
	return (pow(Hyperbolas::l3(D), 2) - alpha() * Hyperbolas::l3(D) + beta()) / (Hyperbolas::Dl31(D) * Hyperbolas::Dl32(D));
}
double c13msq(double D)
{
	return 1 - s13msq(D);
}
// See eq. 49
double s23msq(double D, double delta)
{
	return (pow(E(D), 2) * s23sq + pow(F(D), 2) * c23sq + 2 * E(D) * F(D) * c23 * s23 * cos(delta)) /
			(pow(E(D), 2) + pow(F(D), 2));
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
	A += Talphabetai(alpha, beta, 0, a, delta) * std::complex<double>(cos(Hyperbolas::l1(a) * L2E), -sin(Hyperbolas::l1(a) * L2E));
	A += Talphabetai(alpha, beta, 1, a, delta) * std::complex<double>(cos(Hyperbolas::l2(a) * L2E), -sin(Hyperbolas::l2(a) * L2E));
	A += Talphabetai(alpha, beta, 2, a, delta) * std::complex<double>(cos(Hyperbolas::l3(a) * L2E), -sin(Hyperbolas::l3(a) * L2E));

	return std::norm(A);
}
} // namespace Hyperbolas
