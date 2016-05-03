/*
This basis is that after two rotations and is the same as that presented in arXiv:1505.01826.
Equation numbers in parantheses refer to arXiv:1505.01826.
*/
#include <cmath>
#include <complex>
#include <cassert>

#include "Hat.h"
#include "Constants.h"
#include "Warning.h"

namespace Hat
{
// See eq. 2.3.5 (2.5)
double c2phi(double a)
{
	return (Dmsqee * cos(2 * t13) - a) / (lp(a) - lm(a));
}
// See eq. 2.3.4 (2.5)
double s2phi(double a)
{
	return Dmsqee * sin(2 * t13) / (lp(a) - lm(a));
}
double phi(double a)
{
	return atan2(s2phi(a), c2phi(a)) / 2.;
}
// See eq. 2.3.3 (3.1)
double lm(double a)
{
	double x = Dmsqee + a;
	double y = normal_sign * sqrt(pow(Dmsqee - a, 2) + 4 * s13sq * a * Dmsqee);
	double z = eps * Dmsqee * s12sq;
	return (x - y) / 2. + z;
}
// See eq. 2.3.3 (3.1)
double l0(double a)
{
	return c12sq * eps * Dmsqee;
}
// See eq. 2.3.3 (3.1)
double lp(double a)
{
	double x = Dmsqee + a;
	double y = normal_sign * sqrt(pow(Dmsqee - a, 2) + 4 * s13sq * a * Dmsqee);
	double z = eps * Dmsqee * s12sq;
	return (x + y) / 2. + z;
}
// From eq. (1.1)
double Dl0m(double a)
{
	return l0(a) - lm(a);
}
double Dlpm(double a)
{
	return lp(a) - lm(a);
}
double Dlp0(double a)
{
	return lp(a) - l0(a);
}
double Pee(double a, double LE)
{
	double L4E = km_per_GeV_to_per_eV2 * LE / 4; // in eV^-2 now
	double Deltapm = (lp(a) - lm(a)) * L4E;
	double P = 1 - pow(s2phi(a) * sin(Deltapm), 2);
	ProbabilityWarning W(&P);
	return P;
}
// From eq. (2.14)
double Pemu(double a, double LE, double delta)
{
	double L4E = km_per_GeV_to_per_eV2 * LE / 4; // in eV^-2 now

	double Dlpm = lp(a) - lm(a);
	double Dlp0 = lp(a) - l0(a);
	double Dlm0 = lm(a) - l0(a);

	double Deltapm = Dlpm * L4E;
	double Deltap0 = Dlp0 * L4E;
	double Deltam0 = Dlm0 * L4E;

	double A = s23sq * pow(sin(2 * t13), 2) + 4 * eps * Jr * cos(delta) * (Dlpm - (Dmsqee - a)) / Dlp0;
	double B = 8 * eps * Jr * pow(Dmsqee, 3) / (Dlpm * Dlp0 * Dlm0);
	double C = sin(Deltapm) * sin(Deltam0) * cos(delta - Deltap0);

	double P = A * pow(Dmsqee * sin(Deltapm) / Dlpm, 2) + B * C;
	ProbabilityWarning W(&P);
	return P;
}
Matrix<double> W(double a, int order)
{
	assert (order >= 0 && order <= 1); // order valid for 0, 1

	Matrix<double> W(3, 3);
	double c = cos(phi(a) - t13);
	double s = sin(phi(a) - t13);
	double coef = eps * c12 * s12 * Dmsqee;

	switch (order)
	{
		case 1:
			W(0, 1) += coef * c / Dl0m(a);
			W(1, 0) -= coef * c / Dl0m(a);
			W(1, 2) += coef * s / Dlp0(a);
			W(2, 1) -= coef * s / Dlp0(a);

		case 0:
			// zeroth order
			for (int i = 0; i < 3; i++)
				W(i, i) += 1;
	}
	return W;
}
Matrix<std::complex<double> > UMNSm(double a, double delta)
{
	return U23(delta) * double2complex(U13(phi(a)));
}
Matrix<std::complex<double> > V(double a, double delta, int order)
{
	assert (order >= 0 && order <= 1); // order valid for 0, 1
	Matrix<double> W = Hat::W(a, order);

	return UMNSm(a, delta) * double2complex(W);
}
std::complex<double> Talphabetai(int alpha, int beta, int i, double a, double delta, int order)
{
	assert (alpha >= 0 && alpha <= 2); // alpha valid for 0, 1, 2
	assert (beta >= 0 && beta <= 2); // beta valid for 0, 1, 2
	assert (i >= 0 && i <= 2); // i valid for 0, 1, 2
	assert (order >= 0 && order <= 1); // order valid for 0, 1

	Matrix<std::complex<double> > Vtmp = V(a, delta, order);
	return std::conj(Vtmp(alpha, i)) * Vtmp(beta, i);
}
double Palphabeta(int alpha, int beta, double a, double LE, double delta, int order)
{
	assert (alpha >= 0 && alpha <= 2); // alpha valid for 0, 1, 2
	assert (beta >= 0 && beta <= 2); // beta valid for 0, 1, 2
	assert (order >= 0 && order <= 1); // order valid for 0, 1

	double L2E = km_per_GeV_to_per_eV2 * LE / 2; // in eV^-2 now

	double lm = Hat::lm(a);
	double l0 = Hat::l0(a);
	double lp = Hat::lp(a);

	std::complex<double> A(0, 0);
	A += Talphabetai(alpha, beta, 0, a, delta, order) * std::complex<double>(cos(lm * L2E), -sin(lm * L2E));
	A += Talphabetai(alpha, beta, 1, a, delta, order) * std::complex<double>(cos(l0 * L2E), -sin(l0 * L2E));
	A += Talphabetai(alpha, beta, 2, a, delta, order) * std::complex<double>(cos(lp * L2E), -sin(lp * L2E));

	double P = std::norm(A);
	ProbabilityWarning W(&P);
	return P;
}
} // namespace Hat

