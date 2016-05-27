#include <cmath>
#include <complex>
#include <cassert>

#include "Check.h"
#include "Constants.h"
#include "Hat.h"
#include "Matrix.h"
#include "Warning.h"

namespace Check
{
double l1(double a, int order)
{
	assert (order >= 0 && order <= 2); // order valid for 0, 1, 2
	double l1 = l10(a);
	if (order == 2)
		l1 += l12(a);
	return l1;
}
double l2(double a, int order)
{
	assert (order >= 0 && order <= 2); // order valid for 0, 1, 2
	double l2 = l20(a);
	if (order == 2)
		l2 += l22(a);
	return l2;
}
double l3(double a, int order)
{
	assert (order >= 0 && order <= 2); // order valid for 0, 1, 2
	double l3 = l30(a);
	if (order == 2)
		l3 += l32(a);
	return l3;
}
// See eq. 2.4.5
double l10(double a)
{
	double c = cos(Hat::phi(a) - t13);
	double x = Hat::l0(a) + Hat::lm(a);
	double y = pow(Hat::l0(a) - Hat::lm(a), 2) + pow(2 * eps * c12 * s12 * Dmsqee * c, 2);
	return (x - sqrt(y)) / 2.;
}
// See eq. 2.4.5
double l20(double a)
{
	double c = cos(Hat::phi(a) - t13);
	double x = Hat::l0(a) + Hat::lm(a);
	double y = pow(Hat::l0(a) - Hat::lm(a), 2) + pow(2 * eps * c12 * s12 * Dmsqee * c, 2);
	return (x + sqrt(y)) / 2.;
}
// See eq. 2.4.5
double l30(double a)
{
	return Hat::lp(a);
}
// See eq. 3.1.3
double l12(double a)
{
	return -pow(epsp(a) * Dmsqee * sin(psi(a)), 2) / (l30(a) - l10(a));
}
// See eq. 3.1.3
double l22(double a)
{
	return -pow(epsp(a) * Dmsqee * cos(psi(a)), 2) / (l30(a) - l20(a));
}
// See eq. 3.1.3
double l32(double a)
{
	return pow(epsp(a) * Dmsqee, 2) * (pow(sin(psi(a)), 2) / (l30(a) - l10(a)) + pow(cos(psi(a)), 2) / (l30(a) - l20(a)));
}
double Dl21(double a, int order)
{
	assert (order >= 0 && order <= 2); // order valid for 0, 1, 2
	if (order < 2)
		return l20(a) - l10(a);
	else
		return l2(a) - l1(a);
}
double Dl31(double a, int order)
{
	assert (order >= 0 && order <= 2); // order valid for 0, 1, 2
	if (order < 2)
		return l30(a) - l10(a);
	else
		return l3(a) - l1(a);
}
double Dl32(double a, int order)
{
	assert (order >= 0 && order <= 2); // order valid for 0, 1, 2
	if (order < 2)
		return l30(a) - l20(a);
	else
		return l3(a) - l2(a);
}
// See eq. 2.4.9
double c2psi(double a)
{
	return (Hat::l0(a) - Hat::lm(a)) / (l20(a) - l10(a));
}
// See eq. 2.4.7
double s2psi(double a)
{
	double c = cos(Hat::phi(a) - t13);
	return 2 * eps * Dmsqee * c12 * s12 * c / (l20(a) - l10(a));
}
double psi(double a)
{
	return atan2(s2psi(a), c2psi(a)) / 2.;
}
// See table 1
double Jrm(double a)
{
	return c23 * s23 * pow(cos(Hat::phi(a)), 2) * sin(Hat::phi(a)) * cos(psi(a)) * sin(psi(a));
}
Matrix<std::complex<double> > W(double a, int order)
{
	assert (order >= 0 && order <= 2); // order valid for 0, 1, 2

	Matrix<std::complex<double> > W(3, 3);
	double e31, e32, coef;

	switch (order)
	{
		case 2:
			// Second order from eq. 3.2.6
			coef = -pow(epsp(a) * Dmsqee, 2) / 2.;
			W(0, 0) += coef * pow(sin(psi(a)) / Dl31(a), 2);
			W(1, 1) += coef * pow(cos(psi(a)) / Dl32(a), 2);
			W(2, 2) += coef * (pow(sin(psi(a)) / Dl31(a), 2) + pow(cos(psi(a)) / Dl32(a), 2));
			W(1, 0) += coef * s2psi(a) / (Dl31(a) * Dl21(a));
			W(0, 1) += -coef * s2psi(a) / (Dl32(a) * Dl21(a));

		case 1:
			// First order from eq. 3.2.5
			e31 = epsp(a) * Dmsqee * sin(psi(a)) / Dl31(a);
			e32 = epsp(a) * Dmsqee * cos(psi(a)) / Dl32(a);
			W(0, 2) += -e31;
			W(1, 2) += e32;
			W(2, 0) += e31;
			W(2, 1) += -e32;

		case 0:
			// zeroth order
			for (int i = 0; i < 3; i++)
				W(i, i) += 1;
	}
	return W;
}
// See eq. 2.5.2
Matrix<std::complex<double> > UMNSm(double a, double delta)
{
	return U23(delta) * U13(Hat::phi(a)) * U12(psi(a));
}
// See eq. 3.2.3
Matrix<std::complex<double> > V(double a, double delta, int order)
{
	assert (order >= 0 && order <= 2); // order valid for 0, 1, 2
	Matrix<std::complex<double> > W = Check::W(a, order);

	return UMNSm(a, delta) * W;
}
std::complex<double> Talphabetai(int alpha, int beta, int i, double a, double delta, int order)
{
	assert (alpha >= 0 && alpha <= 2); // alpha valid for 0, 1, 2
	assert (beta >= 0 && beta <= 2); // beta valid for 0, 1, 2
	assert (i >= 0 && i <= 2); // i valid for 0, 1, 2
	assert (order >= 0 && order <= 2); // order valid for 0, 1, 2

	Matrix<std::complex<double> > Vtmp = V(a, delta, order);
	return std::conj(Vtmp(alpha, i)) * Vtmp(beta, i);
}
double Palphabeta(int alpha, int beta, double a, double LE, double delta, int order)
{
	assert (alpha >= 0 && alpha <= 2); // alpha valid for 0, 1, 2
	assert (beta >= 0 && beta <= 2); // beta valid for 0, 1, 2
	assert (order >= 0 && order <= 2); // order valid for 0, 1, 2

	double L2E = km_per_GeV_to_per_eV2 * LE / 2; // in eV^-2 now

	double l1, l2, l3;
	if (order < 2) // first order is same as zeroth order
	{
		l1 = l10(a);
		l2 = l20(a);
		l3 = l30(a);
	}
	else
	{
		l1 = Check::l1(a);
		l2 = Check::l2(a);
		l3 = Check::l3(a);
	}

	// See eq. 4.0.1 or eq. 3.12 in arXiv:1505.01826
	std::complex<double> A(0, 0);
	A += Talphabetai(alpha, beta, 0, a, delta, order) * std::complex<double>(cos(l1 * L2E), -sin(l1 * L2E));
	A += Talphabetai(alpha, beta, 1, a, delta, order) * std::complex<double>(cos(l2 * L2E), -sin(l2 * L2E));
	A += Talphabetai(alpha, beta, 2, a, delta, order) * std::complex<double>(cos(l3 * L2E), -sin(l3 * L2E));

	double P = std::norm(A);
	ProbabilityWarning W(&P);
	return P;
}
} // namespace Check
