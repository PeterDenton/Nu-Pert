#include <cassert>

#include "Tilde.h"
#include "Constants.h"
#include "Warning.h"

namespace Tilde
{
// See eq. 2.2.5
double la(double a)
{
	return a + (s13sq + eps * s12sq) * Dmsqee;
}
// See eq. 2.2.5
double lb(double a)
{
	return eps * c12sq * Dmsqee;
}
// See eq. 2.2.5
double lc(double a)
{
	return (c13sq + eps * s12sq) * Dmsqee;
}
double Dlba(double a)
{
	return lb(a) - la(a);
}
double Dlca(double a)
{
	return lc(a) - la(a);
}
double Dlcb(double a)
{
	return lc(a) - lb(a);
}
Matrix<std::complex<double> > W(double a, int order)
{
	assert (order >= 0 && order <= 1); // order valid for 0, 1

	Matrix<std::complex<double> > W(3, 3);
	double coef;

	switch (order)
	{
		case 1:
			W(0, 2) += c13 * s13 * Dmsqee / Dlca(a);
			W(2, 0) -= c13 * s13 * Dmsqee / Dlca(a);

			coef = eps * c12 * s12 * Dmsqee;
			W(0, 1) += coef * c13 / Dlba(a);
			W(1, 0) -= coef * c13 / Dlba(a);
			W(1, 2) -= coef * s13 / Dlcb(a);
			W(2, 1) += coef * s13 / Dlcb(a);

		case 0:
			// zeroth order
			for (int i = 0; i < 3; i++)
				W(i, i) += 1;
	}
	return W;
}
Matrix<std::complex<double> > UMNSm(double a, double delta)
{
	return U23(delta);
}
Matrix<std::complex<double> > V(double a, double delta, int order)
{
	assert (order >= 0 && order <= 1); // order valid for 0, 1
	Matrix<std::complex<double> > W = Tilde::W(a, order);

	return UMNSm(a, delta) * W;
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

	double la = Tilde::la(a);
	double lb = Tilde::lb(a);
	double lc = Tilde::lc(a);

	// See eq. 4.0.1 or eq. 3.12 in arXiv:1505.01826
	std::complex<double> A(0, 0);
	A += Talphabetai(alpha, beta, 0, a, delta, order) * std::complex<double>(cos(la * L2E), -sin(la * L2E));
	A += Talphabetai(alpha, beta, 1, a, delta, order) * std::complex<double>(cos(lb * L2E), -sin(lb * L2E));
	A += Talphabetai(alpha, beta, 2, a, delta, order) * std::complex<double>(cos(lc * L2E), -sin(lc * L2E));

	double P = std::norm(A);
	ProbabilityWarning W(&P);
	return P;
}
} // namespace Tilde
