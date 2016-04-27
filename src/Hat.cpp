/*
This basis is that after two rotations and is the same as that presented in arXiv:1505.01826.
Equation numbers in parantheses refer to arXiv:1505.01826.
*/
#include <cmath>
#include <complex>

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
}; // namespace Hat

