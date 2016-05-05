#include <cmath>
#include <complex>

#include "Constants.h"
#include "Hat.h"

const double s12sq = 0.308;
const double s23sq = 0.446;
const double s13sq = 0.0237;
const double Dmsq21 = 7.54e-5; // eV^2
const double Dmsq31_base = 2.43e-3; // eV^2
double Dmsq31 = Dmsq31_base; // default is normal hierarchy

const double c12sq = 1 - s12sq;
const double c23sq = 1 - s23sq;
const double c13sq = 1 - s13sq;

const double s12 = sqrt(s12sq);
const double s23 = sqrt(s23sq);
const double s13 = sqrt(s13sq);

const double c12 = sqrt(c12sq);
const double c23 = sqrt(c23sq);
const double c13 = sqrt(c13sq);

const double t12 = asin(s12);
const double t23 = asin(s23);
const double t13 = asin(s13);

// See eq. A.7.2
const double Jr = c12 * s12 * c23 * s23 * c13sq * s13;

// See eq. 2.2.3
double Dmsqee = Dmsq31 - s12sq * Dmsq21;
// See eq. 2.2.4
double eps = Dmsq21 / Dmsqee;
// See eq. 2.5.4
double epsp(double a)
{
	return eps * sin(Hat::phi(a) - t13) * c12 * s12;
}

const double km_per_GeV_to_per_eV2 = 1e-6 / (1.97327e-7);
const double Y_to_a = 1.52e-4;

bool normal = true;
int normal_sign = 1;
void set_ordering(bool _normal)
{
	normal = _normal;
	normal_sign = normal ? 1 : -1;

	Dmsq31 = normal_sign * Dmsq31_base;

	Dmsqee = Dmsq31 - s12sq * Dmsq21;
	eps = Dmsq21 / Dmsqee;
}
// See eq. 2.1.4
Matrix<std::complex<double> > U23(double delta, double angle)
{
	std::complex<double> eid(cos(delta), sin(delta));
	Matrix<std::complex<double> > ret(3, 3);
	ret(0, 0) = 1;
	ret(1, 1) = cos(angle);
	ret(1, 2) = sin(angle) * eid;
	ret(2, 1) = -sin(angle) * std::conj(eid);
	ret(2, 2) = cos(angle);
	return ret;
}
// See eq. 2.1.4
Matrix<std::complex<double> > U13(double angle)
{
	Matrix<std::complex<double> > ret(3, 3);
	ret(0, 0) = cos(angle);
	ret(0, 2) = sin(angle);
	ret(1, 1) = 1;
	ret(2, 0) = -sin(angle);
	ret(2, 2) = cos(angle);
	return ret;
}
// See eq. 2.1.4
Matrix<std::complex<double> > U12(double angle)
{
	Matrix<std::complex<double> > ret(3, 3);
	ret(0, 0) = cos(angle);
	ret(0, 1) = sin(angle);
	ret(1, 0) = -sin(angle);
	ret(1, 1) = cos(angle);
	ret(2, 2) = 1;
	return ret;
}
// See before eq. 2.1.4
Matrix<std::complex<double> > UMNS(double delta)
{
	return U23(delta) * U13() * U12();
}
