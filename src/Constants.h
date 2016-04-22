#ifndef Constants_H
#define Constants_H

#include <complex>

#include "Matrix.h"

// from the pdg
extern const double s12sq, s23sq, s13sq;
extern const double Dmsq21, Dmsq31_base;
extern double Dmsq31;

// calculated from the above
extern const double c12sq, c23sq, c13sq;
extern const double c12, s12, c23, s23, c13, s13;
extern const double t12, t23, t13;

extern const double Jr; // reduced Jarlskog

// for the perturbative calculation
extern double Dmsqee;
extern double eps;
double epsp(double a);

// units
extern const double km_per_GeV_to_per_eV2; // x (km/GeV) * km_per_GeV_to_per_eV2 = y eV^-2
extern const double Y_to_a; // x Y_e*rho*E (g*GeV/cc) * Y_to_a = y eV^2

extern bool normal;
extern int normal_sign;
void set_ordering(bool _normal); // default is set to true (normal)

Matrix<std::complex<double> > U23(double delta, double angle = t23);
Matrix<double> U13(double angle = t13);
Matrix<double> U12(double angle = t12);

Matrix<std::complex<double> > UPMNS(double delta);

#endif
