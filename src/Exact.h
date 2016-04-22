#ifndef Exact_H
#define Exact_H

#include <complex>

#include "Matrix.h"


namespace Exact
{
double M1sq(double D);
double M2sq(double D);
double M3sq(double D);

double DM21sq(double D);
double DM31sq(double D);
double DM32sq(double D);

Matrix<std::complex<double> > UMNSm(double a, double delta);
std::complex<double> Talphabetai(int alpha, int beta, int i, double a, double delta);

// alpha, beta = 0, 1, 2 for e, mu, tau
double Palphabeta(int alpha, int beta, double a, double LE, double delta);
} // namespace Exact

#endif
