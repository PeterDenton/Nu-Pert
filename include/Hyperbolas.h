#ifndef Hyperbolas_H
#define Hyperbolas_H

#include <complex>

#include "Matrix.h"

namespace Hyperbolas
{
double l1(double a);
double l2(double a);
double l3(double a);

Matrix<std::complex<double> > UMNSm(double a, double delta);
std::complex<double> Talphabetai(int alpha, int beta, int i, double a, double delta);

// alpha, beta = 0, 1, 2 for e, mu, tau
// order = 0, (1), 2 is the order for the eigenvalues ONLY
double Palphabeta(int alpha, int beta, double a, double LE, double delta);
} // namespace Hyperbolas

#endif
