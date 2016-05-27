#ifndef Mixed_H
#define Mixed_H

#include <complex>

#include "Matrix.h"

namespace Mixed
{
Matrix<std::complex<double> > UMNSm(double a, double delta, int order);
std::complex<double> Talphabetai(int alpha, int beta, int i, double a, double delta, int order);

// alpha, beta = 0, 1, 2 for e, mu, tau
// order = 0, (1), 2 is the order for the eigenvalues ONLY
double Palphabeta(int alpha, int beta, double a, double LE, double delta, int order);
} // namespace Mixed

#endif
