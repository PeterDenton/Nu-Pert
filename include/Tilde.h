#ifndef Tilde_H
#define Tilde_H

#include <complex>

#include "Matrix.h"

namespace Tilde
{
double la(double a);
double lb(double a);
double lc(double a);

double Dlba(double a);
double Dlca(double a);
double Dlcb(double a);

Matrix<std::complex<double> > W(double a, int order = 1);
Matrix<std::complex<double> > V(double a, double delta, int order = 1);
Matrix<std::complex<double> > UMNSm(double a, double delta);

double Palphabeta(int alpha, int beta, double a, double LE, double delta, int order);
} // namespace Tilde

#endif
