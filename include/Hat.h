#ifndef Hat_H
#define Hat_H

#include <complex>

#include "Matrix.h"

namespace Hat
{
double c2phi(double a);
double s2phi(double a);
double phi(double a);

double lm(double a);
double l0(double a);
double lp(double a);

double Dl0m(double a);
double Dlpm(double a);
double Dlp0(double a);

// LE = L/E in km/GeV
double Pee(double a, double LE);
double Pemu(double a, double LE, double delta);

Matrix<double> W(double a, int order = 1);
Matrix<std::complex<double> > V(double a, double delta, int order = 1);
Matrix<std::complex<double> > UMNSm(double a, double delta);

double Palphabeta(int alpha, int beta, double a, double LE, double delta, int order);
} // namespace Hat

#endif
