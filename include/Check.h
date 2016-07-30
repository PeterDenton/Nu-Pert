#ifndef Check_H
#define Check_H

#include <complex>

#include "Matrix.h"

namespace Check
{
double l1(double a, int order = 2);
double l2(double a, int order = 2);
double l3(double a, int order = 2);

double l10(double a);
double l20(double a);
double l30(double a);

double l12(double a);
double l22(double a);
double l32(double a);

double Dl21(double a, int order = 0);
double Dl31(double a, int order = 0);
double Dl32(double a, int order = 0);

double c2psi(double a);
double s2psi(double a);
double psi(double a);

double Jrm(double a);
double Jrrm(double a);

Matrix<std::complex<double> > W(double a, int order = 1);
Matrix<std::complex<double> > V(double a, double delta, int order = 1);
Matrix<std::complex<double> > UMNSm(double a, double delta);

double Palphabeta(int alpha, int beta, double a, double LE, double delta, int order);
} // namespace Check

#endif
