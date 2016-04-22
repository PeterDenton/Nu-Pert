#ifndef Hat_H
#define Hat_H

namespace Hat
{
double c2phi(double a);
double s2phi(double a);
double phi(double a);

double lm(double a);
double l0(double a);
double lp(double a);

// LE = L/E in km/GeV
double Pee(double a, double LE);
double Pemu(double a, double LE, double delta);
} // namespace Hat

#endif
