/* From the Li, Wang, Xing paper arXiv:1605.00900 using formula 18.
Equation numbers throughout are in reference to the LWX paper */

#include "LWX.h"
#include "Constants.h"

namespace LWX
{
// eq 13
double P0(double LE)
{
	return pow(sin(2 * t12), 2) * pow(c13, 4) * pow(sin(F21(LE)), 2);
}
double Pstar(double LE)
{
	return 0.5 * pow(sin(2 * t13), 2) * (1 - cos(Fstar(LE)) * cos(F21(LE)) + cos(2 * t12) * sin(Fstar(LE)) * sin(F21(LE)));
}
// after eq 2 and after eq 14
double Dstar()
{
	return (2 * Dmsq31 - Dmsq21);
}
// after eq 13
double F21(double LE)
{
	return 1.267 * LE * Dmsq21;
}
// eq 14
double Fstar(double LE)
{
	return 1.267 * LE * Dstar();
}
// eq 17
double D21M(double a, double LE)
{
	return Dmsq21 + a * cos(2 * t12);
}
// eq 18
double P0M(double a, double LE)
{
	return P0(LE) + a * pow(sin(2 * t12), 2) * cos(2 * t12) * pow(cos(t13), 4) * (1.267 * LE * sin(2 * F21(LE)) - 2 * pow(sin(F21(LE)), 2) / Dmsq21);
}
// eq 18
double PstarM(double a, double LE)
{
	return Pstar(LE) + 0.5 * a * pow(sin(2 * t13), 2) * (1.267 * LE * ((1 + pow(cos(2 * t12), 2)) * sin(Fstar(LE)) * cos(F21(LE)) + 2 * cos(2 * t12) * cos(Fstar(LE)) * sin(F21(LE))) + s12sq * sin(Fstar(LE)) * sin(F21(LE)) / Dmsq21);
}
// above eq 13
// units are in km / GeV as usual
double Pee(double a, double LE)
{
	return 1 - P0M(a, LE) - PstarM(a, LE);
}
} // namespace LWX
