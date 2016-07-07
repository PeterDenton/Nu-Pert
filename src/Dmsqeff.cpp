/* See arXiv:1601.07464 which pulls from arXiv:hep-ph/0503283 */

#include "Dmsqeff.h"
#include "Hat.h"
#include "Constants.h"

namespace Dmsqeff
{
double Dlee(double a)
{
	return Hat::Dlpm(a);
}
double Dlmumu(double a, double delta)
{
	return Hat::Dlp0(a) + cos(delta) * sin(Hat::phi(a)) * cos(Hat::phi(a) - t13) * sin(2 * t12) * tan(t23) * Dmsq21;
}
double Pee(double a, double LE)
{
	double L4E = km_per_GeV_to_per_eV2 * LE / 4; // in eV^-2 now

	double Delta21 = Dmsq21 * L4E;
	double Deltaee = Dlee(a) * L4E;

	return 1 - pow(c13sq * sin(2 * t12) * sin(Delta21), 2) - pow(sin(2 * t13) * sin(Deltaee), 2);
}
} // namespace Dmsqeff
