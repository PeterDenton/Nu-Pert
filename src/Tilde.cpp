#include "Tilde.h"
#include "Constants.h"

namespace Tilde
{
// See eq. 2.2.5
double la(double a)
{
	return a + (s13sq + eps * s12sq) * Dmsqee;
}
// See eq. 2.2.5
double lb(double a)
{
	return eps * c12sq * Dmsqee;
}
// See eq. 2.2.5
double lc(double a)
{
	return (c13sq + eps * s12sq) * Dmsqee;
}
} // namespace Tilde
