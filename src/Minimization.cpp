#include <cmath>

#include "Minimization.h"
#include "Warning.h"

#include <iostream>

// from https://en.wikipedia.org/wiki/Golden_section_search
// a<b
// p is the parameters to pass into f
template <class T>
T gss_min(T (*f)(T, const T*), T a_initial, T b_initial, T tol, T *p)
{
	T gr = (sqrt(5.) - 1) / 2;
	T a = a_initial;
	T b = b_initial;
	T c = b - gr * (b - a);
	T d = a + gr * (b - a);
	while (std::abs(c - d) > tol)
	{
		if (f(c, p) < f(d, p))
		{
			b = d;
			d = c;
			c = b - gr * (b - a);
		}
		else
		{
			a = c;
			c = d;
			d = a + gr * (b - a);
		}
	}
	T ret = (b + a) / 2.;

	// Warns user in case final minimum is close to an edge
	if (std::abs(ret - a_initial) < 2 * tol or std::abs(b_initial - ret) < 2 * tol)
		EdgeCaseWarning W(a_initial, b_initial, ret);
	return ret;
}

template double gss_min<double>(double (*)(double, const double*), double, double, double, double*);

