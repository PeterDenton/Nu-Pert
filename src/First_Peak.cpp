#include <math.h>
#include <iostream>
#include <fstream>
#include <assert.h>

#include "First_Peak.h"
#include "Exact.h"
#include "Constants.h"
#include "Experiment.h"
#include "Hat.h"
#include "Check.h"
#include "GF.h"
#include "Minimization.h"
#include "Progress.h"
#include "Linear_Fit.h"

double tol = 1e-5;

// P as a function of L, inverted for appearance so that the first max is always a minimum
// params p: alpha, beta, E, Yrho, delta
double PL(double L, const double * p)
{
	int alpha = p[0];
	int beta = p[1];
	double E = p[2];
	double Yrho = p[3];
	double delta = p[4];
	int sign = p[5];
	return sign * Exact::Palphabeta(alpha, beta, Yrho * E * Y_to_a, L / E, delta);
}
double First_Extremum_Exact(int alpha, int beta, double E, double Yrho, double delta, bool max)
{
	double L_range_min, L_range_max;
	double p[6];
	p[0] = alpha;
	p[1] = beta;
	p[2] = E;
	p[3] = Yrho;
	p[4] = delta;
	p[5] = (alpha == beta) ? 1 : -1;
	p[5] *= (max) ? 1 : -1; // switch sign for max/min

	// numbers here are selected so that the first max/min are contained in this region, but not the second
	if (max)
	{
		L_range_min = E * 300;
		L_range_max = E * 1000;
	}
	else
	{
		L_range_min = E * 500;
		L_range_max = E * 1800;
	}

	return std::abs(gss_min(PL, L_range_min, L_range_max, tol, p));
}
double FEE_x(double E, const double *p)
{
	int alpha = p[0];
	int beta = p[1];
	double Yrho = p[2];
	double delta = p[3];
	bool max = p[4] > 0;

	return First_Extremum_Exact(alpha, beta, E, Yrho, delta, max) / E;
}
double FEE_y(double E, const double *p)
{
	int alpha = p[0];
	int beta = p[1];
	double Yrho = p[2];
	double delta = p[3];
	bool max = p[4] > 0;

	return First_Extremum_Exact(alpha, beta, E, Yrho, delta, max);
}
void FP_line(int alpha, int beta, double Yrho, double delta, bool max, double *m, double *b)
{
	double p[5];
	p[0] = alpha;
	p[1] = beta;
	p[2] = Yrho;
	p[3] = delta;
	p[4] = 1; // look for a max

	double rsq;
	Linear_Fit(FEE_x, FEE_y, 0.1, 5., 1e2, p, m, b, &rsq);
}

double LE_SP(double delta)
{
	double k = 8 * Jr * eps / (s23sq * pow(sin(2 * t13), 2));
	double sd = sin(delta);
	double cd = cos(delta);
	return (4 / (Dmsqee + Jr * Dmsq21)) * ((M_PI / 2) - (k / 2) * (sd + (M_PI / 2) * cd) * (1 - k * (cd - (M_PI / 2) * sd))) / km_per_GeV_to_per_eV2;
}
