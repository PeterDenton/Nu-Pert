#ifndef First_Peak_H
#define First_Peak_H

double PL(double L, const double * p);
// at fixed E
// max = true calculates the max, max = false calculates the min
double First_Extremum_Exact(int alpha, int beta, double E, double Yrho, double delta, bool max); // fixed E

double FEE_x(double E, const double *p);
double FEE_y(double E, const double *p);

// finds the slope and y-int in the L vs L/E plane for a given CP phase
void FP_line(int alpha, int beta, double Yrho, double delta, bool max, double *m, double *b);

double LE_SP0(double delta, int order); // SP's solution attempt at L/E in vac
double LE_SP1(double delta, int order); // SP's solution attempt at L/E in vac

#endif
