#include <cmath>

#include "Linear_Fit.h"

#include <iostream>

template <class T>
void Linear_Fit(T (*x)(T t, const T *p), T (*y)(T t, const T *p), T t_min, T t_max, int n, T *p, T *m, T *b, T *rsq)
{
	T sum_x = 0;
	T sum_xsq = 0;
	T sum_y = 0;
	T sum_ysq = 0;
	T sum_xy = 0;

	T x_tmp, y_tmp;
	for (T t = t_min; t <= t_max; t += (t_max - t_min) / (n - 1))
	{
		x_tmp = x(t, p);
		y_tmp = y(t, p);
		sum_x += x_tmp;
		sum_xsq += pow(x_tmp, 2);
		sum_y += y_tmp;
		sum_ysq += pow(y_tmp, 2);
		sum_xy += x_tmp * y_tmp;
	}

	T x_bar = sum_x / n;
	T y_bar = sum_y / n;

	T ss_xx = sum_xsq - n * pow(x_bar, 2);
	T ss_yy = sum_ysq - n * pow(y_bar, 2);
	T ss_xy = sum_xy - n * x_bar * y_bar;

	*m = ss_xy / ss_xx;
	*b = y_bar - *m * x_bar;
	*rsq = pow(ss_xy, 2) / (ss_xx * ss_yy);
}

template void Linear_Fit(double (*x)(double t, const double *p), double (*y)(double t, const double *p), double t_min, double t_max, int n, double *p, double *m, double *b, double *rsq);
