#ifndef Linear_Fit_H
#define Linear_Fit_H

template <class T>
void Linear_Fit(T (*x)(T t, const T *p), T (*y)(T t, const T *p), T t_min, T t_max, int n, T *p, T *m, T *b, T *rsq);

#endif
