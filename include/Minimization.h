#ifndef Minimization_H
#define Minimization_H

template <class T>
T gss_min(T (*f)(T, const T*), T a, T b, T tol, T *p);

#endif
