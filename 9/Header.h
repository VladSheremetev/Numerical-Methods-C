#pragma once
#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED

# include <time.h>
# include <math.h>
# include <iostream>
# include <fstream>
# include <string>
# include <iomanip>

//using namespace std;

double f1(double);
double f2(double, double, double);
double f2_temp(double, double, double);
double p(double);
double q(double);
double r(double);
double f2_dy(double, double, double);
double f2_dy1(double, double, double);
double f2_temp2(double, double, double, double, double);
double l2(double*, size_t);
double l1(double*, size_t);
double l0(double*, size_t);
void matrix_print(double**, size_t, size_t);
double* vector_shift(double*, double*, size_t);
void rk4(double*, double*, double*, size_t);
void rk4_tmp(double*, double*, double*, size_t);
void rk4_tmp2(double*, double*, double*, double*, double*, size_t);
void ballistic(double*, double*, size_t);
void newton(double*, double*, size_t, double, double, double, double);
void chordes(double*, double*, size_t);
void finite_subs(double*, double*, size_t);
double f2_nonlinear(double, double, double);
void rk4_nonlin(double*, double*, double*, size_t);
double* zeros(size_t);
void solveMatrix (int, double*, double *, double *, double *, double *);
void print_array(double*, size_t);


#endif // HEADER_H_INCLUDED
