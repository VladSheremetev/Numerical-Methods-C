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

double f(double, double);
double psi0(double);
double psi1(double);
double phi(double);
double l2(double*, size_t);
double l1(double*, size_t);
double l0(double*, size_t);
void matrix_print(double**, size_t, size_t);
double**  grid_shift(double**, double**, size_t, size_t, size_t);
double** explicit_schema(double*, double*, size_t, size_t,  double, double);
double** weighted_schema(double*, double*, size_t, size_t,  double, double, double);
double*  matrix_to_vector(double**, size_t, size_t);
void print_array(double *, size_t);
void solveMatrix (int n, double *, double *, double *, double *, double *);
void print_error(double*, double*, double **, double **, int, int, int, int);


#endif // HEADER_H_INCLUDED
