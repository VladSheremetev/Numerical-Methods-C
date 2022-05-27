#pragma once
#ifndef PRAK6_H_INCLUDED
#define PRAK6_H_INCLUDED

# include <time.h>
# include <math.h>
# include <iostream>
# include <Eigen/Dense>
# include <fstream>
# include <string>
# include <iomanip>

using namespace Eigen;

double*		x_gen			 	(double,  size_t, double);
double* 	y_gen			 	(double*, size_t);
void        print_array		 	(double*, size_t);
void        error_calculation	(double*, double*,size_t);
double      *x_rand_gen		 	(double,  double, size_t);
bool 		search_point	 	(double*, size_t, double);
double      *lagrangebasic	 	(double*, size_t, double,  size_t);
int         slau			 	(double*, size_t, size_t,  size_t);
double      scalar_multiplication(double*, double*, size_t);
void        matrix_multiplication(double*, double**, double*, size_t, int);
void        matrix_print        (double**, size_t, size_t);
double*     gauss               (double**, double*, size_t, int);
double*     lu                  (double**, double*, size_t);
void        decompose           (double**, double**, double**, size_t, int);
double*     relax               (double**, double*, double, size_t, int);
double      l0                  (double* , size_t);
double*     grad                (double**, double*, size_t, int);
double      l2                  (double* , size_t);
double      l1                  (double* , size_t);
void        error_nev           (double* , double* , size_t);


#endif // PRAK6_H_INCLUDED
