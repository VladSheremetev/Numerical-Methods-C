#pragma once
#ifndef PRAK4_H
# define PRAK4_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

double   *y_function	   (double*, size_t);
double   *y_derivative	   (double*, size_t);
double   *x_gen			   (double,  size_t, double);
double   *oneh_derivative  (double*, size_t, double);
double   *halfh_derivative (double*, size_t, double);
void     print_array       (double*, size_t);
void 	 print_file        (std::string, double*, double*, size_t);
void     error_calculation (double*, double*, size_t);
double   *runge_calculation(double*, double*, size_t);

#endif
