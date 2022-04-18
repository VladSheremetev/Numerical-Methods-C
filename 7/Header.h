#pragma once
#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED

# include <time.h>
# include <math.h>
# include <iostream>
# include <fstream>
# include <string>
# include <iomanip>


double f1(double, double);
double f2(double, double);
double f_x(double);
double f_y(double);
double f1_deriv_x(double, double);
double f2_deriv_x(double, double);
double f1_deriv_y(double, double);
double f2_deriv_y(double, double);
double l2(double*, size_t);
double l1(double*, size_t);
double l0(double*, size_t);
double* simple_it(double, double);
double* newton(double, double);
double* newton_mod(double, double);
double* newton_disc(double, double);
double* lu(double**, double**, double*, size_t);
void  decompose(double**, double**, double**, size_t);
void matrix_print(double**, size_t, size_t);
double* gauss(double**, double*, size_t);
double f1_deriv_num(double, double, double);
double f2_deriv_num(double, double, double);
void print_array(double*, size_t);


#endif // HEADER_H_INCLUDED
