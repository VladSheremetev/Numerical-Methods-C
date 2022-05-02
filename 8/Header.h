#pragma once

#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED

# include <time.h>
# include <math.h>
# include <iostream>
# include <fstream>
# include <string>
# include <iomanip>

using namespace std;

double f1(double);
double f2(double, double, double);
double l2(double*, size_t);
double l1(double*, size_t);
double l0(double*, size_t);
void euler(double*, double*, double*, size_t);
double* sub(double*, double*, size_t, size_t);
void euler_predict(double*, double*, double*, size_t);
void runge_kutt_2(double*, double*, double*, size_t);
void runge_kutt_4(double*, double*, double*, size_t);
void adams_3(double*, double*, double*, size_t);
void matrix_print(double**, size_t, size_t);


#endif
