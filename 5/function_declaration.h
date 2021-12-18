#pragma once
#ifndef function_declaration_h
#define function_declaration_h

#include <math.h>
#include <iostream>

double  *x_gen		   (double,  size_t, double);
double  *yfun		   (double*, size_t);
double  yintegral	   (double,  double);
double  square_formula (double*, size_t, double);
double  trapeze_formula(double*, size_t, double);
double  simpson_formula(double*, size_t, double);
double  newton_formula (double*, size_t, double);
double  gauss_formula  (double*, size_t, double);

#endif
