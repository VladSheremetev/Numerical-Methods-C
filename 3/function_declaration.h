#pragma once
#ifndef function_declaration_h
#define function_declaration_h

# include <time.h>
# include <math.h>
# include <iostream>
# include <Eigen/Dense>
# include <fstream>
# include <string>

using namespace Eigen;

double*		x_gen			 	(double,  size_t, double);
double* 	y_gen			 	(double*, size_t);
void        quick_sort		 	(double*, size_t, size_t);
void        print_array		 	(double*, size_t);
void        error_calculation	(double*, double*,size_t);
double      *x_rand_gen		 	(double,  double, size_t);
bool 		search_point	 	(double*, size_t, double);
double      *lagrangebasic	 	(double*, size_t, double,  size_t);
double      *best_fit_func_rand	(double*, size_t, double*, size_t,  size_t, double*);
double      *best_fit_func		(double*, size_t, double*, size_t,  size_t, double*);
void 	     print_file		 	(std::string,     double*, double*, size_t);
double*      slau			 	(double*, size_t, size_t,  size_t,  double*, double*);

#endif
