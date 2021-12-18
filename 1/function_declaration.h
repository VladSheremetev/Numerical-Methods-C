#pragma once
#ifndef function_declaration_h
#define function_declaration_h

#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

double*     x_rand_gen(double, double, size_t);
double*     y_gen(double*, size_t);
void        print_array(double*, size_t);
double  	lagrange(double*, double*, size_t, double);
double*  	x_gen(double, size_t, double);
double*   	l_gen(double*, double*, double*, size_t, size_t);
void 		print_file(std::string, double*, double*, size_t);
bool 		search_point(double*, size_t, double);
void 		quick_sort(double*, size_t, size_t);

#endif
