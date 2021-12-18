#pragma once
#ifndef function_declaration_h
#define function_declaration_h

#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

double*		x_gen(double, size_t, double);
double*   	y_gen(double*, size_t);
double* 	l_gen(double*, double*, double*, size_t, size_t);
double      lagrange(double*, double*, size_t, double);
void        print_array(double*, size_t);
void        error_calculation(double*, double*, size_t);
void 		print_file(std::string, double*, double*, size_t);
int 		count_step_between_point(double, double, double);

#endif
