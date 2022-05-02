#include "Header.h"

#define N 1
#define M 1


double f1(double y1){
	return y1;
}


double f2(double x, double y, double y1){
	return pow(y1, N) * log(x + 2) + pow(y, M) * exp(x) + pow(x, 3);
}


double l0(double* x, size_t n)
{
	double max_elem = 0.0;
	for (size_t i = 0; i < n; ++i)
		if (max_elem <= abs(x[i]))
			max_elem = abs(x[i]);
	return max_elem;
}


double l2(double* x, size_t n)
{
	double sum = 0.0;
	for (size_t i = 0; i < n; ++i)
		sum += pow(abs(x[i]), 2);
	return sqrt(sum);
}


double l1(double* x, size_t n){
    	double sum = 0.0;
	for (size_t i = 0; i < n; ++i)
		sum += abs(x[i]);
	return sum;
}


void matrix_print(double** a, size_t n, size_t m)
{
	for (size_t i = 0; i < n; ++i)
	    {
	        for (size_t j = 0; j < m; ++j)
	            std::cout<<std::setprecision(2)<<a[i][j]<< "\t";
	        std::cout << std::endl;
	    }
	    std::cout << std::endl;
}

double* sub(double* a, double* b, size_t k, size_t m){
    	size_t n = min(m, k);
	double* c = new double [n];
	if (k > m)
		for (size_t i = 0; i < n; ++i)
			c[i] = b[i] - a[2 * i];
	else if (k < m)
		for (size_t i = 0; i < n; ++i)
			c[i] = b[2 * i] - a[i];
	else
		for (size_t i = 0; i < n; ++i)
			c[i] = b[i] - a[i];
	return c;
}


void euler(double* x, double* y, double* y1, size_t n){
	double h = (x[n - 1] - x[0]) / (n - 1);
	for (size_t i = 1; i < n; ++i)
	{
		y[i] = y[i - 1] + h * f1(y1[i - 1]);
		y1[i] = y1[i - 1] + h * f2(x[i - 1], y[i - 1], y1[i - 1]);
	}
}


void euler_predict(double* x, double* y, double* y1, size_t n){
	double h = (x[n - 1] - x[0]) / (n - 1);
	double y_pred, y1_pred;
	double h2 = h * 0.5;
	for (size_t i = 1; i < n; ++i)
	{
		y_pred = y[i - 1] + h * f1(y1[i - 1]);
		y1_pred = y1[i - 1] + h * f2(x[i - 1], y[i - 1], y1[i - 1]);
		y[i] = y[i - 1] + h2 * (f1(y1[i - 1]) + f1(y1_pred));
		y1[i] = y1[i - 1] + h2 * (f2(x[i - 1], y[i - 1], y1[i - 1]) + f2(x[i], y_pred, y1_pred));
	}
}


void runge_kutt_2(double* x, double* y, double* y1, size_t n){
	double h = (x[n - 1] - x[0]) / (n - 1);
	double k1, k2;					
	double k11, k21;	
	for (size_t i = 1; i < n; ++i)
	{
	    k1 = h*f1(y1[i - 1]);
	    k11 = h * f2(x[i - 1], y[i - 1], y1[i - 1]);
	    k2 = h*f1(y1[i - 1] + (2*k11)/3);
	    k21 = h * f2(x[i - 1] + (2*h)/3, y[i - 1] + (2*k1)/3, y1[i - 1] + (2*k11)/3);
	    y[i] = y[i - 1] + k1*0.25 + k2*0.75;
	    y1[i] = y1[i - 1] + k11*0.25 + k21*0.75;
	}
}


void runge_kutt_4(double* x, double* y, double* y1, size_t n){
	double h = (x[n - 1] - x[0]) / (n - 1);
	double h2 = h / 2.0;
	double h6 = h / 6.0;
	double k1, k2, k3, k4;					
	double k11, k21, k31, k41;
	for (size_t i = 1; i < n; ++i)
	{
		k1 = f1(y1[i - 1]);
		k11 = f2(x[i - 1], y[i - 1], y1[i - 1]);
		k2 = f1(y1[i - 1] + h2 * k11);
		k21 = f2(x[i - 1] + h2, y[i - 1] + h2 * k1, y1[i - 1] + h2 * k11);
		k3 = f1(y1[i - 1] + h2 * k21);
		k31 = f2(x[i - 1] + h2, y[i - 1] + h2 * k2, y1[i - 1] + h2 * k21);
		k4 = f1(y1[i - 1] + h * k31);
		k41 = f2(x[i - 1] + h, y[i - 1] + h * k3, y1[i - 1] + h * k31);
		y[i] = y[i - 1] + h6 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
		y1[i] = y1[i - 1] + h6 * (k11 + 2.0 * k21 + 2.0 * k31 + k41);
	}
}


void adams_3(double* x, double* y, double* y1, size_t n){
	double h = (x[n - 1] - x[0]) / (n - 1);
	double k1 = 0.0, k2 = 0.0, k3 = 0.0;
	double k11 = 0.0, k21 = 0.0, k31 = 0.0;
	double delta_y, delta_y1;

	runge_kutt_4(x, y, y1, 3);

	for (size_t i = 3; i < n; ++i)
	{
		k1 = h * f1(y1[i - 1]);
		k11 = h * f2(x[i - 1], y[i - 1], y1[i - 1]);
		k2 = h * f1(y1[i - 2]);
        	k21 = h * f2(x[i - 2], y[i - 2], y1[i - 2]);
		k3 = h * f1(y1[i - 3]);
        	k31 = h * f2(x[i - 3], y[i - 3], y1[i - 3]);
		delta_y = (23.0 * k1 - 16.0 * k2 + 5.0 * k3) / 12.0;
		delta_y1 = (23.0 * k11 - 16.0 * k21 + 5.0 * k31) / 12.0;
		y[i] = y[i - 1] + delta_y;
		y1[i] = y1[i - 1] + delta_y1;
	}
}
