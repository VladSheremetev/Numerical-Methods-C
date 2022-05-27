#include "Header.h"

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


double* zeros(size_t n){
    double* x = new double[n];
    for (size_t i = 0; i < n; ++i) x[i] = 0.0;
    return x;
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


double* vector_shift(double* a, double* b, size_t n){
	double* c = new double [n];
    for (size_t i = 0; i < n; ++i)
        c[i] = b[2 * i] - a[i];
	return c;
}


double f1(double y1){
	return y1;
}

double f2(double x, double y, double y1) {
	return y1 * log(x + 2) + y * exp(x) + pow(x, 3);
}

double f2_nonlinear(double x, double y, double y1){
	return y1 * y1 * log(x + 2) + y * y * exp(x) + pow(x, 3);
}

double f2_temp(double x, double y, double y1){
	return y1 * log(x + 2) + y * exp(x);
}

double f2_temp2(double x, double y, double y1, double u1, double u2){
	return f2_dy(x,y,y1) * u1 + f2_dy1(x,y,y1) * u2;
}

double p(double x){
	return log(x + 2);
}

double q(double x){
	return exp(x);
}

double r(double x){
	return pow(x, 3);
}

double f2_dy(double x, double y, double y1){
	return 2 * y * exp(x);
}

double f2_dy1(double x, double y, double y1){
	return 2 * y1 * log(x + 2);
}

void ballistic(double* x, double* y, size_t n){
	double* ERK4_1  = new double [n];
	double* ERK4_2  = new double [n];
	double* temp_y1 = new double [n];

    //начальные условия для подзадачи
	temp_y1[0] = 0;
	ERK4_1[0] = y[0];

	rk4(x, ERK4_1, temp_y1, n);

	//начальные условия для подзадачи
	temp_y1[0] = 1;
	ERK4_2[0] = 0;

	rk4_tmp(x, ERK4_2, temp_y1, n);

	double koef = (y[n - 1] - ERK4_1[n - 1]) / ERK4_2[n - 1];
	for (size_t i = 1; i < n; ++i)
		y[i] = ERK4_1[i] + koef * ERK4_2[i];

    delete [] ERK4_1;
    delete [] ERK4_2;
    delete [] temp_y1;
}

void finite_subs(double* x, double* y, size_t n){
	double h = (x[n - 1] - x[0]) / (n - 1);
	double* gf = new double [n];
	double* a1 = new double [n]; //верхняя диагональ
	double* a2 = new double [n]; //средняя диагональ
	double* a3 = new double [n]; //нижняя диагональ

	for (size_t i = 1; i < n; ++i){
		a2[i] = - (4.0 + h * h * 2.0 * q(x[i]));
		a3[i] = 2.0 + h * p(x[i]);
		a1[i - 1] = 2.0 - h * p(x[i]);
	}

    a2[0] = 1.0; a2[n - 1] = 1.0;
	a1[0] = 0.0; a3[n - 1] = 0.0;
	gf[0] = y[0];
	gf[n - 1] = y[n - 1];

	for (size_t i = 1; i < n - 1; ++i)
        gf[i] = 2 * r(x[i]) * h * h;

	solveMatrix(n, a3, a2, a1, gf, y); //метод прогонки для 3х диагональной матрицы

    delete [] gf;
    delete [] a1;
    delete [] a2;
    delete [] a3;
}

void print_array(double *arr, size_t n)
{
    for (size_t i = 0; i < n; ++i)
        std::cout << arr[i] << " ";
    std::cout<<std::endl;
}

//метод прогонки
/*
	 * n - число уравнений (строк матрицы)
	 * b - диагональ, лежащая над главной (нумеруется: [0;n-2])
	 * c - главная диагональ матрицы A (нумеруется: [0;n-1])
	 * a - диагональ, лежащая под главной (нумеруется: [1;n-1])
	 * f - правая часть (столбец)
	 * x - решение, массив x будет содержать ответ
*/
void solveMatrix (int n, double *a, double *c, double *b, double *f, double *x)
{
	double m;
	for (size_t i = 1; i < n; ++i)
	{
		m = a[i] / c[i-1];
		c[i] = c[i] - m*b[i-1];
		f[i] = f[i] - m*f[i-1];
	}

	x[n-1] = f[n-1] / c[n-1];

	for (int i = n - 2; i >= 0; --i)
    {
		x[i] = (f[i] - b[i]*x[i+1]) / c[i];
    }
}


void newton(double* x, double* y, size_t n, double q0, double q1, double a, double b){
	double s = (q1 - q0) / (b - a), delta = 100.0;
	double eps = 0.0001;
	size_t num_it = 0;
	size_t MAX_ITER = 500;
	double* y1 = new double [n];
	double* u1 = new double [n];
	double* u2 = new double [n];

	for ( ; num_it < MAX_ITER; ++num_it){
        //начальные условия
        y1[0] = s;
        u1[0] = 0;
		u2[0] = 1;

		rk4_nonlin(x, y, y1, n);
		rk4_tmp2(x, y, y1, u1, u2, n);

		delta = (y[n - 1] - q1) / u1[n - 1];
		if (fabs(delta) < eps) break;

		s -= delta;
	}

	std::cout << "newton_num_it=" << num_it << std::endl;

	delete [] y1;
	delete [] u1;
	delete [] u2;
}


void chordes(double* x, double* y, size_t n){
	double u1 = 0.0, u2 = 0.5;
	double* ERK4_y11 = new double [n];
	double* ERK4_y12 = new double [n];
	double* temp_y = new double [n];
	size_t MAX_ITER = 500;
	double q1 = y[n - 1];
	double eps = 1e-4, temp;
	temp_y[0] = y[0];
	size_t num_it = 0;

	for ( ;num_it < MAX_ITER; ++num_it){
		ERK4_y11[0] = u1;
		ERK4_y12[0] = u2;

		rk4_nonlin(x, y, ERK4_y12, n);
		rk4_nonlin(x, temp_y, ERK4_y11, n);
		if (fabs((y[n - 1] - q1) / q1) < eps) break;
		temp = u2;
		u2 -= (u2 - u1) * (y[n - 1] - q1) / (y[n - 1] - temp_y[n - 1]);
		u1 = temp;
	}
	std::cout << "secant_num_it=" << num_it << std::endl;

	delete [] ERK4_y11;
	delete [] ERK4_y12;
	delete [] temp_y;
}

void rk4(double* x, double* y, double* y1, size_t n)
{
	double h = (x[n - 1] - x[0]) / (n - 1);
	double h2 = h * 0.5, h6 = h / 6.0;
	double k1, k2, k3, k4;
	double k11, k21, k31, k41;
	for (size_t i = 1; i < n; ++i)
	{
		k1  = f1(y1[i - 1]);
		k11 = f2(x[i - 1], y[i - 1], y1[i - 1]);
		k2  = f1(y1[i - 1] + h2 * k11);
		k21 = f2(x[i - 1]  + h2, y[i - 1] + h2 * k1, y1[i - 1] + h2 * k11);
		k3  = f1(y1[i - 1] + h2 * k21);
		k31 = f2(x[i - 1]  + h2, y[i - 1] + h2 * k2, y1[i - 1] + h2 * k21);
		k4  = f1(y1[i - 1] + h * k31);
		k41 = f2(x[i - 1]  + h, y[i - 1] + h * k3, y1[i - 1] + h * k31);
		y[i]  = y[i - 1]  + h6 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
		y1[i] = y1[i - 1] + h6 * (k11 + 2.0 * k21 + 2.0 * k31 + k41);
	}
}

void rk4_tmp(double* x, double* y, double* y1, size_t n){
	double h = (x[n - 1] - x[0]) / (n - 1);
	double h2 = h * 0.5, h6 = h / 6.0;
	double k1, k2, k3, k4;
	double k11, k21, k31, k41;
	for (size_t i = 1; i < n; ++i)
	{
		k1  = f1(y1[i - 1]);
		k11 = f2_temp(x[i - 1], y[i - 1], y1[i - 1]);
		k2  = f1(y1[i - 1] + h2 * k11);
		k21 = f2_temp(x[i - 1] + h2, y[i - 1] + h2 * k1, y1[i - 1] + h2 * k11);
		k3  = f1(y1[i - 1] + h2 * k21);
		k31 = f2_temp(x[i - 1] + h2, y[i - 1] + h2 * k2, y1[i - 1] + h2 * k21);
		k4  = f1(y1[i - 1] + h * k31);
		k41 = f2_temp(x[i - 1] + h, y[i - 1] + h * k3, y1[i - 1] + h * k31);
		y[i]  = y[i - 1] + h6 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
		y1[i] = y1[i - 1] + h6 * (k11 + 2.0 * k21 + 2.0 * k31 + k41);
	}
}


void rk4_nonlin(double* x, double* y, double* y1, size_t n)
{
	double h = (x[n - 1] - x[0]) / (n - 1);
	double h2 = h * 0.5, h6 = h / 6;
	double k1, k2, k3, k4;
	double k11, k21, k31, k41;
	for (size_t i = 1; i < n; ++i)
	{
		k1  = f1(y1[i - 1]);
		k11 = f2_nonlinear(x[i - 1], y[i - 1], y1[i - 1]);
		k2  = f1(y1[i - 1] + h2 * k11);
		k21 = f2_nonlinear(x[i - 1] + h2, y[i - 1] + h2 * k1, y1[i - 1] + h2 * k11);
		k3  = f1(y1[i - 1] + h2 * k21);
		k31 = f2_nonlinear(x[i - 1] + h2, y[i - 1] + h2 * k2, y1[i - 1] + h2 * k21);
		k4  = f1(y1[i - 1] + h*k31);
		k41 = f2_nonlinear(x[i - 1] + h, y[i - 1] + h*k3, y1[i - 1] + h*k31);
		y[i]  = y[i - 1] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * h6;
		y1[i] = y1[i - 1] + (k11 + 2.0 * k21 + 2.0 * k31 + k41) * h6;
	}
}

void rk4_tmp2(double* x, double* y, double* y1, double* u1, double* u2, size_t n)
{
	double h = (x[n - 1] - x[0]) / (n - 1);
	double h2 = h * 0.5, h6 = h/6;
	double k1, k2, k3, k4;
	double k11, k21, k31, k41;
	for (size_t i = 1; i < n; ++i)
	{
		k1  = f1(u2[i - 1]);
		k11 = f2_temp2(x[i - 1], y[i - 1], y1[i - 1], u1[i - 1], u2[i - 1]);
		k2  = f1(u2[i - 1] + h2 * k11);
		k21 = f2_temp2(x[i - 1] + h2, y[i - 1], y1[i - 1], u1[i - 1] + k1 * h2, u2[i - 1] + k11 * h2);
		k3  = f1(u2[i - 1] + h2 * k21);
		k31 = f2_temp2(x[i - 1] + h2, y[i - 1], y1[i - 1], u1[i - 1] + k2 * h2, u2[i - 1] + k21 * h2);
		k4  = f1(u2[i - 1] + h*k31);
		k41 = f2_temp2(x[i - 1] + h, y[i - 1], y1[i - 1], u1[i - 1] + h*k3, u2[i - 1] + h*k31);
		u1[i] = u1[i - 1] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * h6;
		u2[i] = u2[i - 1] + (k11 + 2.0 * k21 + 2.0 * k31 + k41) * h6;
	}
}
