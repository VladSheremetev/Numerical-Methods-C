#include "Header.h"
#define _USE_MATH_DEFINES

double f(double x, double t){
	return log(1 + t) * pow(x, 3);
}

double psi1(double t){
	return 2 * sin(M_PI*t);
}

double psi0(double t){
	return 0.5 - cos(M_PI*t);
}

double phi(double x){
	return 4*pow(x - 0.5, 2) - 1;
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


double** grid_shift(double** u1, double** u2, size_t n_x0, size_t n_t, size_t n_x1){
    size_t n = n_x0*n_t;
	double** c = new double* [n_x0];
	for (size_t i = 0; i < n_x0; ++i) c[i] = new double [n_t];
    for (size_t j = 0; j < n_t; ++j){
        for (size_t i = 0; i < n_x0; ++i)
        c[i][j] = u2[(2 * i)][j] - u1[i][j];
    }
	return c;
}

void print_error(double* err1h1, double* err2h1, double **u_1_h1, double **u_2_h1, int n_x0, int n_t, int n_x01, int n_t2){
    std::cout << std::setw(10)<< "absolute error:\t\t\t" << "l2:\t\t" << "l1:\t\t" << "l0:" << std::endl;
	std::cout <<std::scientific<< "explicit_schema:\t\t" << l2(err1h1, n_x0*n_t) << '\t' << l1(err1h1, n_x0*n_t) << '\t' << l0(err1h1, n_x0*n_t) << std::endl;
	std::cout << "weighted_schema:\t\t" << l2(err2h1, n_x01*n_t2) << '\t' << l1(err2h1, n_x01*n_t2) << '\t' << l0(err2h1, n_x01*n_t2) << std::endl;
	std::cout << "relative error:" << std::endl;
	std::cout << "explicit_schema:\t\t" << l2(err1h1, n_x0*n_t) / l2(matrix_to_vector(u_1_h1, n_x0, n_t), n_x0*n_t) << '\t' << l1(err1h1, n_x0*n_t) / l1(matrix_to_vector(u_1_h1, n_x0, n_t), n_x0*n_t) << '\t' << l0(err1h1, n_x0*n_t) / l0(matrix_to_vector(u_1_h1, n_x0, n_t), n_x0*n_t) << std::endl;
	std::cout << "weighted_schema:\t\t" << l2(err2h1, n_x01*n_t2) / l2(matrix_to_vector(u_2_h1, n_x01, n_t2), n_x01*n_t2) << '\t' << l1(err2h1, n_x01*n_t2) / l1(matrix_to_vector(u_2_h1, n_x01, n_t2), n_x01*n_t2) << '\t' << l0(err2h1, n_x01*n_t2) / l0(matrix_to_vector(u_2_h1, n_x01, n_t2), n_x01*n_t2) << std::endl;

}

double* matrix_to_vector(double** a, size_t n, size_t m) {
	double* res = new double [n * m];
	for (size_t j = 0; j < m; ++j)
		for (size_t i = 0; i < n; ++i)
		{
			res[j * n + i] = a[i][j];
		}
	return res;
}

double** explicit_schema(double* x, double* t, size_t n_x, size_t n_t, double h_x, double h_t){
	double th = h_t / (h_x * h_x);
	double** u =  new double* [n_x];
	for (size_t i = 0; i < n_x; ++i) u[i] = new double [n_t];
	for (size_t i = 0; i < n_x; ++i){
        for (size_t j = 0; j < n_t; ++j) u[i][j] = 0.0;
	}

	for (size_t i = 0; i < n_x; ++i)
		u[i][0] = phi(x[i]);

	for (size_t j = 0; j < n_t - 1; ++j)
	{
		for (size_t i = 1; i < n_x - 1; ++i)
		{
			u[i][j + 1] = u[i][j] + th * (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) + h_t * f(x[i], t[j]);
		}
		u[0][j + 1] = psi0(t[j + 1]);
		u[n_x - 1][j + 1] = psi1(t[j + 1]);
	}

	return u;
}

void print_array(double *arr, size_t n)
{
    for (size_t i = 0; i < n; ++i)
        std::cout <<std::scientific<< arr[i] << " ";
    std::cout<<std::endl;
}

double** weighted_schema(double* x, double* t, size_t n_x, size_t n_t, double sigma, double h_x, double h_t){

	double th = h_t / (h_x * h_x), l = sigma * th, c = -(1 + 2 * sigma * th);

	double** u =  new double* [n_x];
	for (size_t i = 0; i < n_x; ++i) u[i] = new double [n_t];
	for (size_t i = 0; i < n_x; ++i){
        for (size_t j = 0; j < n_t; ++j) u[i][j] = 0.0;
	}

	for (int i = 0; i < n_x; ++i)
		u[i][0] = phi(x[i]);

	double* gf = new double [n_x];
	double* temp = new double [n_x];

	double* a1 = new double [n_x]; //верхняя диагональ
	double* a2 = new double [n_x]; //средняя диагональ
	double* a3 = new double [n_x]; //нижняя диагональ

	for (size_t i = 0; i < n_x; i++){
		a2[i] = c;
        a3[i] = l;
		a1[i] = l;
    }

	a2[0] = 1.0;
	a2[n_x - 1] = 1.0;
	a1[0] = 0.0;
	a3[n_x - 1] = 0.0;

	double *a1_tmp = new double [n_x];
    double *a2_tmp = new double [n_x];
    double *a3_tmp = new double [n_x];

	for (size_t j = 0; j < n_t - 1; ++j)
	{
	    for (size_t i = 0; i < n_x; ++i){
            a1_tmp[i] = a1[i];
            a2_tmp[i] = a2[i];
            a3_tmp[i] = a3[i];
	    }

		for (size_t i = 1; i < n_x - 1; ++i)
			gf[i] = -(u[i][j] + (1 - sigma) * th * (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) + h_t * f(x[i], t[j]));
		gf[0] = u[0][j];
		gf[n_x - 1] = u[n_x - 1][j];
        solveMatrix(n_x, a3_tmp, a2_tmp, a1_tmp, gf, temp);
		for (size_t i = 1; i < n_x - 1; ++i)
		{
			u[i][j + 1] = temp[i];
		}
		u[0][j + 1] = psi0(t[j + 1]);
		u[n_x - 1][j + 1] = psi1(t[j + 1]);
	}

	delete [] a1_tmp;
	delete [] a2_tmp;
	delete [] a3_tmp;
	delete [] gf;
	delete [] temp;
	delete [] a1;
	delete [] a2;
	delete [] a3;

	return u;
}

/*
     * Метод прогонки
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
