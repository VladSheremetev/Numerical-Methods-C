#include "Header.h"

double f1(double x, double y){
	return x - pow(y, 3);
}


//печать вектора (массива)
void print_array(double *arr, size_t n)
{
    for(size_t i = 0; i < n; ++i) std::cout << arr[i] << " ";
        std::cout<<std::endl;
}


double f2(double x, double y){
	return tan(x) - y;
}

/*
double f2_y(double x) {
	return tan(x);
}

double f1_x(double y) {
	return pow(y, 3);
}
*/

double f2_y(double x) {
	return pow(x, 1.0/3.0);
}


double f1_x(double y) {
	return atan(y);
}


double f1_deriv_x(double x, double y){
	return 1.0;
}


double f1_deriv_y(double x, double y){
	return -3.0 * pow(y, 2);
}


double f2_deriv_x(double x, double y){
	return (1.0/(pow(cos(x), 2)));
}


double f2_deriv_y(double x, double y){
	return -1.0;
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


double* simple_it(double x0, double y0){
	int num_it = 0;
	double eps = 1e-10;
	double* err = new double [2];
	double* iter_xy = new double [2];
	double* ans = new double [3];

	iter_xy[0] = x0, iter_xy[1] = y0;

	for ( ; ; )
	{
		++num_it;

		if (num_it > 1000)
			break;

		ans[0] = f1_x(iter_xy[1]);
		ans[1] = f2_y(iter_xy[0]);
		err[0] = iter_xy[0] - ans[0];
		err[1] = iter_xy[1] - ans[1];

        if (l2(err, 2) / l2(ans, 2) < eps) break;

		iter_xy[0] = ans[0];
		iter_xy[1] = ans[1];

	}
	ans[2] = num_it;

	delete [] iter_xy;
	delete [] err;

	return ans;
}


double* newton(double x0, double y0){
	double** w = new double* [2];
	for (int i = 0; i < 2; ++i) {w[i] = new double [2];}

	double* gf    = new double [2];
	double* delta = new double [2];
	double* ans   = new double [3];

	int num_it = 0;
	double eps = 1e-3;

	ans[0] = x0, ans[1] = y0;
	for ( ; ; )
	{
		++num_it;
		if (num_it > 1000)
			break;

		w[0][0] = f1_deriv_x(ans[0], ans[1]);
		w[0][1] = f1_deriv_y(ans[0], ans[1]);
		w[1][0] = f2_deriv_x(ans[0], ans[1]);
		w[1][1] = f2_deriv_y(ans[0], ans[1]);
		gf[0] = f1(ans[0], ans[1]);
		gf[1] = f2(ans[0], ans[1]);
		delta = gauss(w, gf, 2);
		ans[0] -= delta[0];
		ans[1] -= delta[1];
		if (l2(delta, 2) / l2(ans, 2) < eps) break;
	}

	for (size_t i = 0; i < 2; ++i) {delete [] w[i];}
	delete [] gf;
	delete [] delta;
	delete [] w;

	ans[2] = num_it;
	return ans;
}


double* newton_mod(double x0, double y0){
	double** w = new double* [2];
	double** l = new double* [2];
	double** u = new double* [2];
	for (int i = 0; i < 2; ++i) {w[i] = new double [2]; l[i] = new double [2]; u[i] = new double [2];}

	//зануляю элементы
	for (size_t i = 0; i < 2; ++i){
        for(size_t j = 0; j < 2; ++j) {l[i][j] = 0.0; u[i][j] = 0.0;}
    }

	double* gf    = new double [2];
	double* delta = new double [2];
	double* ans   = new double [3];

	int num_it = 0;
	double eps = 1e-6;

	ans[0] = x0, ans[1] = y0;

	w[0][0] = f1_deriv_x(ans[0], ans[1]);
    w[0][1] = f1_deriv_y(ans[0], ans[1]);
    w[1][0] = f2_deriv_x(ans[0], ans[1]);
    w[1][1] = f2_deriv_y(ans[0], ans[1]);

    decompose(w, l, u, 2);

	for ( ; ; )
	{
		++num_it;
		if (num_it > 1000)
			break;

		gf[0] = f1(ans[0], ans[1]);
		gf[1] = f2(ans[0], ans[1]);
		delta = lu(l, u, gf, 2);
		ans[0] -= delta[0];
		ans[1] -= delta[1];
		if (l2(delta, 2) / l2(ans, 2) < eps) break;
	}

	ans[2] = num_it;
	for (size_t i = 0; i < 2; ++i) {delete [] l[i]; delete [] u[i]; delete [] w[i];}
	delete [] l;
	delete [] u;
	delete [] gf;
	delete [] delta;
	delete [] w;

	return ans;
}


double* lu(double** l, double** u, double* b, size_t n){

    double* x = new double [n];
	double* y = new double [n];

	for(size_t i = 0; i < n; ++i) {x[i] = 0.0; y[i] = 0.0;}

	double sum = 0.0;

	for (size_t i = 0; i < n; ++i)
	{
		sum = 0.0;

		for (size_t j = 0; j < i; ++j)
			sum += l[i][j] * y[j];
		y[i] = (b[i] - sum) / l[i][i];

	}

	for (int i = n - 1; i >= 0; --i)
	{
	    sum = 0.0;

		for (int j = n - 1; j >= i; --j)
			sum += u[i][j] * x[j];
		x[i] = (y[i] - sum) / u[i][i];
	}

	delete [] y;
	return x;
}


void decompose(double** a, double** l, double** u, size_t n){

    for(size_t i = 0; i < n; ++i){

        for(size_t j = i; j < n; ++j){
            l[j][i] = a[j][i];
            for(size_t k = 0; k < i; ++k){
               l[j][i] = l[j][i] - l[j][k] * u[k][i];
            }
         }

        for(size_t j = i; j < n; ++j) {
            if(j == i)
                u[i][j] = 1.0;
            else{
                u[i][j] = a[i][j] / l[i][i];
                for(size_t k = 0; k < i; ++k)
                    u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
            }
        }
    }
}


double f1_deriv_num(double x, double y, double h){
    return x ? (f1(x + h, y) - f1(x - h, y)) / (2*h) : (f1(x, y + h) - f1(x, y - h)) / (2*h);
}


double f2_deriv_num(double x, double y, double h){
    return x ? (f2(x + h, y) - f2(x - h, y)) / (2*h) : (f2(x, y + h) - f2(x, y - h)) / (2*h);
}


double* newton_disc(double x0, double y0){
	double** w = new double* [2];
	for (int i = 0; i < 2; ++i) {w[i] = new double [2];}

	double* gf    = new double [2];
	double* delta = new double [2];
	double* ans   = new double [3];

	int num_it = 0;
	double eps = 1e-3;

	ans[0] = x0, ans[1] = y0;
	for ( ; ; )
	{
		++num_it;
		if (num_it > 1000)
			break;

        double h = 1e-3;

		w[0][0] = f1_deriv_num(ans[0] + h, 0.0, h);
		w[0][1] = f1_deriv_num(0.0, ans[1], h);
		w[1][0] = f2_deriv_num(ans[0], 0.0, h);
		w[1][1] = f2_deriv_num(0.0, ans[1], h);
		gf[0] = f1(ans[0], ans[1]);
		gf[1] = f2(ans[0], ans[1]);
		delta = gauss(w, gf, 2);
		ans[0] -= delta[0];
		ans[1] -= delta[1];
		if (l2(delta, 2) / l2(ans, 2) < eps) break;
	}
	ans[2] = num_it;

	for (size_t i = 0; i < 2; ++i) delete [] w[i];
	delete [] gf;
	delete [] delta;
	delete [] w;

	return ans;
}


double* gauss(double** a, double* y, size_t n)
{
	double* x = new double [n];
	double	max;
	size_t index;
	double eps = 1e-5;  // точность
	for (size_t k = 0; k < n; ++k)
	{
	    // Поиск строки с максимальным a[i][k]
		max = abs(a[k][k]);
		index = k;
		for (size_t i = k + 1; i < n; ++i)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			std::cout << "error ";
			std::cout << index << std::endl;
			return x;
		}

		//swap(...)
		for (size_t j = 0; j < n; ++j)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double tmp = y[k];
		y[k] = y[index];
		y[index] = tmp;

		// Нормализация уравнений
        // обработка k-ой строчки
        tmp = a[k][k];
        for (size_t j = 0; j < n; ++j)
            a[k][j] = a[k][j] / tmp;
        y[k] = y[k] / tmp;

		for (size_t i = k + 1; i < n; ++i)
		{
			tmp = a[i][k];

			if (abs(tmp) < eps) continue; // для нулевого коэффициента пропустить

			for (size_t j = 0; j < n; ++j)
				a[i][j] = a[i][j] / tmp - a[k][j];

			y[i] = y[i] / tmp  - y[k];
		}
	}

	// обратная подстановка
	for (int k = n - 1; k >= 0; --k)
	{
		x[k] = y[k];
		for (int i = 0; i < k; ++i)
			y[i] = y[i] - a[i][k] * x[k];
	}

	return x;
}
