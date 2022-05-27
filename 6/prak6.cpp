#include "prak6.h"
#define _USE_MATH_DEFINES


double*	x_gen(double a, size_t count_uniform_points, double step)
{
    double*	x = new double [count_uniform_points];

    for (size_t i = 0; i < count_uniform_points; ++i)
    {
    	x[i] = a + step * i;
    }

    return x;
}


double*   x_rand_gen(double a, double b, size_t count_point)
{
    double*		x = new double [count_point];
    double 		point;

    for (size_t i = 0 ; i < count_point; )
    {
    	point = a + (b - a) * (static_cast<double>(rand()) / (static_cast<double>(RAND_MAX)));
    	if (!search_point(x, i, point)) {
    		x[i] = point;
    		++i;
    	}
    }
    return x;
}


bool search_point(double* array, size_t size_array, double point){
	double eps = 0.001;
	for (size_t i = 0; i < size_array; ++i){
		if (abs(point - array[i]) < eps) return true;
	}
	return false;
}


double* y_gen(double* x, size_t count_point)
{
    double*	y = new double [count_point];

    for (size_t i = 0; i < count_point; ++i)
    {
        y[i] = pow(x[i], 3); //sin(2*M_PI*x[i]); //sin(x[i]);
    }

    return y;
}


double      *lagrangebasic(double *x_arr, size_t n, double x, size_t p)
{
    double  *basicfun = new double [n];
    int ielem = p * (n - 1);

    	for (size_t i = 0; i < n; ++i)
        {
    		double basic = 1.0;
            for (size_t j = 0; j < n; ++j)
            {
                if (j != i)
                	basic *= ((x - x_arr[j + ielem]) / (x_arr[i + ielem] - x_arr[j + ielem]));
            }
            basicfun[i] = basic;
        }

    return basicfun;
}


double scalar_multiplication(double* a, double* b, size_t n)		// скалярное умножение
{
	double res = 0.0;
	for (size_t i = 0; i < n; ++i) res += a[i] * b[i];
	return res;
}


// перемножение матриицы и вектора диагональных матриц
void matrix_multiplication(double* res, double** A, double* b, size_t n, int N)
{
	for (size_t i = 0; i < n; ++i){
        res[i] = 0.0;
        size_t k = i / (N - 1), s = i % (N - 1) != 0 || i == 0 ? (N - 1) * k : (N - 1) * (k - 1);
		for (size_t j = s; j < n; ++j)
        {
            if (j > s + (N - 1) && i % (N - 1) != 0) break;
            else if (j > s + 2*N - 2) break;
			res[i] += A[i][j] * b[j];
        }
	}
}


double* gauss(double** a, double* y, size_t n, size_t N)
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
		for (size_t i = k + 1; i - k < N && i < n; ++i)
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
		if (index != k)
        {
		for (size_t j = k; (j - k) < 2 * N && j < n; ++j)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double tmp = y[k];
		y[k] = y[index];
		y[index] = tmp;
        }
		// Нормализация уравнений
        // обработка k-ой строчки
		double tmp = a[k][k];
        for (size_t j = k; j - k < 2 * N && j < n; ++j)
            a[k][j] = a[k][j] / tmp;
        y[k] = y[k] / tmp;

		for (size_t i = k + 1; i - k < N && i < n; ++i)
		{
			tmp = a[i][k];

			if (abs(tmp) < eps) continue; // для нулевого коэффициента пропустить

			for (size_t j = k; j - k < 2 * N && j < n; ++j)
				a[i][j] = a[i][j] / tmp - a[k][j];

			y[i] = y[i] / tmp  - y[k];
		}
	}

	// обратная подстановка
	for (int k = n - 1; k >= 0; --k)
	{
		x[k] = y[k];
		for (int i = k - 1; k - i < 2*N && i >= 0; --i)
			y[i] = y[i] - a[i][k] * x[k];
	}

	return x;
}


double* lu(double** a, double* b, size_t n, int N){
	double** l = new double* [n];
	double** u = new double* [n];
	double* x = new double [n];
	double* y = new double [n];
	for (size_t i = 0; i < n; ++i){
        l[i] = new double [n];
        u[i] = new double [n];
        x[i] = 0.0;
        y[i] = 0.0;
	}

	for (size_t i = 0; i < n; ++i){
        for(size_t j = 0; j < n; ++j) {l[i][j] = 0.0; u[i][j] = 0.0;}
    }

    decompose(a, l, u, n, N);

	double sum = 0.0;

	for (size_t i = 0; i < n; ++i)
	{
		sum = 0.0;
		size_t k = i / (N - 1), s = i % (N - 1) != 0 || i == 0 ? (N - 1) * k : (N - 1) * (k - 1);

		for (size_t j = s; j < i; ++j){
            if (j > s + (N - 1) && i % (N - 1) != 0) break;
            else if (j > s + 2*N - 2) break;
			sum += l[i][j] * y[j];
		}
		y[i] = (b[i] - sum) / l[i][i];

	}

	for (int i = n - 1; i >= 0; --i)
	{
		sum = 0.0;

		size_t k = (n - 1 - i) / (N - 1), s = (n - 1 - i) % (N - 1) != 0 || (n - 1 - i) == 0 ? (N - 1) * k : (N - 1) * (k - 1);

		for (int j = n - 1 - s; j >= i; --j)
			sum += u[i][j] * x[j];
		x[i] = (y[i] - sum) / u[i][i];
	}

	for (size_t i = 0; i < n; ++i) {delete [] l[i]; delete [] u[i];}
	delete [] l;
	delete [] u;
	delete [] y;
	return x;
}

// разложение диагональной матрицы в нижнетреугольную и верхнетреугольную
void decompose(double** a, double** l, double** u, size_t n, int N){
    int i = 0, j = 0, k = 0;

    for (i = 0; i < n; ++i) {
        for (j = i; j - i < N && j < n; ++j) {
            l[j][i] = a[j][i];
            for (k = 0; k < i; ++k) {
               l[j][i] = l[j][i] - l[j][k] * u[k][i];
            }
         }

    for (j = i; j - i < N && j < n; ++j) {
        if (j == i)
            u[i][j] = 1;
        else {
            u[i][j] = a[i][j] / l[i][i];
        for (k = 0; k < i; ++k) {
            u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
            }
        }
      }
   }
}


double* relax(double** a, double* b, double w, size_t n, int N)
{
	double* x_new = new double [n];
	double* x_old = new double [n];
	double* diff  = new double [n];  //разность нового и старого
	double err = 1e3;  //погрешность
	double eps = 1e-7; //точность
	double sum = 0.0;
	size_t  num_of_iter = 0;
	for (size_t i = 0; i < n; ++i)
		x_old[i] = b[i] / a[i][i];
	while (err > eps){
		++num_of_iter;
		for (int i = 0; i < n; ++i)
		{
			sum = 0.0;
			x_new[i] = w * b[i] / a[i][i] - x_old[i] * (w - 1);
			for (int j = i - N <= 0 ? 0 : i - N;  j < i; ++j)
				sum += a[i][j] * x_new[j];
			for (int j = i+1; j - i < N && j < n; ++j)
				sum += a[i][j] * x_old[j];
			x_new[i] -= w * sum / a[i][i];
		}

		for (size_t i = 0; i < n; ++i)
			diff[i] = x_new[i] - x_old[i];
		err = l0(diff, n);
		for (size_t i = 0; i < n; ++i) x_old[i] = x_new[i];
	}
    delete [] x_new;
    delete [] diff;
	return x_old;
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


double* grad(double** a, double* b, size_t n, int N)
{
	size_t num_of_iter = 0;					// номер итерации
	double b_norm = l2(b, n);
	double* x = new double[n];
	for (size_t i = 0; i < n; ++i) x[i] = 0.0;
	double* r_old = new double[n];
	for (size_t i = 0; i < n; ++i) r_old[i] = b[i];
    double* r_new = new double[n];			// невязка
	double* p = new double[n];
	for (size_t i = 0; i < n; ++i) p[i] = r_old[i];
	double* p_new = new double[n];
	double* ap  = new double[n];						// вектор произведения
	double r_mult = scalar_multiplication(r_old, r_old, n), r_mult_tmp;
	double alpha, beta;							// константы для сопряженного градиента
	double eps = 1e-16;							// точность
	double err = 1000.0;						// погрешность
	while (err > eps){
		++num_of_iter;
		matrix_multiplication(ap, a, p, n, N);
		if (num_of_iter > 1000)
			break;
		alpha =  r_mult / scalar_multiplication(ap, p, n);
		for (size_t i = 0; i < n; ++i)
			x[i] += alpha * p[i];
		for (size_t i = 0; i < n; ++i)
			r_new[i] = r_old[i] - alpha * ap[i];
        r_mult_tmp = scalar_multiplication(r_new, r_new, n);
		beta = r_mult_tmp / r_mult;
		r_mult = r_mult_tmp;
		for (size_t i = 0; i < n; ++i)
			p_new[i] = r_new[i] + beta * p[i];
		for (size_t i = 0; i < n; ++i){r_old[i] = r_new[i]; p[i] = p_new[i];}
		err = l2(r_old, n) / b_norm;
	}

	delete[] r_old;
	delete[] r_new;
	delete[] p_new;
	delete[] ap;
	delete[] p;
	return x;
}


void  error_calculation(double* y, double* l, size_t size_array)
{
    double  abs_norm_1 = 0.0,
    		abs_norm_2 = 0.0,
			rel_norm_1 = 0.0,
			rel_norm_2 = 0.0,
			tmp_abs_norm_0 = 0.0,
			tmp_rel_norm_0 = 0.0,
    		abs_norm_0,
			rel_norm_0;

    for (size_t i = 0; i < size_array; ++i)
    {
        abs_norm_1 += abs(l[i] - y[i]);
        rel_norm_1 += abs(y[i]);
        abs_norm_2 += pow(abs(l[i] - y[i]), 2);
        rel_norm_2 += pow(abs(y[i]), 2);
        if(tmp_abs_norm_0 < abs(l[i] - y[i])) tmp_abs_norm_0 = abs(l[i] - y[i]);
        if(tmp_rel_norm_0 < abs(y[i])) tmp_rel_norm_0 = abs(y[i]);

    }
    rel_norm_1 = abs_norm_1 / rel_norm_1;
    abs_norm_2 = sqrt(abs_norm_2);
    rel_norm_2 = abs_norm_2 / rel_norm_2;
    abs_norm_0 = tmp_abs_norm_0;
    rel_norm_0 = tmp_rel_norm_0;
    rel_norm_0 = abs_norm_0 / rel_norm_0;
    std::cout<<"Absolute norme 0 := "<<std::scientific<<abs_norm_0<<std::endl;
    std::cout<<"Absolute norme 1 := "<<std::scientific<<abs_norm_1<<std::endl;
    std::cout<<"Absolute norme 2 := "<<std::scientific<<abs_norm_2<<std::endl;
    std::cout<<"Relative norme 0 := "<<std::scientific<<rel_norm_0<<std::endl;
    std::cout<<"Relative norme 1 := "<<std::scientific<<rel_norm_1<<std::endl;
    std::cout<<"Relative norme 2 := "<<std::scientific<<rel_norm_2<<std::endl;
}


void  error_nev(double* y, double* b, size_t size_array){
    std::cout<<"Absolute norme 0 := "<<std::scientific<<l0(y, size_array)<<std::endl;
    std::cout<<"Absolute norme 1 := "<<std::scientific<<l1(y, size_array)<<std::endl;
    std::cout<<"Absolute norme 2 := "<<std::scientific<<l2(y, size_array)<<std::endl;
    std::cout<<"Relative norme 0 := "<<std::scientific<<l0(y, size_array) / l0(b, size_array)<<std::endl;
    std::cout<<"Relative norme 1 := "<<std::scientific<<l1(y, size_array) / l1(b, size_array)<<std::endl;
    std::cout<<"Relative norme 2 := "<<std::scientific<<l2(y, size_array) / l2(b, size_array)<<std::endl;
}


int slau(double *res_x, size_t N, size_t L, size_t K){
	double  *xrand; 					//случайные точки сгенерированные для каждого конечного элемента
	double  *yrand; 					//значение функции в случайных точках, сгенерированных для каждого конечного элемента
	size_t  M = (N - 1) * K + 1;
	double  **phi = new double* [L];	//базисные функции Лагранжа
	double  *koef = new double [M];		//коэффициенты для функции наилучшего приближения
	Eigen::MatrixXf A(M, M);			//матрица скалярного произведения базисных функций
	Eigen::VectorXf B(M), C1;
	double  **a = new double* [M];
	double  *b = new double [M];
	double  *c1, *c2, *c3, *c4, *c5;
	for (size_t i = 0; i < M; ++i){a[i] = new double [M];}

	for (size_t i = 0; i < M; ++i){ b[i] = 0.0;
        for(size_t j = 0; j < M; ++j) a[i][j] = 0.0;
    }

	A.setZero();
	B.setZero();

	srand(static_cast<unsigned int>(time(NULL)));
	for (size_t p = 0; p < K; p++)
	{
				size_t ielem = p * (N - 1);
	    	    xrand = x_rand_gen(res_x[ielem], res_x[ielem + (N - 1)], L);
	    	    yrand = y_gen(xrand, L);

	    	    for (size_t i = 0; i < L; ++i)
	    	        {
	    	        phi[i] = lagrangebasic(res_x, N, xrand[i], p);
	    	        }


	    	    for (size_t i = 0; i < N; ++i)
	    	        {
	    	            for (size_t m = 0; m < N; ++m)
	    	            {
	    	                for (size_t k = 0; k < L; ++k){
	    	                    A(ielem + m, ielem + i) += phi[k][i] * phi[k][m];
	    	                    a[ielem + m][ielem + i] += phi[k][i] * phi[k][m];
	    	                }
	    	            }

	    	            for (size_t k = 0; k < L; ++k) {B(ielem + i) += yrand[k] * phi[k][i]; b[ielem + i] += yrand[k] * phi[k][i];}
	    	        }

	    	    for (size_t i = 0; i < L; ++i) delete [] phi[i];
	 }

		//std::cout<<A<<std::endl;

	    double*   tmp1 = new double [M];
	    double*   tmp2 = new double [M];
	    double*   tmp3 = new double [M];
	    double*   tmp4 = new double [M];
	    double*   tmp5 = new double [M];
        double** a_tmp = new double*[M];
        double*  b_tmp = new double [M];
        for (size_t i = 0; i < M; ++i) {a_tmp[i] = new double[M]; b_tmp[i] = b[i];}
        for (size_t i = 0; i < M; ++i) {for (size_t j = 0; j < M; ++j) a_tmp[i][j] = a[i][j];}

        C1 = A.colPivHouseholderQr().solve(B);
	    for(size_t i = 0; i < M; ++i) koef[i] = C1(i);
        matrix_multiplication(tmp1, a, koef, M, N);

	    c1 = gauss(a_tmp, b_tmp, M, N);
	    matrix_multiplication(tmp2, a, c1, M, N);
	    std::cout <<"Comparison of solutions for Gauss" << std::endl;
        error_calculation(koef, c1, M);
        std::cout <<std::endl;
        delete [] c1;

	    c2 = lu(a, b, M, N);
	    matrix_multiplication(tmp3, a, c2, M, N);
	    std::cout <<"Comparison of solutions for LU" << std::endl;
        error_calculation(koef, c2, M);
        std::cout <<std::endl;
        delete [] c2;

	    double w = 1.1;
        c3 = relax(a, b, w, M, N);
        matrix_multiplication(tmp4, a, c3, M, N);
        std::cout <<"Comparison of solutions for Relax upper" << std::endl;
        error_calculation(koef, c3, M);
        std::cout <<std::endl;
        delete [] c3;

        c4 = grad(a, b, M, N);
        matrix_multiplication(tmp5, a, c4, M, N);
        std::cout <<"Comparison of solutions for Grad" << std::endl;
        error_calculation(koef, c4, M);
        std::cout <<std::endl;
        delete [] c4;

        for (size_t i = 0; i < M; ++i)
        {
		tmp1[i] -= b[i];
		tmp2[i] -= b[i];
		tmp3[i] -= b[i];
		tmp4[i] -= b[i];
		tmp5[i] -= b[i];
        }

        std::cout <<"Comparing Residuals of Solutions for a Library Function" << std::endl;
        error_nev(tmp1, b, M);
        std::cout <<std::endl;
        std::cout <<"Comparing Residuals of Solutions for a Gauss" << std::endl;
        error_nev(tmp2, b, M);
        std::cout <<std::endl;
        std::cout <<"Comparing Residuals of Solutions for a LU" << std::endl;
        error_nev(tmp3, b, M);
        std::cout <<std::endl;
        std::cout <<"Comparing Residuals of Solutions for a Relax upper" << std::endl;
        error_nev(tmp4, b, M);
        std::cout <<std::endl;
        std::cout <<"Comparing Residuals of Solutions for a Grad" << std::endl;
        error_nev(tmp5, b, M);
        std::cout <<std::endl;

        for (size_t i = 0; i < M; ++i) {delete [] a[i];}
        for (size_t i = 0; i < M; ++i) {delete [] a_tmp[i];}
        delete [] a;
        delete [] tmp1;
        delete [] tmp2;
        delete [] tmp3;
        delete [] tmp4;
        delete [] tmp5;
        delete [] a_tmp;
        delete [] b;
        delete [] phi;
        delete [] xrand;
        delete [] yrand;
	    return 1;

}

//печать вектора (массива)
void print_array(double *arr, size_t n)
{
    for (size_t i = 0; i < n; ++i)
        std::cout << arr[i] << " ";
        std::cout<<std::endl;
}


// печать матрицы на экран
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
