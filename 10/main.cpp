#include "Header.h"


int main()
{
	double  a     = 0.0,b    = 1.0;
	double  t0    = 0.0,T    = 2.0;
	size_t  n_x0  = 10, n_t  = n_x0 * n_x0 * n_x0 * n_x0;
	size_t  n_x01 = 10, n_t2 = 1000;
	double  h_x   = (b - a)  / (n_x0 - 1);
	double  h_x1  = (b - a)  / (n_x01 - 1);
	double  h_t   = (T - t0) / (n_t - 1);
	double  h_t2  = (T - t0) / (n_t2 - 1);

    //равномерная сетка
	double* x0  = new double [n_x0];
	double* x01 = new double [n_x01];
	double* t   = new double [n_t];
	double* t2  = new double [n_t2];

	double sigma = 0.5;

	for (size_t i = 0; i < n_x0; ++i)
		x0[i] = a + i * h_x;
	for (size_t i = 0; i < n_x01; ++i)
		x01[i] = a + i * h_x1;
	for (size_t i = 0; i < n_t; ++i)
		t[i] = a + i * h_t;
	for (size_t i = 0; i < n_t2; ++i)
		t2[i] = a + i * h_t2;

    //численное значение функции в узлах сетки
	double **u_1_h, **u_2_h;

	u_1_h = explicit_schema(x0, t, n_x0, n_t, h_x, h_t);
	u_2_h = weighted_schema(x01, t2, n_x01, n_t2, sigma, h_x1, h_t);

    //сетка h/2
	size_t n_x1  = 2 * n_x0 - 1;
	size_t n_x11 = 2 * n_x01 - 1;
	h_x  /= 2.0;
	h_x1 /= 2.0;
	double* x1 = new double[n_x1];
	for (size_t i = 0; i < n_x1; ++i)
		x1[i] = a + i * h_x;

	double* x11 = new double[n_x11];
	for (size_t i = 0; i < n_x11; ++i)
		x11[i] = a + i * h_x1;

    //численное значение функции в узлах сетки h/2
	double **u_1_h2, **u_2_h2;
	u_1_h2 = explicit_schema(x1, t, n_x1, n_t, h_x, h_t);
	u_2_h2 = weighted_schema(x11, t2, n_x11, n_t2, sigma, h_x1, h_t);


	double  *err1, *err2;
	err1 = matrix_to_vector(grid_shift(u_1_h, u_1_h2, n_x0, n_t, n_x1), n_x0, n_t);
	err2 = matrix_to_vector(grid_shift(u_2_h, u_2_h2, n_x01, n_t2, n_x11), n_x01, n_t2);
	print_error(err1, err2, u_1_h, u_2_h, n_x0, n_t, n_x01, n_t2);

	std::ofstream fout("explicit.txt");
	for (size_t i = 0; i < n_x0; i++)
		fout << x0[i] << ' ';
	fout << std::endl;
	for (size_t i = 0; i < n_t; i++)
		fout << t[i] << ' ';
	fout << std::endl;
	for (size_t i = 0; i < n_x0; i++)
	{
		for (size_t j = 0; j < n_t; j++)
		{
			fout << u_1_h[i][j] << ' ';
		}
		fout << std::endl;
	}
	fout.close();

	std::ofstream fout1("weighted.txt");
	for (size_t i = 0; i < n_x01; i++)
		fout1 << x01[i] << ' ';
	fout1 << std::endl;
	for (size_t i = 0; i < n_t2; i++)
		fout1 << t2[i] << ' ';
	fout1 << std::endl;
	for (size_t i = 0; i < n_x01; i++)
	{
		for (size_t j = 0; j < n_t2; j++)
			fout1 << u_2_h[i][j] << ' ';
		fout1 << std::endl;
	}
	fout1.close();

	//визуализация
	system("python \"vis.py\" \"explicit.txt\"");
	system("python \"vis.py\" \"weighted.txt\"");

	delete [] x0;
	delete [] x01;
	delete [] t;
	delete [] t2;
	delete [] x1;
	delete [] x11;
	delete [] err1;
	delete [] err2;
	for (size_t i = 0; i < n_x0; ++i) delete [] u_1_h[i];
	delete [] u_1_h;
	for (size_t i = 0; i < n_x01; ++i) delete [] u_2_h[i];
	delete [] u_2_h;
	for (size_t i = 0; i < n_x1; ++i) delete [] u_1_h2[i];
	delete [] u_1_h2;
	for (size_t i = 0; i < n_x11; ++i) delete [] u_2_h2[i];
	delete [] u_2_h2;

	return 0;
}
