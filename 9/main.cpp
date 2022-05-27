#include "Header.h"

int main()
{
	double a = -1, b = 1.0;
	size_t n = 100;
	double q0 = 0.0, q1 = 0.95;
	double h = (b - a) / (n - 1);
	double* x;
	double* y;
	double* y_1;
	double* y_2;
	double* y_3;
	double* y1;
	x = zeros(n);
	y = zeros(n);
	y_1 = zeros(n);
	y_2 = zeros(n);
	y_3 = zeros(n);
	for (size_t i = 0; i < n; ++i) x[i] = a + i * h;
	y[0] = q0, y[n - 1] = q1;
	y_1[0] = q0, y_1[n - 1] = q1;
	y_2[0] = q0, y_2[n - 1] = q1;
	y_3[0] = q0, y_3[n - 1] = q1;
	ballistic(x, y, n);
	finite_subs(x, y_3, n);
	newton(x, y_1, n, q0, q1, a, b);
	chordes(x, y_2, n);


	std::ofstream fout("linear.txt");
	for (size_t i = 0; i < n; ++i)
		fout << x[i] << ' ';
	fout << std::endl;
	for (size_t i = 0; i < n; ++i)
		fout << y[i] << ' ';
	fout << std::endl;
	for (size_t i = 0; i < n; ++i)
		fout << y_3[i] << ' ';
	fout << std::endl;
	fout.close();
	std::ofstream fout1("nonlinear.txt");
	for (size_t i = 0; i < n; ++i)
		fout1 << x[i] << ' ';
	fout1 << std::endl;
	for (size_t i = 0; i < n; ++i)
		fout1 << y_1[i] << ' ';
	fout1 << std::endl;
	for (size_t i = 0; i < n; ++i)
		fout1 << y_2[i] << ' ';
	fout1 << std::endl;
	fout1.close();

	size_t n2 = 2 * n;
	double h2 = (b - a) / (n2 - 1);
	double* x2;
	double* y2;
	double* y2_1;
	double* y2_2;
	double* y2_3;
	double* y21;
	x2 = zeros(n2);
	y2 = zeros(n2);
	y2_1 = zeros(n2);
	y2_2 = zeros(n2);
	y2_3 = zeros(n2);
	for (int i = 0; i < n2; ++i) x2[i] = a + i * h2;

	y2[0] = q0, y2[n2 - 1] = q1;
	y2_1[0] = q0, y2_1[n2 - 1] = q1;
	y2_2[0] = q0, y2_2[n2 - 1] = q1;
	y2_3[0] = q0, y2_3[n2 - 1] = q1;
	ballistic(x2, y2, n2);
	newton(x2, y2_1, n2, q0, q1, a, b);
	chordes(x2, y2_2, n2);
	finite_subs(x2, y2_3, n2);

	std::cout<<std::endl;
	std::cout << std::setw(10)<<"Absolute error:" <<std::setw(12)<< "l2:" <<std::setw(24)<< "l1:" <<std::setw(16)<< "l0:" << std::endl;
	std::cout << "ballistic:\t\t" <<std::scientific<< l2(vector_shift(y, y2, n), n) << "\t\t" << l1(vector_shift(y, y2, n), n) << "\t" << l0(vector_shift(y, y2, n), n) << std::endl;
	std::cout << "finite:\t\t\t" << l2(vector_shift(y_3, y2_3, n), n) << "\t\t" << l1(vector_shift(y_3, y2_3, n), n)  << "\t" << l0(vector_shift(y_3, y2_3, n), n) << std::endl;
	std::cout << "newton:\t\t\t" << l2(vector_shift(y_1, y2_1, n), n) << "\t\t"<< l1(vector_shift(y_1, y2_1, n), n) << "\t" << l0(vector_shift(y_1, y2_1, n), n) << std::endl;
	std::cout << "chordes:\t\t" << l2(vector_shift(y_2, y2_2, n), n) << "\t\t" <<l1(vector_shift(y_2, y2_2, n), n) << "\t" << l0(vector_shift(y_2, y2_2, n), n)  << std::endl;
    std::cout << std::setw(10)<<"Otnosit error:" <<std::setw(13)<< "l2:" <<std::setw(24)<< "l1:" <<std::setw(16)<< "l0:" << std::endl;
    std::cout << "ballistic:\t\t" <<std::scientific<< l2(vector_shift(y, y2, n), n) / l2(y, n) << "\t\t" << l1(vector_shift(y, y2, n), n) / l1(y, n) << "\t" << l0(vector_shift(y, y2, n), n) / l0(y, n) << std::endl;
	std::cout << "finite:\t\t\t" << l2(vector_shift(y_3, y2_3, n), n) / l2(y_3, n)<< "\t\t" << l1(vector_shift(y_3, y2_3, n), n) / l1(y_3, n) << "\t" << l0(vector_shift(y_3, y2_3, n), n) / l0(y_3, n) << std::endl;
	std::cout << "newton:\t\t\t" << l2(vector_shift(y_1, y2_1, n), n) / l2(y_1, n) << "\t\t"<< l1(vector_shift(y_1, y2_1, n), n) / l1(y_3, n) << "\t" << l0(vector_shift(y_1, y2_1, n), n) / l0(y_1, n) << std::endl;
	std::cout << "chordes:\t\t" << l2(vector_shift(y_2, y2_2, n), n) / l2(y_2, n) << "\t\t" <<l1(vector_shift(y_2, y2_2, n), n) / l1(y_2, n) << "\t" << l0(vector_shift(y_2, y2_2, n), n) / l0(y_2, n)  << std::endl;

	system("python \"vis.py\" ");
	system("python \"vis2.py\" ");

	delete [] x;
	delete [] y;
	delete [] y_1;
	delete [] y_2;
	delete [] y_3;
	delete [] y1;
	delete [] x2;
	delete [] y2;
	delete [] y2_1;
	delete [] y2_2;
	delete [] y2_3;
	delete [] y21;
	return 0;
}
