#include "Header.h"

int main()
{
	double a = 0.0, b = 1.0;
	size_t n = 200;
	double q0 = 1.0, q1 = 0.0;
	double h = (b - a) / (n - 1);

	double* x = new double [n];
	for (size_t i = 0; i < n; ++i) x[i] = a + i * h;
	double* y = new double [n];
	double*y_1 = new double [n];
	double*y_2 = new double [n];
	double*y_3 = new double [n];
	double*y_4 = new double [n];
	double*y1 = new double [n];
	double*y1_1 = new double [n];
	double*y1_2 = new double [n];
	double*y1_3 = new double [n];
	double*y1_4 = new double [n];

	y[0] = q0, y1[0] = q1;
	y_1[0] = q0, y1_1[0] = q1;
	y_2[0] = q0, y1_2[0] = q1;
	y_3[0] = q0, y1_3[0] = q1;
	y_4[0] = q0, y1_4[0] = q1;

	euler(x, y, y1, n);
	euler_predict(x, y_1, y1_1, n);
	runge_kutt_2(x, y_2, y1_2, n);
	runge_kutt_4(x, y_3, y1_3, n);
	adams_3(x, y_4, y1_4, n);

	size_t n2 = 2 * n;
	double h2 = (b - a) / (n2 - 1);
	double* x2 = new double [n2];
	for (size_t i = 0; i < n2; ++i) x2[i] = a + i * h2;
	double* y2 = new double[n2];
	double*y2_1 = new double[n2];
	double*y2_2 = new double[n2];
	double*y2_3 = new double[n2];
	double*y2_4 = new double[n2];
	double* y21 = new double[n2];
	double*y21_1 = new double[n2];
    double*y21_2 = new double[n2];
    double*y21_3 = new double[n2];
    double*y21_4 = new double[n2];
	y2[0] = q0, y21[0] = q1;
	y2_1[0] = q0, y21_1[0] = q1;
	y2_2[0] = q0, y21_2[0] = q1;
	y2_3[0] = q0, y21_3[0] = q1;
	y2_4[0] = q0, y21_4[0] = q1;

	euler(x2, y2, y21, n2);
	euler_predict(x2, y2_1, y21_1, n2);
	runge_kutt_2(x2, y2_2, y21_2, n2);
	runge_kutt_4(x2, y2_3, y21_3, n2);
	adams_3(x2, y2_4, y21_4, n2);


	size_t n3 = n * 0.5;
	double h3 = (b - a) / (n3 - 1);
	double* x3 = new double[n3];
	for (size_t i = 0; i < n3; ++i) x3[i] = a + i * h3;
	double* y3 = new double[n3];
	double*y3_1 = new double[n3];
	double*y3_2 = new double[n3];
	double*y3_3 = new double[n3];
	double*y3_4 = new double[n3];
	double* y31 = new double[n3];
    double*y31_1 = new double[n3];
    double*y31_2 = new double[n3];
    double*y31_3 = new double[n3];
    double*y31_4 = new double[n3];
	y3[0] = q0, y31[0] = q1;
	y3_1[0] = q0, y31_1[0] = q1;
	y3_2[0] = q0, y31_2[0] = q1;
	y3_3[0] = q0, y31_3[0] = q1;
	y3_4[0] = q0, y31_4[0] = q1;

	euler(x3, y3, y31, n3);
	euler_predict(x3, y3_1, y31_1, n3);
	runge_kutt_2(x3, y3_2, y31_2, n3);
	runge_kutt_4(x3, y3_3, y31_3, n3);
	adams_3(x3, y3_4, y31_4, n3);

	cout << std::setw(10)<<"Absolute error:" <<std::setw(16)<< "l2:" <<std::setw(32)<< "l1:" <<std::setw(32)<< "l0:" << endl;
	cout << "euler:\t\t" <<scientific<< l2(sub(y, y2, n, n2), n) << ", " << l2(sub(y, y3, n, n3), n3) << '\t' << l1(sub(y, y2, n, n2), n) << ", " << l1(sub(y, y3, n, n3), n3) << '\t' << l0(sub(y, y2, n, n2), n) << ", " << l0(sub(y, y3, n, n3), n3) << endl;
	cout << "euler_predict:\t" << l2(sub(y_1, y2_1, n, n2), n) << ", " << l2(sub(y_1, y3_1, n, n3), n3) << '\t' << l1(sub(y_1, y2_1, n, n2), n) << ", " << l1(sub(y_1, y3_1, n, n3), n3) << '\t' << l0(sub(y_1, y2_1, n, n2), n) << ", " << l0(sub(y_1, y3_1, n, n3), n3) << endl;
	cout << "runge_kutt_2:\t" << l2(sub(y_2, y2_2, n, n2), n) << ", " << l2(sub(y_2, y3_2, n, n3), n3) << '\t' << l1(sub(y_2, y2_2, n, n2), n) << ", " << l1(sub(y_2, y3_2, n, n3), n3) << '\t' << l0(sub(y_2, y2_2, n, n2), n) << ", " << l0(sub(y_2, y3_2, n, n3), n3) << endl;
	cout << "runge_kutt_4:\t" << l2(sub(y_3, y2_3, n, n2), n) << ", " << l2(sub(y_3, y3_3, n, n3), n3) << '\t' << l1(sub(y_3, y2_3, n, n2), n) << ", " << l1(sub(y_3, y3_3, n, n3), n3) << '\t' << l0(sub(y_3, y2_3, n, n2), n) << ", " << l0(sub(y_3, y3_3, n, n3), n3) << endl;
	cout << "adams_3:\t" << l2(sub(y_4, y2_4, n, n2), n) << ", " << l2(sub(y_4, y3_4, n, n3), n3) << '\t' << l1(sub(y_4, y2_4, n, n2), n) << ", " << l1(sub(y_4, y3_4, n, n3), n3) << '\t' << l0(sub(y_4, y2_4, n, n2), n) << ", " << l0(sub(y_4, y3_4, n, n3), n3) << endl;
	cout << std::setw(10)<<"Otnosit error:" <<std::setw(17)<< "l2:" <<std::setw(32)<< "l1:" <<std::setw(32)<< "l0:" << endl;
	cout << "euler:\t\t" <<scientific<< l2(sub(y, y2, n, n2), n) / l2(y, n) << ", " << l2(sub(y, y3, n, n3), n3) / l2(y, n3) << '\t' << l1(sub(y, y2, n, n2), n) / l1(y, n) << ", " << l1(sub(y, y3, n, n3), n3) / l1(y, n3) << '\t' << l0(sub(y, y2, n, n2), n) / l0(y, n) << ", " << l0(sub(y, y3, n, n3), n3) / l0(y, n3) << endl;
	cout << "euler_predict:\t" << l2(sub(y_1, y2_1, n, n2), n) / l2(y_1, n) << ", " << l2(sub(y_1, y3_1, n, n3), n3) / l2(y_1, n3) << '\t' << l1(sub(y_1, y2_1, n, n2), n) / l1(y_1, n) << ", " << l1(sub(y_1, y3_1, n, n3), n3) / l1(y_1, n3) << '\t' << l0(sub(y_1, y2_1, n, n2), n) / l0(y_1, n) << ", " << l0(sub(y_1, y3_1, n, n3), n3) / l0(y_1, n3) << endl;
	cout << "runge_kutt_2:\t" << l2(sub(y_2, y2_2, n, n2), n) / l2(y_2, n) << ", " << l2(sub(y_2, y3_2, n, n3), n3) / l2(y_2, n3) << '\t' << l1(sub(y_2, y2_2, n, n2), n) / l1(y_2, n) << ", " << l1(sub(y_2, y3_2, n, n3), n3) / l1(y_2, n3) << '\t' << l0(sub(y_2, y2_2, n, n2), n) / l0(y_2, n) << ", " << l0(sub(y_2, y3_2, n, n3), n3) / l0(y_2, n3) << endl;
	cout << "runge_kutt_4:\t" << l2(sub(y_3, y2_3, n, n2), n) / l2(y_3, n) << ", " << l2(sub(y_3, y3_3, n, n3), n3) / l2(y_3, n3) << '\t' << l1(sub(y_3, y2_3, n, n2), n) / l1(y_3, n) << ", " << l1(sub(y_3, y3_3, n, n3), n3) / l1(y_3, n3) << '\t' << l0(sub(y_3, y2_3, n, n2), n) / l0(y_3, n) << ", " << l0(sub(y_3, y3_3, n, n3), n3) / l0(y_3, n3) << endl;
	cout << "adams_3:\t" << l2(sub(y_4, y2_4, n, n2), n) / l2(y_4, n) << ", " << l2(sub(y_4, y3_4, n, n3), n3) / l2(y_4, n3) << '\t' << l1(sub(y_4, y2_4, n, n2), n) / l1(y_4, n) << ", " << l1(sub(y_4, y3_4, n, n3), n3) / l1(y_4, n3) << '\t' << l0(sub(y_4, y2_4, n, n2), n) / l0(y_4, n) << ", " << l0(sub(y_4, y3_4, n, n3), n3) / l0(y_4, n3) << endl;

	ofstream fout("y.txt");
	for (size_t i = 0; i < n; i++)
		fout << x[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < n; i++)
		fout << y[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < n; i++)
		fout << y_1[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < n; i++)
		fout << y_2[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < n; i++)
		fout << y_3[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < n; i++)
		fout << y_3[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < n; i++)
		fout << y_4[i] << ' ';
	fout << endl;

	fout.close();
	ofstream fout1("y1.txt");
	for (size_t i = 0; i < n; i++)
		fout1 << x[i] << ' ';
	fout1 << endl;
	for (size_t i = 0; i < n; i++)
		fout1 << y1[i] << ' ';
	fout1 << endl;
	for (size_t i = 0; i < n; i++)
		fout1 << y1_1[i] << ' ';
	fout1 << endl;
	for (size_t i = 0; i < n; i++)
		fout1 << y1_2[i] << ' ';
	fout1 << endl;
	for (size_t i = 0; i < n; i++)
		fout1 << y1_3[i] << ' ';
	fout1 << endl;
	for (size_t i = 0; i < n; i++)
		fout1 << y1_4[i] << ' ';
	fout1 << endl;
	fout1.close();
	system("python \"vis.py\" ");
	system("python \"vis2.py\" ");

	delete [] y;
	delete [] y_1;
	delete [] y_2;
	delete [] y_3;
	delete [] y_4;
	delete [] y1;
	delete [] y1_1;
	delete [] y1_2;
	delete [] y1_3;
	delete [] y1_4;
	delete [] x;
	delete [] y2;
	delete [] y2_1;
	delete [] y2_2;
	delete [] y2_3;
	delete [] y2_4;
	delete [] y21;
	delete [] y21_1;
	delete [] y21_2;
	delete [] y21_3;
	delete [] y21_4;
	delete [] x2;
	delete [] y3;
	delete [] y3_1;
	delete [] y3_2;
	delete [] y3_3;
	delete [] y3_4;
	delete [] y31;
	delete [] y31_1;
	delete [] y31_2;
	delete [] y31_3;
	delete [] y31_4;
	delete [] x3;
	return 0;
}
