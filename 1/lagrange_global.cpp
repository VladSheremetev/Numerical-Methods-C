#include "function_declaration.h"


int main(void)
{
	size_t  M;    // количество точек визуализации
    size_t 	N; 	  // количество узлов
    double  a, b; // координаты отрезка
    double  step, *x, *y, *xrand, *res_l, *yrand;

    a = 0;
    b = 1;
    M = 500;
    N = 40;
    step = (b - a) / (M - 1); //шаг для сетки визуализации

    xrand = x_rand_gen(a, b, N); //генератор сетки из случайных точек
    yrand = y_gen(xrand, N);

    x = x_gen(a, M, step); //генератор точек из равномерной сетки для построения интерполяции

    y = y_gen(x, M);

    res_l = l_gen(xrand, yrand, x, M, N); //функция для построения интерполирования методом Лагранжа

    print_file("C:\\Users\\Владислав\\lagrange.csv", x, res_l, M);
    print_file("C:\\Users\\Владислав\\func.csv", x, y, M);
    print_file("C:\\Users\\Владислав\\knots.csv", xrand, yrand, N);
    system("C:\\Users\\Владислав\\lagrange.py");

    delete [] x;
    delete [] y;
    delete [] res_l;
    delete [] xrand;
    delete [] yrand;

    return 0;
}
