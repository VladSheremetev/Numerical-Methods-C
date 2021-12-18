#include "function_declaration.h"

#define N 5		//количество узлов
#define K 3		//количество конечных элементов
#define L 10	//количество внутренних (на конечном элементе) точек


int     main(void)
{
	size_t  M_viz;					//количество точек для визуализации
    double  a, b;					//координаты концов отрезка
    size_t  M = (N - 1) * K + 1;	//общее количество узлов
    double  *x;						//равномерная сетка с шагом h
    double  *y;						//значение функции на равномерной сетке
    double  *g;						//значения функции наилучшего приближения в точках
    double  *xviz;					//равномерная сетка из M_viz точек
    double  *yviz;					//значение на равномерной сетке из M_viz точек
    double  iter;					//шаг сетки для визуализации
    double  *koef;					//массив коэффициентов для функции наилучшего приближения
    double  *yrand;					//значения функции на случайных точках внутри конечных элементов
    double  *grand;					//значение функции наилучшего приближения в случайных точках внутри конечных элементов
    double  *xrand;					//координаты случайных точек внутри конечных элементов
    double  *xh;					//равномерная сетка с шагом h/100
    double  *yh;					//значение функции на равномерной сетке с шагом h/100
    double  *gh;					//значение функции наилучшего приближения на равномерной сетке с шагом h/100
    double  h;						//шаг равномерной сетки в случае M точек

    a = 0.0;
    b = 1.0;
    M_viz = 500; //(M - 1) * 100 + 1; //1500;
    iter = (b - a) / (M_viz - 1);
    h = (b - a) / (M - 1);

    yrand = new double [K * L];
    grand = new double [K * L];
    xrand = new double [K * L];
    //Равномерная сеткая с шагом h
    x = x_gen(a, M, h);
    y = y_gen(x, M);
    //Равномерная сетка с шагом h/100
    size_t M_h = (M - 1) * 100 + 1;
    double H   = h * 0.01;
    xh = x_gen(a, M_h, H);
    yh = y_gen(xh, M_h);
    //print_array(xh, M_h);

    //Равномерная сетка для визуализации
    xviz = x_gen(a, M_viz, iter);
    yviz = y_gen(xviz, M_viz);
    //print_array(xviz, M_viz);
    /*for (size_t i = 0; i < M_h; ++i)
            std::cout << yh[i] - yviz[i]<< std::endl;
    std::cout << std::endl;*/

    //Вычисление коэффициентов для функции наилучшего приближения
    koef = slau(x, N, L, K, xrand, yrand);
    //Вычисление значений функции наилучшего приближения
    grand = best_fit_func_rand  (xrand, K, koef, N, L, x);
    g 	  = best_fit_func		(xviz, M_viz, x, K, N, koef);
    gh    = best_fit_func		(xh,   M_h,   x, K, N, koef);

    /*for (size_t i = 0; i < M_h; ++i)
          std::cout << gh[i] - g[i]<< std::endl;*/

    std::cout<<"Error for rand points\n\n";
    error_calculation(yrand, grand, K * L);
    std::cout<<"\n\nError for h/100 \n\n";
    error_calculation(  yh, gh, M_h);
    std::cout<<"\n\nError for M_viz \n\n";
    error_calculation(yviz, g, M_viz);

    print_file("C:\\Users\\Владислав\\lagrange.csv", xviz,    g,  M_viz);
    print_file("C:\\Users\\Владислав\\func.csv", 	 xviz,  yviz, M_viz);
    print_file("C:\\Users\\Владислав\\knots.csv", 	xrand, yrand, K * L);
    system("C:\\Users\\Владислав\\lagrange.py");

    delete [] x;
    delete [] y;
    delete [] g;
    delete [] xviz;
    delete [] yviz;
    delete [] koef;
    delete [] grand;
    delete [] xrand;
    delete [] yrand;
    return 0;
}
