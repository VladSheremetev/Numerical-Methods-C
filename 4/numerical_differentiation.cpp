#include "prak4.h"

int     main(void)
{
    double  a, b; 	//координаты концов отрезка
    size_t  M;		//количество узлов сетки
    double  h;		//шаг сетки узлов
    double  *oneh;	//равномерная сетка с шагом h
    double  *halfh; //равномерная сетка с шагом h/2
    double  *yder;  //значения производной функции для сетки с h
    double  *yderh; //значения производной функции для сетки с h/2
    double  *runge; //массив с поправкой Рунге
    size_t  M_viz;	//количество элементов в равномерной сетке для визуализации
    double  iter;	//шаг сетки для визуализации
    double	*x, *y; //сетка для визуализации
    double  *onex, *oney, *halfx, *halfy; //генерация значений функции с узлами с шагом one == h или half == h/2

    a = 0.0;
    b = 1;
    M = 20;
    M_viz = 500;
    iter = (b - a) / (M_viz - 1);
    h = (b - a) / (M - 1);

    x = x_gen(a, M_viz, iter);
    y = y_derivative(x, M_viz);

    onex = x_gen(a, M, h);
    oney = y_function(onex, M);
    yder = y_derivative(onex, M);
    oneh = oneh_derivative(oney, M, h); //функция для численного дифференцирования методом центральных разностей

    double H   = h * 0.5;
    size_t M_H = 2.0 * M;
    halfx = x_gen(a, M_H, H);
    halfy = y_function(halfx, M_H);
    yderh = y_derivative(halfx, M_H);
    halfh = halfh_derivative(halfy, M_H, H); //функция для численного дифференцирования методом центральных разностей

    printf("error for one h\n\n");
    error_calculation(yder, oneh, M);
    printf("\nerror for half h\n\n");
    error_calculation(yderh, halfh, M_H);
    printf("\nrunge error\n\n");
    //Правило Рунге
    runge = runge_calculation(oneh, halfh, M);
    error_calculation(yder, runge, M);

    print_file("C:\\Users\\Владислав\\analytical_derivative.csv",   x,    y, M_viz);
    print_file("C:\\Users\\Владислав\\derivative_step_h.csv",    onex,   oneh, M);
    print_file("C:\\Users\\Владислав\\derivative_step_half.csv", halfx, halfh, M_H);
    print_file("C:\\Users\\Владислав\\derivative_runge.csv",     onex,  runge, M);
    system("C:\\Users\\Владислав\\derivative.py");

    delete [] yder;
    delete [] yderh;
    delete [] oneh;
    delete [] halfh;
    delete [] runge;
    delete [] x;
    delete [] y;
    delete [] onex;
    delete [] oney;
    delete [] halfx;
    delete [] halfy;

    return 0;
}
