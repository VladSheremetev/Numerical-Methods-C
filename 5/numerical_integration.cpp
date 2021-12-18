#include "function_declaration.h"

int     main(void)
{
    double  a, b;		//координаты концов отрезка
    double	K;			//количество интервалов для разбиения отрезка [a, b]
    double  h;			//шаг равномерной сетки
    double  *x;			//равномерная сетка
    double  yint;		//значение интеграла на отрезке [a, b]

    double  square1;
    double  trapeze1;
    double  simpson1;
    double  newton1;
    double  gauss1;

    double  square2;
    double  trapeze2;
    double  simpson2;
    double  newton2;
    double  gauss2;

    double  err_square1;
    double  err_square2;
    double  err_trapeze1;
    double  err_trapeze2;
    double  err_simpson1;
    double  err_simpson2;
    double  err_newton1;
    double  err_newton2;
    double  err_gauss1;
    double  err_gauss2;

    a = 0.0;
    b = 1.0;
    K = 10;
    h = (b - a) / K;

    yint = yintegral(a, b);

    x = x_gen(a, K + 1, h);

    square1 = square_formula (x, K, h);
    trapeze1= trapeze_formula(x, K, h);
    simpson1= simpson_formula(x, K, h);
    newton1 = newton_formula (x, K, h);
    gauss1  = gauss_formula  (x, K, h);

    err_square1  = abs((yint - square1)) / abs(yint);
    err_trapeze1 = abs((yint - trapeze1))/ abs(yint);
    err_simpson1 = abs((yint - simpson1))/ abs(yint);
    err_newton1  = abs((yint - newton1)) / abs(yint);
    err_gauss1 	 = abs((yint - gauss1))  / abs(yint);

    h = h * 0.5;
    K = K * 2;

    x = x_gen(a, K + 1, h);

    square2  = square_formula (x, K, h);
    trapeze2 = trapeze_formula(x, K, h);
    simpson2 = simpson_formula(x, K, h);
    newton2  = newton_formula (x, K, h);
    gauss2   = gauss_formula  (x, K, h);

    err_square2 = abs((yint - square2))  / abs(yint);
    err_trapeze2= abs((yint - trapeze2)) / abs(yint);
    err_simpson2= abs((yint - simpson2)) / abs(yint);
    err_newton2 = abs((yint - newton2))  / abs(yint);
    err_gauss2  = abs((yint - gauss2))   / abs(yint);


    printf("    h    |   h/2 \n");
    printf("---------------------\n");
    printf("  %.2lf  |  %.2lf      mesh size\n", K / 2, K);
    printf("---------------------\n");
    printf("%.2e | %.2e    square\n",  err_square1,  err_square2);
    printf("---------------------\n");
    printf("%.2e | %.2e    trapeze\n", err_trapeze1, err_trapeze2);
    printf("---------------------\n");
    printf("%.2e | %.2e    simpson\n", err_simpson1, err_simpson2);
    printf("---------------------\n");
    printf("%.2e | %.2e    newton\n",  err_newton1,  err_newton2);
    printf("---------------------\n");
    printf("%.2e | %.2e    gauss\n",   err_gauss1,   err_gauss2);


    delete [] x;

    return 0;
}
