#include "function_declaration.h"

#define _USE_MATH_DEFINES


double   *x_gen(double a, size_t count_point, double h)
{
	 double*   x = new double [count_point];

	    for (size_t i = 0; i < count_point; ++i)
	    {
	    	x[i] = a + h * i;
	    }

	 return x;
}


double  fun(double x)
{
    return pow(x, 1);
}


double  integral(double x)
{
    return pow(x, 2) / 2.0;
}


double   yintegral(double a, double b)
{
    return integral(b) - integral(a);
}


double  square_formula(double *x, size_t K, double h)
{
    double  res = 0;

    for (size_t i = 1; i <= K; ++i)
    {
    	res += fun((x[i] - h / 2.0)) * h;
    }

    return res;
}


double  trapeze_formula(double *x, size_t K, double h)
{
    double  res = 0;

    for (size_t i = 1; i <= K; ++i)
    {
    	res += (fun(x[i]) + fun(x[i - 1])) / 2.0 * h;
    }

    return res;
}


double  simpson_formula(double *x, size_t K, double h)
{
    double  res = 0;

    for (size_t i = 1; i <= K; ++i)
    {
    	res += (fun(x[i - 1]) + fun(x[i]) + 4.0 * fun((x[i-1] + x[i]) * 0.5)) / 6.0 * h;
    }

    return res;
}


double  newton_formula(double *x, size_t K, double h)
{
    double  res = 0, c[6] = {19.0, 75.0, 50.0, 50.0, 75.0, 19.0};


    for (size_t j = 0; j < K; ++j)
    {
        for (size_t i = 0; i < 6; ++i)
        {
            res += c[i] * fun(x[j] + i * h / 5.0);
        }
    }

    res = res / 288.0 * h;

    return res;
}


double  gauss_formula(double *x, size_t K, double h)
{
    double  res = 0;

    for (size_t i = 0; i < K; ++i)
    {
    	res += 5.0 * fun(x[i] + (x[i + 1] - x[i]) * (1.0 - sqrt(3.0 / 5.0)) / 2.0) + 5.0 * fun(x[i] + (x[i + 1] - x[i]) * (1.0 + sqrt(3.0 / 5.0)) / 2.0) + 8.0 * fun(x[i] + h / 2.0);
    }
    res = res * h / 2.0 / 9.0;

    return res;
}
