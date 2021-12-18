#include "prak4.h"
#define _USE_MATH_DEFINES

double   *x_gen(double a, size_t count_uniform_points, double step)
{
    double*	x = new double [count_uniform_points];

    for (size_t i = 0; i < count_uniform_points; ++i)
    {
    	x[i] = a + step * i;
    }

    return x;
}

double   *oneh_derivative(double *y, size_t M, double h)
{
    double  *res = new double [M];

    res[0] = (y[1] - y[0]) / (h);
    res[M - 1] = (y[M - 1] - y[M - 2]) / (h);

    for (size_t i = 1 ; i < M - 1 ; ++i)
    {
        res[i] = ((y[i + 1]) - (y[i - 1])) / (2.0 * h);
    }

    return res;
}

double   *halfh_derivative(double *y, size_t M, double h)
{
	double  *res = new double [M];
	res[0] = (y[1] - y[0]) / (h);
	res[M - 1] = (y[M - 1] - y[M - 2]) / (h);

	for (size_t i = 1 ; i < M - 1 ; ++i)
	{
		res[i] = ((y[i + 1]) - (y[i - 1])) / (2.0 * h);
	}

    return res;
}

double   *y_function(double *x, size_t count_point)
{
	double*	y = new double [count_point];

	for (size_t i = 0; i < count_point; ++i)
	{
        y[i] = sin(2.0*M_PI*x[i]); //pow(mass_x[i], 2);
	}

    return y;
}


double   *y_derivative(double *x, size_t count_point)
{
	double*	y = new double [count_point];

	for (size_t i = 0; i < count_point; ++i)
	{
		y[i] = 2.0*M_PI*cos(2.0*M_PI*x[i]); //2.0 * mass_x[i];
	}

    return y;
}


double   *runge_calculation(double *one, double *half, size_t M)
{
    double  *runge = new double [M];;

    runge[0] = (half[0] - one[0]) + half[0];
    runge[M - 1] = (half[2 * (M - 1)] - one[M - 1]) + half[2 * (M - 1)];
    for (size_t i = 1 ; i < M - 1; ++i)
    {
        runge[i] = (half[2 * i] - one[i]) / 3.0 + half[2 * i];
    }
    return runge;
}


void print_file(std::string name_file, double* x, double* y, size_t size_array_points){
	std::ofstream out;
	out.open(name_file);
	out<<"x_coordinate,y_coordinate"<<std::endl;
	for(size_t i = 0; i < size_array_points; ++i)
	{
	   out<<x[i]<<",";
	   out<<y[i]<<std::endl;
	}
	out.close();
}


void  error_calculation(double* y, double* l, size_t size_array)
{
    double  abs_norm_1 = 0,
    		abs_norm_2 = 0,
			rel_norm_1 = 0,
			rel_norm_2 = 0,
			tmp_abs_norm_inf = 0,
			tmp_rel_norm_inf = 0,
    		abs_norm_inf,
			rel_norm_inf;

    for (size_t i = 0; i < size_array; ++i)
    {
        abs_norm_1 += abs(l[i] - y[i]);
        rel_norm_1 += abs(y[i]);
        abs_norm_2 += pow(abs(l[i] - y[i]), 2);
        rel_norm_2 += pow(abs(y[i]), 2);
        if(tmp_abs_norm_inf < abs(l[i] - y[i])) tmp_abs_norm_inf = abs(l[i] - y[i]);
        if(tmp_rel_norm_inf < abs(y[i])) tmp_rel_norm_inf = abs(l[i]);

    }
    rel_norm_1   = abs_norm_1 / rel_norm_1;
    abs_norm_2   = sqrt(abs_norm_2);
    rel_norm_2   = abs_norm_2 / rel_norm_2;
    abs_norm_inf = tmp_abs_norm_inf;
    rel_norm_inf = tmp_rel_norm_inf;
    rel_norm_inf = abs_norm_inf / rel_norm_inf;
    std::cout<<"Absolute norme 0 := "<<abs_norm_inf<<std::endl;
    std::cout<<"Absolute norme 1 := "<<abs_norm_1<<std::endl;
    std::cout<<"Absolute norme 2 := "<<abs_norm_2<<std::endl;
    std::cout<<"Relative norme 0 := "<<rel_norm_inf<<std::endl;
    std::cout<<"Relative norme 1 := "<<rel_norm_1<<std::endl;
    std::cout<<"Relative norme 2 := "<<rel_norm_2<<std::endl;
}

void print_array(double *array, size_t n)
{
    for (size_t i = 0; i < n; ++i)
        std::cout << array[i] << std::endl;
}

void quick_sort(double*	array, size_t low, size_t high)
{
	size_t i = low, j = high;
    double pivot = array[(i + j) / 2], tmp;

    for ( ; i <= j ; )
{
        for ( ; array[i] < pivot ; ++i);
        for ( ; array[j] > pivot ; --j);
        if (i <= j)
        {
        	tmp = array[i];
            array[i] = array[j];
            array[j] = tmp;
            ++i;
            --j;
        }
    }
    if (j > low) quick_sort(array, low, j);
    if (i < high) quick_sort(array, i, high);
}
