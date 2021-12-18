#include "function_declaration.h"

#define _USE_MATH_DEFINES


double*   x_rand_gen(double a, double b, size_t count_point)
{
    double*		arrx = new double [count_point];
    int         count_point_tmp = count_point;
    double 		point;

    srand(static_cast<unsigned int>(time(NULL)));
    for (size_t index = 0 ; count_point > 0 ; )
    {
    	point = a + (b - a) * (static_cast<double>(rand()) / (static_cast<double>(RAND_MAX)));
    	if (!search_point(arrx, index, point)){
    		arrx[index] = point;
    		++index,
			--count_point;
    	}
    }
    quick_sort(arrx, 0, count_point_tmp);
    arrx[0] = a;
    arrx[count_point_tmp - 1] = b;
    return arrx;
}


bool search_point(double* array, size_t size_array, double point){
	double eps = 0.01;
	for (size_t i = 0; i < size_array; ++i){
		if (abs(point - array[i]) < eps) {return true;}
	}
	return false;
}


double*   y_gen(double *x, size_t count_point)
{
    double*	y = new double [count_point];

    for (size_t index = 0; count_point > 0; ++index, --count_point)
    {
    	y[index] = cos(2.0*M_PI*x[index]);
    }

    return y;
}


double  lagrange(double *knots_x, double *knots_y, size_t count_point, double point)
{
    double  basic, lagrange = 0;

        for (size_t i = 0; i < count_point; ++i)
        {
            basic = 1;
            for (size_t j = 0; j < count_point; ++j)
            {
                if (j != i)
                    basic *= ((point - knots_x[j]) / (knots_x[i] - knots_x[j]));
            }
            lagrange += basic * knots_y[i];
        }
    return lagrange;
}


double*  	x_gen(double a, size_t count_uniform_points, double step)
{
    double*	uniform_x = new double [count_uniform_points];

    for (size_t i = 0; i < count_uniform_points; ++i)
    {
    	uniform_x[i] = a + step * i;
    }

    return uniform_x;
}


double*   l_gen(double *knots_x, double *knots_y, double *uniform_points, size_t count_uniform_points, size_t count_point)
{
    double*	l = new double[count_uniform_points];

    for (size_t i = 0; i < count_uniform_points; ++i)
        l[i] = lagrange(knots_x, knots_y, count_point, uniform_points[i]);

    return l;
}


void print_array(double *array, size_t size)
{
    for (size_t i = 0; i < size; ++i) std::cout << array[i] << std::endl;
    std::cout <<std::endl;
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
