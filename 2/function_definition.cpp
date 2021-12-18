#include "function_declaration.h"

#define _USE_MATH_DEFINES


double*	x_gen(double a, size_t count_uniform_points, double step)
{
    double*	x = new double [count_uniform_points];

    for (size_t i = 0; i < count_uniform_points; ++i)
    {
    	x[i] = a + step * i;
    }

    return x;
}


double* y_gen(double* mass_x, size_t count_point)
{
    double*	y = new double [count_point];

    for (size_t index = 0; count_point > 0; ++index, --count_point)
    {
        y[index] = cos(2.0*M_PI*mass_x[index]);
    }

    return y;
}


double* l_gen(double* mass_x, double* mass_y, double* new_x, size_t count_new_x, size_t count_point)
{
    double*	l = new double[count_new_x];
    for (size_t i = 0; i < count_new_x; ++i)
        l[i] = lagrange(mass_x, mass_y, count_point, new_x[i]);

    return l;
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
    rel_norm_1 = abs_norm_1 / rel_norm_1;
    abs_norm_2 = sqrt(abs_norm_2);
    rel_norm_2 = abs_norm_2 / rel_norm_2;
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


void print_array(double* array, size_t size)
{
    for (size_t i = 0; i < size; ++i)
        std::cout << array[i] << std::endl;
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


double lagrange(double* mass_x, double*	mass_y, size_t count_point, double point)
{
    double  basic, lagrange = 0;

        for (size_t i = 0; i < count_point; ++i)
        {
            basic = 1;
            for (size_t j = 0; j < count_point; ++j)
            {
            	    if(j != i)
                    basic *= ((point - mass_x[j]) / (mass_x[i] - mass_x[j]));
            }
            lagrange += basic * mass_y[i];
        }
    return lagrange;
}


int count_step_between_point(double point_1, double point_2, double step){
	int count_step = 0;
	for(double eps = 0.0000000000000000001 ; point_1 - point_2 < eps ; ++count_step, point_1 += step);
	return count_step;
}
