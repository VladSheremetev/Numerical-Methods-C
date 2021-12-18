#include "function_declaration.h"


int main(void)
{
    double  a; 			//координата x левой точки отрезка [a, b]
    double	b; 			//координата x правой точки отрезка [a, b]
    double	*res_x; 	//координаты x узлов на интервалах
    double	*res_y;		//координаты y узлов на интервалах
    double	*xint;		//координаты x равновномерной сетки на интервалах
    double	*y;			//координаты y равновномерной сетки на отрезке [a, b]
    double	*knots_x;	//координаты x узлов на отрезке [a, b]
    double	*knots_y;	//координаты y узлов на отрезке [a, b]
    double	*res_l;		//координаты y интерполяции Лагранжа на интервалах
	double 	u_g_s;		//шаг равномерной сетки M_viz
	double	*x;			//координаты x равновномерной сетки на отрезке [a, b]
	double	*l;			//координаты y интерполяции Лагранжа на отрезке [a, b]
	double	h;			//шаг между узлами интерполяции
	double	*x_er_new;	//координаты x равновномерной сетки из h/100 элементов на интервалах
	double	*l_er;		//координаты y интерполяции Лагранжа на равномерной сетке h/100 на интервалах
	double	*x_e;		//координаты x равновномерной сетки из h/100 элементов на отрезке [a, b]
	double	*l_e;		//координаты y интерполяции Лагранжа на равномерной сетке h/100 на отрезке [a, b]
	double	*y_e;		//координаты y на равномерной сетке h/100 на отрезке [a, b]
	size_t  M_viz;		//количество элементов в равномерной сетке
	size_t  K;			//количество интервалов
	size_t	M;			//количество узлов сетки
	size_t	N;			//степень многочлена Лагранжа на интервалах
	size_t	c_s_u;		//количество шагов равномерной сетки N_viz на интервале (от одного конца до другого)
	size_t	i_l = 0;	//индекс у массива с координатами интерполяции на сетке M_viz
	size_t	i_x = 0;	//индекс у массива с координатами равномерной сетки M_viz
	size_t	i_k_x = 0;	//индекс у массива с координатами узлов
	size_t	i_x_e = 0;	//индекс у массива с координатами равномерной сетки h/100
	size_t	i_l_e = 0;	//индекс у массива с координатами интерполяции на сетке h/100

    a = 0.0;
    b = 1.0;
    M_viz = 500;
    N = 2;
    K = 10;
    M = N * K;
    h = (b - a) / M;
    u_g_s   = (b - a) / M_viz;
    c_s_u   = count_step_between_point(a, a + (b - a)/K, u_g_s);
    knots_x = new double [M + 1];
    l       = new double [M_viz + 1];
    x       = new double [M_viz + 1];
    x_e     = new double[M*100 + 1];
    l_e     = new double[M*100 + 1];
    double lenght_int = (b - a) / K;//длина интервала
    size_t count_point_int = N + 1; //количество узлов на интервале
    double H = h / 100;				//шаг сетки h/100
    double cpih = N * 100 + 1;		//количество точек равномерной сетки h/100 на интервале
    size_t M_H = M*100;				//количество шагов на интервале при сетке h/100

    for (size_t i = 0; i < K; ++i)
    {
    	res_x    = x_gen(a + i * lenght_int, count_point_int, h);
        res_y    = y_gen(res_x, count_point_int);
        xint     = x_gen(a + i * lenght_int, c_s_u, u_g_s);
        x_er_new = x_gen(a + i * lenght_int, cpih, H);
        res_l    = l_gen(res_x, res_y, xint, c_s_u, count_point_int);
        l_er     = l_gen(res_x, res_y, x_er_new, cpih, count_point_int);

        for(size_t j = 0; j < c_s_u; ++j, ++i_x, ++i_l){
        	l[i_l] = res_l[j];
        	x[i_x] = xint[j];
        	if(j < N){
        		knots_x[i_k_x] = res_x[j];
        		++i_k_x;
        		}
        }

        for(size_t j = 0; j < cpih; ++j, ++i_x_e, ++i_l_e){
        	x_e[i_x_e] = x_er_new[j];
        	l_e[i_l_e] = l_er[j];
        }

    }

    y       = y_gen(x, M_viz);
    y_e     = y_gen(x_e, M_H);
    knots_y = y_gen(knots_x, M);
    error_calculation(y_e, l_e, M_H);

    print_file("C:\\Users\\Владислав\\lagrange.csv", x, l, M_viz);
    print_file("C:\\Users\\Владислав\\func.csv", x, y, M_viz);
    print_file("C:\\Users\\Владислав\\knots.csv", knots_x, knots_y, M);
    system("C:\\Users\\Владислав\\lagrange.py");

    delete [] res_x;
    delete [] res_y;
    delete [] xint;
    delete [] knots_x;
    delete [] knots_y;
    delete [] res_l;
    delete [] x;
    delete [] y;
    delete [] l;
    return 0;
}
