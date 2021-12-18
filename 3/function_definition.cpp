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


double*   x_rand_gen(double a, double b, size_t count_point)
{
    double*		x = new double [count_point];
    double 		point;

    srand(static_cast<unsigned int>(time(NULL)));
    for (size_t i = 0 ; i < count_point; )
    {
    	point = a + (b - a) * (static_cast<double>(rand()) / (static_cast<double>(RAND_MAX)));
    	if (!search_point(x, i, point)){
    		x[i] = point;
    		++i;
    	}
    }
    quick_sort(x, 0, count_point - 1);
    return x;
}


bool search_point(double* array, size_t size_array, double point){
	double eps = 0.01;
	for (size_t i = 0; i < size_array; ++i){
		if (abs(point - array[i]) < eps) {return true;}
	}
	return false;
}


double* y_gen(double* x, size_t count_point)
{
    double*	y = new double [count_point];

    for (size_t i = 0; i < count_point; ++i)
    {
        y[i] = cos(2.0*M_PI*x[i]);
    }

    return y;
}


double      *lagrangebasic(double *x_arr, size_t N, double x, size_t p)
{
    double  *basicfun = new double [N];
    int ielem = p * (N - 1);

    	for (size_t i = 0; i < N; ++i)
        {
    		double basic = 1.0;
            for (size_t j = 0; j < N; ++j)
            {
                if (j != i)
                	basic *= ((x - x_arr[j + ielem]) / (x_arr[i + ielem] - x_arr[j + ielem]));
            }
            basicfun[i] = basic;
        }

    return basicfun;
}


double      *best_fit_func_rand(double *xrand, size_t K, double *koef, size_t N, size_t L, double *x)
{

    double  *res;
    double  *grand = new double [K *  L];

    for (size_t i = 0; i < K; ++i)
    {
        for (size_t j = 0; j < L; ++j)
        {
        	double lag = 0.0;
            if (x[i * (N - 1)] <= xrand[j + i * L] && xrand[j + i * L] <= x[(i + 1) * (N - 1)])
                {
                    res = lagrangebasic(x, N, xrand[j + i * L], i);
                    for (size_t k = 0; k < N; k++) lag += res[k] * koef[i * (N - 1) + k];
                    grand[i * L + j] = lag;
                }
        }
    }
    delete [] res;
    return grand;
}


double      *best_fit_func(double *xviz, size_t M_viz, double *x, size_t K, size_t N, double *koef)
{
    double  *g = new double [M_viz];
    double  *res;

    for (size_t i = 0; i < M_viz; ++i)
    {
    	double lag = 0.0;
        for (size_t j = 0; j < K; ++j)
        {
        	int ielem = j * (N - 1);

            if (x[ielem] <= xviz[i] && xviz[i] <= x[ielem + (N - 1)])
            {
                res = lagrangebasic(x, N, xviz[i], j);
                for (size_t k = 0; k < N; ++k) lag += res[k] * koef[ielem + k];
                g[i] = lag;
            }
        }
    }
    delete [] res;
    return g;
}


void  error_calculation(double* y, double* l, size_t size_array)
{
    double  abs_norm_1 = 0.0,
    		abs_norm_2 = 0.0,
			rel_norm_1 = 0.0,
			rel_norm_2 = 0.0,
			tmp_abs_norm_0 = 0.0,
			tmp_rel_norm_0 = 0.0,
    		abs_norm_0,
			rel_norm_0;

    for (size_t i = 0; i < size_array; ++i)
    {
        abs_norm_1 += abs(l[i] - y[i]);
        rel_norm_1 += abs(y[i]);
        abs_norm_2 += pow(abs(l[i] - y[i]), 2);
        rel_norm_2 += pow(abs(y[i]), 2);
        if(tmp_abs_norm_0 < abs(l[i] - y[i])) tmp_abs_norm_0 = abs(l[i] - y[i]);
        if(tmp_rel_norm_0 < abs(y[i])) tmp_rel_norm_0 = abs(y[i]);

    }
    rel_norm_1 = abs_norm_1 / rel_norm_1;
    abs_norm_2 = sqrt(abs_norm_2);
    rel_norm_2 = abs_norm_2 / rel_norm_2;
    abs_norm_0 = tmp_abs_norm_0;
    rel_norm_0 = tmp_rel_norm_0;
    rel_norm_0 = abs_norm_0 / rel_norm_0;
    std::cout<<"Absolute norme 0 := "<<std::scientific<<abs_norm_0<<std::endl;
    std::cout<<"Absolute norme 1 := "<<std::scientific<<abs_norm_1<<std::endl;
    std::cout<<"Absolute norme 2 := "<<std::scientific<<abs_norm_2<<std::endl;
    std::cout<<"Relative norme 0 := "<<std::scientific<<rel_norm_0<<std::endl;
    std::cout<<"Relative norme 1 := "<<std::scientific<<rel_norm_1<<std::endl;
    std::cout<<"Relative norme 2 := "<<std::scientific<<rel_norm_2<<std::endl;
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


double* slau(double *res_x, size_t N, size_t L, size_t K, double *xrandsave, double *yrandsave){
	double  *xrand; 					//случайные точки сгенерированные для каждого конечного элемента
	double  *yrand; 					//значение функции в случайных точках, сгенерированных для каждого конечного элемента
	size_t  M = (N - 1) * K + 1;
	double  **phi = new double* [L];	//базисные функции Лагранжа
	for (size_t i = 0; i < L; ++i) phi[i] = new double [N];
	double  *koef = new double [M];		//коэффициенты для функции наилучшего приближения
	Eigen::MatrixXf A(M, M);			//матрица скалярного произведения базисных функций
	Eigen::VectorXf B(M), C1;
	A.setZero();
	B.setZero();

	for (size_t p = 0; p < K; p++)
	{
				size_t ielem = p * (N - 1);

	    	    xrand = x_rand_gen(res_x[ielem], res_x[ielem + (N - 1)], L);
	    	    yrand = y_gen(xrand, L);

	    	    for (size_t i = 0; i < L; ++i)
	    	        {
	    	        yrandsave[p*L + i] = yrand[i];
	    	        xrandsave[p*L + i] = xrand[i];
	    	        phi[i] = lagrangebasic(res_x, N, xrand[i], p);
	    	        }


	    	    for (size_t i = 0; i < N; ++i)
	    	        {
	    	            for (size_t m = 0; m < N; ++m)
	    	            {
	    	                for (size_t k = 0; k < L; ++k)
	    	                    A(ielem + m, ielem + i) += phi[k][i] * phi[k][m];
	    	            }

	    	            for (size_t k = 0; k < L; ++k) B(ielem + i) += yrand[k] * phi[k][i];
	    	        }
	 }

	    //cout<<A<<endl;

	    //GAUSS
	    C1 = A.colPivHouseholderQr().solve(B);
	    //cout<<C1<<endl;

	    for(size_t i = 0; i < M; ++i) koef[i] = C1(i);

	    return koef;

}
