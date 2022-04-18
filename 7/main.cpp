#include "Header.h"

int main(){

    double* ans_s_it;
    double* ans_newton;
    double* ans_newton_mod;
    double* ans_newton_disc;

	//double x_0 = 0.72, y_0 = 0.89; //начальные точки
	double x_0 = 1.5, y_0 = 1.5;

	ans_s_it        = simple_it(x_0, y_0);
	ans_newton      = newton(x_0, y_0);
	ans_newton_mod  = newton_mod(x_0, y_0);
	ans_newton_disc = newton_disc(x_0, y_0);

	std::cout <<std::setw(10)<<"Simple itr: "<<std::setw(10)<<ans_s_it[2]<<std::setw(10)<<" Solution: "; print_array(ans_s_it, 2);
	std::cout <<std::setw(10)<<"Newton itr: " <<std::setw(10)<<ans_newton[2]<<std::setw(10)<<" Solution: "; print_array(ans_newton, 2);
	std::cout <<std::setw(10)<<"Newton modiv itr: "<<std::setw(4)<<ans_newton_mod[2]<<std::setw(10)<<" Solution: "; print_array(ans_newton_mod, 2);
	std::cout <<std::setw(10)<<"Newton discret itr: "<<std::setw(2)<<ans_newton_disc[2]<<std::setw(10)<<" Solution: "; print_array(ans_newton_disc, 2);

	system("python viz.py");
	
	delete [] ans_s_it;
	delete [] ans_newton;
	delete [] ans_newton_mod;
	delete [] ans_newton_disc;
	
	return 0;
}
