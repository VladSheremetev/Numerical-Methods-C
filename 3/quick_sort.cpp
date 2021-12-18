#include "function_declaration.h"


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
