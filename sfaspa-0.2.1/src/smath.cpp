#include "smath.h"


double math::sum(const double *data, size_t n) 
{
    double init = 0.0;
    return std::accumulate(data, data+n, init);
}

double math::mean(const double *data, size_t n) 
{
    return math::sum(data,n)/n;
}

double math::max(const double *data, size_t n)  
{
    return *std::max_element(data, data+n);
}
	
double math::min(const double *data, size_t n)  
{
    return *std::min_element(data, data+n);
}

double math::median(double *data, size_t n, bool sorted)  
{
    if (!sorted) math::sort(data, n);
    return *(data+n/2);
}

double math::quantile(double *data, size_t n, double f, bool sorted)  
{
    assert( f>=0.0 && f<=1.0 );

    if (!sorted) std::sort(data, data+n);

    int k = f==1.0 ? n-1 : int(n*f)-1;
    return data[k];
}

void math::sort(double *data, size_t n)  
{
    std::sort(data, data+n);
}

