#include "interval_array.h"

IntervalArray::IntervalArray()
{
    leaf_ints = NULL;
    intervals = NULL;
}

IntervalArray::IntervalArray(ValueType *ints, size_t n) 
{
    build(ints, n);
}

IntervalArray::~IntervalArray()
{
    clear();
}

void IntervalArray::clear()
{
    if ( intervals != NULL )
        delete[] intervals;
    intervals = NULL;
}

void IntervalArray::build(ValueType *ints, size_t n)
{
    leaf_ints = ints;
    leaf = n;

    size = n-1;
    intervals = new ValueType[size];
    memset(intervals, MAX_INTERVAL_VALUE, size*sizeof(ValueType));

    fill( 1, n );
    intervals[0] = 0;
}

ValueType IntervalArray::fill( size_t l, size_t r )
{
    if ( l == r ) 
        return MAX_INTERVAL_VALUE;

    size_t m = getIndex(l,r);
    assert(m>0);

    if ( l == r-1 ) {
        return leaf_ints[l];
    }

    ValueType lv = fill( l, m );
    ValueType rv = fill( m, r );
    ValueType min = std::min( lv, rv );
    
    intervals[m-1] = min;
    return min;
}

void IntervalArray::print( std::ostream &out )
{
    for ( size_t i = 0; i < size; i++ )
        out << (i+1) << "\t" << (int)intervals[i] << "\n";
}

ValueType IntervalArray::getValue(size_t i, size_t j)
{
    if ( i >= j ) return 0;

    long k = getIndex(i, j);
    if ( k < 0 ) {
        std::cout << "Invalid index\n";
        return 0;
    }
    return intervals[k-1];
}

long IntervalArray::getIndex(size_t l, size_t r)
{
    if ( r < l ) return -1;
    return (l+r)/2;
}

void IntervalArray::dump( const char *filename )
{
    fio::write<ValueType>( intervals, size, filename );
}

void IntervalArray::load( const char *filename, size_t filesize )
{
    clear(); // drop array if exists.

    size = filesize;
    int n = fio::read<ValueType>( intervals, filename, 0 );
    assert( (size_t)n == size );
}
