#include "ReadStartCount.h"

ReadStartCount::ReadStartCount()
{
    size = 0;
    cumulatives = NULL;
}

ReadStartCount::ReadStartCount( const GSA *gsa, const unsigned &n )
{
    init(gsa,n);
}

void ReadStartCount::init( const GSA *gsa, const unsigned &n )
{
    assert(n>0);
    size = n;

    unsigned *starts = new unsigned[size];
    memset( starts, 0, sizeof(unsigned) * size );

    for ( unsigned i = 0; i < size; i++ ) {
        GsaType item = gsa->getAt(i);
        //-----------
        // Read start
        //-----------
        if ( item.pos == 0 ) starts[i]++;
    }

    cumulatives = new unsigned[size];
    memcpy( cumulatives, starts, size*sizeof(unsigned));

    for ( unsigned i = 1; i < size; i++ ) 
        cumulatives[i] += cumulatives[i-1];
    
    delete[] starts;
}

ReadStartCount::~ReadStartCount()
{
    if ( cumulatives != NULL ) delete[] cumulatives;
    cumulatives = NULL;
}

unsigned ReadStartCount::getCount(int a, int b)
{
    assert(b>=a && b>=0 && b<(int)size);

    unsigned count = a > 0 ? 
        cumulatives[b] - cumulatives[a-1] : 
        cumulatives[b];
    return count;
}

unsigned ReadStartCount::get(int a)
{
    assert(a>=0 && a<(int)size);

    return cumulatives[a];
}

