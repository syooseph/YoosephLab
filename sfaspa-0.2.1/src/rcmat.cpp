#include "rcmat.h"

RangeCountMatrix::RangeCountMatrix()
{
}

RangeCountMatrix::RangeCountMatrix(int r, int c)
{
    init(r, c);
}


RangeCountMatrix::RangeCountMatrix(int len, int min, int max)
{
    init(len, min, max);
}

RangeCountMatrix::~RangeCountMatrix()
{
    clear();
}

void RangeCountMatrix::clear()
{
    if ( matrix.size() ) matrix.clear();
}

void RangeCountMatrix::init(int len, int min, int max) 
{
    init(len-min+1, max-min+1);
}

void RangeCountMatrix::init(int r, int c)
{
    nrow = r;
    ncol = c;

    matrix = Matrix(nrow, std::vector<int>(ncol,0));
}

void RangeCountMatrix::increment(int pos, int len, int min, int max)
{
    assert(len>=min);
    assert(pos<nrow);
    assert(pos+len<nrow+min);
    int iter = len-min+1;
    for ( int i = 0; i < iter; i++ ) {
        if ( len < min ) break;

        int cnum = len > max ? max-min+1 : len-min+1;
        increment(pos, 0, cnum-1);

        len--;
        pos++;
    }
}

void RangeCountMatrix::increment(int row, int beg, int end )
{
    for ( int i = beg; i <= end; i++ )
        matrix[row][i]++;
}

int RangeCountMatrix::get(int left, int right, int min, int max)
{
    int diff = right-left+1;
    assert(diff>=min && diff<=max);

    int row = left;
    int col = right-(left+min)+1;
    return get(row,col);
}

int RangeCountMatrix::get(int row, int col)
{
    assert(row<nrow);
    assert(col<ncol);
    return matrix[row][col];
}

void RangeCountMatrix::print()
{
    printf("%10s\t", "");
    for ( int i = 0; i < ncol; i++ ) printf("%4d\t",i);
    printf("\n");
    
    for ( int i = 0; i < nrow; i++ ) {
        printf("%10d\t", i);
        for ( int j = 0; j < ncol; j++ ) printf("%4d\t", matrix[i][j]);
        printf("\n");
    }
}

void RangeCountMatrix::trim(int s, int e)
{
    if ( s == 0 && e == nrow-1 ) return;

    if ( e < nrow-1 ) matrix = Matrix( matrix.begin(), matrix.begin()+(e+1) );
    if ( s > 0 ) matrix = Matrix( matrix.begin()+s, matrix.end() );

    nrow = matrix.size();

    if ( Param::verbose ) std::cout << "Matrix trimmed - nrow:" << nrow << "\n";
}
