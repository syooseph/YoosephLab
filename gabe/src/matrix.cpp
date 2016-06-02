#include "matrix.h"

using namespace std;

Matrix::Matrix()
{
}

Matrix::Matrix( const char *csv,
				char delim)
{
    load(csv,delim);
}

Matrix::~Matrix()
{
	matrix.clear();
}

void Matrix::init()
{
	matrix = MatrixType(nrow, ncol);
}

//====================================================================
// data file loader
//====================================================================
void Matrix::load( const char *data, char delim )
{
    ifstream in( data );
    if(!in) {
        cerr << "File open error:" << data << "\n";
        exit(1);
    }

    size_t i = 0, r = 0;
    string line;
    while (getline(in, line)) {
        i++;
        if ( line[0]=='#' ) {
            if ( line.find("#row:") != string::npos ) 
                nrow = stoul(line.substr(5));
            else if ( line.find("#col:") != string::npos ) {
                ncol = stoul(line.substr(5));
                assert( nrow != 0 );
                assert( ncol != 0 );
                init();
            } else if ( line.find("genomes") != string::npos )
                parseHeader(line, delim);
            continue;
        }
        
        istringstream iss(line);
        string token;

        size_t c = 0;
		double v = 0.0;
		while(getline(iss, token, delim)) {
			size_t p = token.find(':');
			if ( p == string::npos ) continue; // rowname
			c = stoul(token.substr(0,p));
			assert(c <= ncol);
			v = stod(token.substr(p+1),nullptr);
			matrix(r,c-1) = v;
        }
        r++;
    }
    assert(r==nrow);
    
    in.close();
}


void Matrix::parseHeader(std::string &line, char delim)
{
    istringstream iss(line);
    string token;
    
    size_t c = 0;
    while(getline(iss, token, delim)) {
        if ( c > 0 ) {
			//colnames.push_back(token);
			size_t pos = token.find(':');
			assert( pos != string::npos );
			string con = token.substr(pos+1);
			colnames.push_back(con);
		}
        c++;
    }
    assert(c==ncol+1);
}
    


//====================================================================
// Sampling with replacement
//====================================================================
void Matrix::resample( MatrixType *sampled, size_t size, double seed )
{
    assert( size <= nrow );
    assert( sampled != NULL );
    srand(seed);
    for ( size_t i = 0; i < size; i++ ) {
        size_t k = rand()%nrow;
		for ( size_t j = 0; j < ncol; j++ ) {
			double v = matrix(k,j);
			if ( v != 0 ) 
				(*sampled)(i,j) = matrix(k,j);
		}
    }
}

//====================================================================
// Subsampling without replacement
//====================================================================
void Matrix::subsample( MatrixType *sampled, size_t size, double seed )
{
    assert( size <= nrow );
    assert( sampled != NULL );
    bool *picked = new bool[nrow];
    fill( picked, picked+nrow, false);
    size_t good = 0;
    
    srand(seed);
    while ( good < size ) {
        size_t k = rand()%nrow;
        if ( ! picked[k] ) {
			for ( size_t j = 0; j < ncol; j++ ) {
				double v = matrix(k,j);
				if ( v != 0 ) 
					(*sampled)(good,j) = matrix(k,j);
			}
            picked[k] = true;
            good++;
        }
    }
    delete[] picked;
}
