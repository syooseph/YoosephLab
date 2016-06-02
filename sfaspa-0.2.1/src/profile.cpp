#include "profile.h"

Profile::Profile()
{
    nrow = ncol = 0;
    matrix = NULL;
    replaced = false;
}

Profile::Profile(MSA *msa, int t)
{
    replaced = false;

    std::string pivot   = msa->getPivot();
    Sequences   members = msa->getSequences();

    init(pivot.size(), t);
    make(members);
}

Profile::~Profile()
{
    clean();
}

void Profile::init(size_t c, int t)
{
    type = t;
    if ( c == 0 ) return;
    assert( type == AA_TYPE || type == DNA_TYPE );

    nrow = c;
    ncol = (type==AA_TYPE) ? naas : ndna;

    matrix = new unsigned*[nrow];
    for ( size_t i = 0; i < nrow; i++ ) 
        matrix[i] = new unsigned[ncol];
    
    for ( size_t i = 0; i < nrow; i++ ) 
        for ( size_t j = 0; j < ncol; j++ )
            matrix[i][j] = 0;
}

void Profile::clean() 
{
    if ( matrix == NULL ) return;
    for ( size_t i = 0; i < nrow; i++ )
        delete[] matrix[i];
    delete[] matrix;
    matrix = NULL;
    nrow = ncol = 0;
}

void Profile::make(Sequences &padded_reads)
{
    for ( size_t i = 0; i < padded_reads.size(); i++ ) {
        for ( size_t j = 0; j < nrow; j++ ) {
            char ch = padded_reads[i][j];
            if ( ch == '.' ) continue;
            int nc;
            if ( type == AA_TYPE ) {
                nc = alpha::AsciiToAA[(int)ch];
                if ( nc < 1 || nc > 27 ) 
                    std::cerr << "Invalid aa:" << ch << "\tnum:" << nc << "\tcol:" << j << "\t" << padded_reads[i] << "\n";
                else
                    matrix[j][nc-1]++;
            } else {
                nc = alpha::AsciiToDNA[(int)ch];
                if ( nc == -1 ) 
                    std::cout << "Invalid dna:" << ch << "\tnum:" << nc << "\tcol:" << j << "\t" << padded_reads[i] << "\n";
                else
                    matrix[j][nc]++;
            }
        }
    }
}

void Profile::printAlphabets( std::ostream &out )
{
    char buf[25];
    out << "\t";    
    if ( type == AA_TYPE ) {
        for ( int i = 0; i < naas; i++ ) {
            sprintf(buf, "%4c ", alpha::AminoAcid[i+1]);
            out << buf;
        }
	} else {
        for ( int i = 0; i < ndna; i++ ) {
            sprintf(buf, "%4c ", alpha::DNA[i]);
            out << buf;
        }
    }
    out << "\n";
}

void Profile::printMatrix(std::ostream &out)
{
    char buf[25];
    for ( size_t i = 0; i < nrow; i++ ) {
        out << i << "\t";
        for ( size_t j = 0; j < ncol; j++ ) {
            sprintf(buf, "%4d ", matrix[i][j]);
            out << buf;
        }
        out << "\n";
    }
}

void Profile::print(std::ostream &out)
{
    printAlphabets(out);
    printMatrix(out);
}

std::string Profile::makeConsensus()
{
    std::string consensus = std::string(nrow, '.');
    for ( size_t i = 0; i < nrow; i++ ) {
        std::multimap<int, int> cmap;
        for ( size_t j = 0; j < (size_t)ncol; j++ ) 
            cmap.insert( std::pair<int, int>(matrix[i][j], j) );
        
        std::multimap<int, int>::reverse_iterator it = cmap.rbegin();
        int max = it->first;
        if ( max == 0 ) continue;

        if ( type == AA_TYPE ) assert(it->second>=0 && it->second<=26);
        if ( type == DNA_TYPE ) assert(it->second>=0 && it->second<=5);
        consensus[i] = (type==AA_TYPE) ? 
            alpha::AminoAcid[(it->second)+1] :
            alpha::DNA[it->second] ;
        
        if ( type == AA_TYPE && consensus[i] == '*' && Param::replace_stop ) {
            if ( i == nrow-1 ) continue; // last one - leave it

            ++it;
            assert(it->second>=0 && it->second<=26);
            if ( it->first > 0 ) {
                consensus[i] = alpha::AminoAcid[(it->second)+1];
                replaced = true;
            }
        }
    }
    return consensus;
}

