#include "mg.h"

using namespace std;

MG::MG( const char *file)
{
    loadGenomeSizes(file);
}

MG::~MG()
{
    /* Default deconstructor is fine */
}

//====================================================================
// Load genome size
//====================================================================
void MG::loadGenomeSizes( const char *file )
{
    ifstream in( file );    
    if(!in) {
        cerr << "File open error:" << file << "\n";
        exit(1);
    }
    
    string line;
    while ( getline(in, line) ) {
        std::string genome;
        uint32_t size;
        
        istringstream iss(line);
        string token;

        iss >> genome >> token;
        size = stoul(token, nullptr, 0);        
        
        size_map.insert( std::pair<std::string, uint32_t>( genome, size ) );
    }
    in.close();
}

//====================================================================
// Compute abundances of genomes from given mixing coefficients
//====================================================================
void MG::computeAbundances(double *coefficient, std::vector<std::string> &names )
{
    size_t ngenome = names.size();

    if ( abud_map.size() ) abud_map.clear();
    
    double total_ratio = 0;
    for ( size_t j = 0; j < ngenome; j++ ) {
        auto it = size_map.find( names[j] );
        assert( it != size_map.end() );
        uint32_t size = it->second;
        double ratio = coefficient[j]/size;
        total_ratio += ratio;
    }

    for ( size_t j = 0; j < ngenome; j++ ) {
        std::string genome = names[j];
        uint32_t size = size_map[genome];
        double abud = coefficient[j]/(size*total_ratio);
        abud_map.insert( std::pair<std::string, double>( genome, abud ) );
    }
}

