#include "genome_matrix.h"
using namespace std;

GenomeMatrix::GenomeMatrix(const char *genome_file)
{
    loadGenomes(genome_file);
}

GenomeMatrix::~GenomeMatrix()
{
    bwa = NULL;
}

void GenomeMatrix::loadGenomes(const char *genome_file)
{
    ifstream in(genome_file, ios_base::in);
    if ( !in ) {
        cerr << "Can't open " << genome_file << "\n"; exit (1);
    }

    string line;
    while(getline(in,line)) {
        istringstream iss(line);
        vector<string> c;
        string tok;
        while(iss>>tok) c.push_back(tok);
        if ( c.size() < 2 ) {
            cerr << "Format error\n"; exit(1);
        }
		string contig = c[0], genome = c[1];
        con2gnm.insert( pair<string,string>(contig, genome) );
    }
    in.close();

	//-------------------
	// Extract genome set
	//-------------------
	std::set<std::string> gnmset;
	genomes.reserve(con2gnm.size());
	for ( auto it = con2gnm.begin(); it != con2gnm.end(); ++it ) 
		gnmset.insert(it->second);


	//-----------------------------------------------------
	// Sort genomes
	// Should be in order since set keeps elements in order
	//-----------------------------------------------------
	genomes = std::vector<std::string>( gnmset.begin(), gnmset.end() );

	//------------------------------
	// Create genome to number index
	//------------------------------
	for ( size_t i = 0; i < genomes.size(); i++ ) 
		gindexs[genomes[i]] = i+1;
}

void GenomeMatrix::generate()
{
    getGenomeLengths();
    writeFiles();
}

void GenomeMatrix::getGenomeLengths()
{
    NumberToContig* contigs = bwa->getContigNames();
    SbjNumLengths*  conlens = bwa->getSbjctLengths();

    for ( auto it = conlens->begin(); it != conlens->end(); ++it ) {
        sbjct_t contig = it->first;
        size_t  conlen = it->second;

        //-----------------
        // Find contig name
        //-----------------
        auto jt = contigs->find(contig);
		if ( jt == contigs->end() ) {
			fprintf(stderr, "[ERROR] contig:%u does not exist\n", contig);
		}
        //assert( jt != contigs->end() );
        string cname = jt->second;

        //-----------------
        // Find genome name
        //-----------------        
        auto kt = con2gnm.find(cname);
        //assert( kt != con2gnm.end() );
		if ( kt == con2gnm.end() ) {
			fprintf(stderr, "[ERROR] genome name for contig:%u (%s) does not exist\n", contig, cname.c_str());
		}
        string gname = kt->second;

        //---------------------
        // Update genome length
        //---------------------
        auto lt = gnmlens.find(gname);
        if ( lt == gnmlens.end() ) 
            gnmlens.insert( pair<string,size_t>(gname, conlen) );
        else lt->second += conlen;
    }
}


void GenomeMatrix::writeFiles()
{
    writeGenomeMatrix();
    writeSummary();
    writeLengths();
}

void GenomeMatrix::writeGenomeMatrix()
{
	double t0 = mytime();
	
    NumberToContig* contigs = bwa->getContigNames(); // unordered map
    QueryMatches*   queries = bwa->getHits();        // unordered map
    
    string mat_file = outdir + "/gabe_matrix";    
    ofstream out(mat_file.c_str(), ios_base::out);
    if ( !out ) {
        cerr << "Can't open " << mat_file << "\n"; exit (1);
    }

    
    out << "#row:" << queries->size() << "\n";
    out << "#col:" << genomes.size() << "\n";
    
    out << "#genomes";
	for ( size_t i = 0; i < genomes.size(); i++ )
		out << "\t" << i+1 << ":" << genomes[i];
    out << "\n";


	//------------------------------
	// Record genome counts in order
	//------------------------------
	size_t nquery = queries->size();
    map<string,double> counts;
	for ( size_t i = 1; i <= nquery; i++ ) {
		auto it = queries->find(i);
		assert( it != queries->end() );
        counts.clear();

        //------------------------------------------------------------
        // Get counts in genomic level by summing up counts in contigs
        //------------------------------------------------------------
        SbjctList conlist = it->second;
        for ( const auto &jt : conlist ) {
            auto kt = contigs->find(jt);
            assert( kt != contigs->end() );
            string  cname = kt->second;

            auto gt = con2gnm.find(cname);
            assert( gt != con2gnm.end() );
            string gname = gt->second;
            
            auto st = counts.find(gname);
            if ( st == counts.end() ) 
                counts.insert( pair<string,double>(gname, 1.0) );
            else st->second++;
        }

        //--------------------------------------------
        // Now write probability (count/genome_length)
        //--------------------------------------------
        out << it->first;
		for ( auto ct = counts.begin(); ct != counts.end(); ++ct ) {
			double count = ct->second;
			string gname = ct->first;
			size_t index = gindexs[gname];
			size_t gmlen = gnmlens[gname];
			out << "\t" << index << ":" << count/gmlen;
		}
		out << "\n";
    }
    out.close();

	log( "Genome matrix written", t0 );
}


void GenomeMatrix::writeLengths()
{
    string len_file = outdir + "/gabe_lengths";    
    ofstream out(len_file.c_str(), ios_base::out);
    if ( !out ) {
        cerr << "Can't open " << len_file << "\n"; exit (1);
    }

	for ( size_t i = 0; i < genomes.size(); i++ ) 
        out << genomes[i] << "\t" << gnmlens[genomes[i]] << "\n";

    out.close();
}


void GenomeMatrix::writeSummary()
{
    string summary_file = outdir + "/gabe_preliminaries";
    ofstream out(summary_file.c_str(), ios_base::out);
    if ( !out ) {
        cerr << "Can't open " << summary_file << "\n"; exit (1);
    }

    unordered_map<string, double> abunds1, abunds2;
	getAbundanceByEvenSplit(abunds1);
	getAbundanceByMultCount(abunds2);
 
    out << "#abundance1: multiple best hits (N) were counted multiple times (N).\n";
    out << "#abundance2: multiple best hits (N) were distributed evenly (1/N).\n";
    out << "# genome\tabundance1\tabundance2\n";
	for ( size_t i = 0; i < genomes.size(); i++ ) {
        string genome = genomes[i];
        double abund1  = abunds1[genome];
        double abund2  = abunds2[genome];
        out << genome << "\t" << abund1 << "\t" << abund2 << "\n";
    }
}

//====================================================================
// 1 over N strategy
//====================================================================
void GenomeMatrix::getAbundanceByEvenSplit( unordered_map<string, double> &abunds )
{
    QueryMatches*   querys  = bwa->getHits();
    NumberToContig* contigs = bwa->getContigNames();

    unordered_map<string, float> ratios;
    for ( auto jt = querys->cbegin(); jt != querys->cend(); ++jt ) {
        SbjctList sbjs = jt->second;
        size_t count = sbjs.size();

        for ( const auto &jt : sbjs ) {
            auto tt = contigs->find(jt);
            assert( tt != contigs->end() );
            string cname = tt->second;

            auto gt = con2gnm.find(cname);
            assert( gt != con2gnm.end() );
            string genome = gt->second;

            auto rt = ratios.find(genome);
            if ( rt == ratios.end() )
                ratios.insert( pair<string, float>( genome, 1.0/count ) );
            else rt->second += 1.0/count;
        }
    }

    double sum_ratio= 0, nquery = (double)querys->size();
    unordered_map<string, double> coeffs;    
    for ( auto it = ratios.cbegin(); it != ratios.cend(); ++it ) {
        string gname = it->first;        
        double  coeff = it->second / nquery;
        size_t length = gnmlens[gname];
        coeffs.insert( pair<string, double>(gname, coeff) );
        sum_ratio += coeff/(double)length;
    }

    for ( auto it = gnmlens.cbegin(); it != gnmlens.cend(); ++it ) {
        string genome = it->first;
        double abd = coeffs[genome] / ( gnmlens[genome]*sum_ratio );
        abunds.insert( pair<string, double>( genome, abd ) );
    }
}

void GenomeMatrix::getAbundanceByMultCount( unordered_map<string, double> &abunds )
{
    SbjctCounts*    sbjcts  = bwa->getSbjHits();
    NumberToContig* contigs = bwa->getContigNames();

    //-------------------------------------------
    // Total count by allowing multiple countings
    //-------------------------------------------
    size_t total = 0;
    for ( auto it = sbjcts->cbegin(); it != sbjcts->cend(); ++it ) 
        total += it->second;

    // Multiple assignment
    double sum_ratio = 0;
    unordered_map<string, double> coeffs;
    for ( auto it = sbjcts->cbegin(); it != sbjcts->cend(); ++it ) {
        sbjct_t sbjct = it->first;
        auto tt = contigs->find(sbjct);
        string cname = tt->second;
        string genome = con2gnm[cname];
        double  coeff = it->second / (double)total;
        auto ct = coeffs.find(genome);
        if ( ct == coeffs.end() ) 
            coeffs.insert( pair<string, double>(genome, coeff) );
        else ct->second += coeff;
    }

    for ( auto it = coeffs.cbegin(); it != coeffs.cend(); ++it ) {
        string gname = it->first;
        double coeff = it->second;
        sum_ratio += coeff/gnmlens[gname];
    } 

    for ( auto it = gnmlens.cbegin(); it != gnmlens.cend(); ++it ) {
        string genome = it->first;
        double abd = coeffs[genome] / ( gnmlens[genome]*sum_ratio );
        abunds.insert( pair<string, double>( genome, abd ) );
    }
}

