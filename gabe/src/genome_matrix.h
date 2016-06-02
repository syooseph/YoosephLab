#ifndef __GENOME_MATRIX__
#define __GENOME_MATRIX__

#include "bwa_reader.h"

using ContigToGenome = std::unordered_map<std::string, std::string>;
using GenomeLengths  = std::unordered_map<std::string, size_t>; 
using GenomeIndices  = std::unordered_map<std::string, size_t>;

class GenomeMatrix
{
 private:
	BwaReader*     bwa;
	ContigToGenome con2gnm;
	GenomeLengths  gnmlens;

	std::vector<std::string> genomes;
	GenomeIndices  gindexs;
	std::string outdir;
	
	void loadGenomes(const char*);
	void getGenomeLengths();
	void writeGenomeMatrix();
	void writeSummary();
	void writeLengths();
	void writeFiles();
	void getAbundanceByEvenSplit( std::unordered_map<std::string, double> & );
	void getAbundanceByMultCount( std::unordered_map<std::string, double> & );
 public:
	GenomeMatrix(const char*);
	~GenomeMatrix();
	
	void generate();
	void setBwaReader( BwaReader *b )   { bwa = b; }
	void setDirectory( std::string& o ) { outdir = o; }
};

#endif
