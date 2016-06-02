#ifndef __REVERSE_TRANSLATE_H__
#define __REVERSE_TRANSLATE_H__

#include <iostream>
#include <string>
#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/copy.hpp>
#include <tr1/unordered_map>
#include "cmdargs.h"
#include "sequence.h"
#include "translation.h"
#include "read.h"
#include "profile.h"
#include "score.h"
#include "biostr.h"

//#include "spa.h"

typedef std::tr1::unordered_map<unsigned, std::vector<int> > IndelMap;
//typedef std::tr1::unordered_map<char*, char*> DnaReadMap;

struct Placement
{
	unsigned	id;
	std::string	seq;
	std::list<ReadId> reads;
	std::list<int> starts;
	IndelMap inss, dels;
	
	void clear() {
		seq.clear();
		reads.clear();
		starts.clear();
		inss.clear();
		dels.clear();
	}
};


const size_t NDNA = 5;
const char DNA[NDNA] = {'A', 'C', 'G', 'T', 'N'};


// Gene calling methods
//enum { FGA, MGA, BOTH };

const int BADPOS = -1000000;

//----------------------
// some global variables
//----------------------
std::string faa_file;
std::string ffn_file;
std::string place_file;
std::string out_file;


void trimReference( char **dnas, 
					char **orfs, 
					int norfs );

void generateDnaSeqs( char **orfs,
					  char **dnas );

void rtranslate( Placement &place,
				 char **orfs,
				 char **dnas,
				 std::fstream &out);

void addRead( std::string &line,
              Placement &place );

std::vector<int> getIndel( std::string &col );

void readParams(int argc, char **argv);

#endif
