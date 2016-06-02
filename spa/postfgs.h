/** 
 * \file      postfgs.h
 * \brief     FragGeneScan output handler
 * \details   This program takes FGS output and insert start and stop
 *            codons.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \pre       Run FragGeneScan beforehand.
 * \bug       None.
 * \warning   Improper use can crash your application.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __POST_FGS_H__
#define __POST_FGS_H__

#include <iostream>
#include <list>
#include <vector>
#include <string>
#include<boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
//#include <sstream>
#include "file.h"
#include "sequence.h"
#include "cmdargs.h"

/**< search direction ? */
enum directions { LEFT, RIGHT };
typedef std::list<std::string> sList;

const size_t MIN_MATCH = 10;

struct FgsGene
{
	int  rseq; // reference sequence id
	int  spos; // start position
	int  epos; // end position
	int  cins; // insertion count
	int  cdel; // deletion count;
	int  fram; // frame
	bool pstr; // positive strand?

	FgsGene()
	{
	}

	FgsGene(int r, int s, int e, int f, int i, int d, bool p) 
	{
		rseq = r; spos = s; epos = e; fram = f; cins = i; cdel = d; pstr = p;
	}

	FgsGene(const FgsGene &source)
	{
		rseq = source.rseq; 
		spos = source.spos;
		epos = source.epos;
		fram = source.fram;
		cins = source.cins;
		cdel = source.cdel;
		pstr = source.pstr;
	}

	FgsGene& operator= (const FgsGene &source ) 
	{
		if (this== &source) return *this;	
		rseq = source.rseq; 
		spos = source.spos;
		epos = source.epos;
		fram = source.fram;
		cins = source.cins;
		cdel = source.cdel;
		pstr = source.pstr;
		return *this;
	}
};


/* Global variable */
std::vector<std::string> input_files;
std::string fgs_out, fgs_ffn, fgs_faa, tmp_dir;
bool verbose = false;
	
void loadout(FgsGene*, char**, int); 
void findStopCodons(FgsGene*, int, char**, char**, bool*);
void writeNewFaa(char*, bool*, FgsGene*, int, char**);
void readParams(int argc, char **argv);
int getSize();

#endif
