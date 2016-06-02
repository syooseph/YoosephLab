#ifndef __PROFILE_H__
#define __PROFILE_H__

#include <iostream>
#include <map>
#include "alpha.h"
#include <stdio.h>
#include <assert.h>

//-------------------------------------------------------------------------------------------------------------
// Amino Acids
// 'A'=1, 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'=20.
// 'B'=21 (Aspartic Acid, Asparagine), 'Z'=22 (Glutamic Acid, Glutamine), 'X'=23 (unknown), 
// 'U'=24, 'O'=25, '*'=26 (terminator)
// '-'=27 for gap
//-------------------------------------------------------------------------------------------------------------
const int naas = 27;

struct Profile
{
	unsigned ncol;
	unsigned** matrix;

	Profile()
	{
		ncol = 0;
		matrix = NULL;
	}
	void init(unsigned n) 
	{
		if ( n == 0 ) return;
		ncol = n;
		matrix = new unsigned*[ncol];
		for ( size_t i = 0; i < ncol; i++ ) 
			matrix[i] = new unsigned[naas];

		for ( size_t i = 0; i < ncol; i++ ) 
			for ( int j = 0; j < naas; j++ )
				matrix[i][j] = 0;
	}

	void clean() 
	{
		if ( ncol == 0 ) return;
		for ( size_t i = 0; i < ncol; i++ )
			delete[] matrix[i];
		delete[] matrix;
		matrix = NULL;
		ncol = 0;
	}
	
	void reset()
	{
		for ( size_t i = 0; i < ncol; i++ ) 
			for ( int j = 0; j < naas; j++ ) 
				matrix[i][j] = 0;
	}

	
	void update(std::vector<std::string> &padded_reads)
	{
		reset();
		for ( size_t i = 0; i < padded_reads.size(); i++ ) {
			for ( size_t j = 0; j < ncol; j++ ) {
				char aa = padded_reads[i][j];
				if ( aa == '.' ) continue;
				int na = alpha::AsciiToAA[(int)aa];
				if ( na < 1 || na > 27 ) 
					std::cerr << "Invalid aa:" << aa << "\tnum:" << na << "\tcol:" << j << "\t" << padded_reads[i] << "\n";
				matrix[j][na-1]++;
			}
		}
	}

	void print(std::ostream &out)
	{
		char buf[25];
		
		//------------
		// amino acids
		//------------
		out << "\t";    
		for ( int i = 0; i < naas; i++ ) {
			sprintf(buf, "%3c ", alpha::AminoAcid[i+1]);
			out << buf;
		}
		out << "\n";
		
		//-------
		// matrix
		//------- 
		for ( size_t i = 0; i < ncol; i++ ) {
			out << i << "\t";
			for ( int j = 0; j < naas; j++ ) {
				sprintf(buf, "%3d ", matrix[i][j]);
				out << buf;
			}
			out << "\n";
		}
	}
	
	/**
	 * Consensus module.
	 * Based on profile, report single most frequet AA in each column.
	 * If there is a tie involing stop codon, use stop as a consensus.
	 */
	std::string makeConsensus()
	{
		std::string consensus = std::string(ncol, '.');
		for ( size_t i = 0; i < ncol; i++ ) {
			std::multimap<int, int> cmap;
			for ( size_t j = 0; j < (size_t)naas; j++ ) 
				cmap.insert( std::pair<int, int>(matrix[i][j], j) );
			
			std::multimap<int, int>::reverse_iterator it = cmap.rbegin();
			int max = it->first;
			if ( max == 0 ) continue;
			
			consensus[i] = alpha::AminoAcid[(it->second)+1];
			if ( consensus[i] != '*' ) {
				++it;
				for ( ; it != cmap.rend(); ++it ) {
					if ( it->first < max ) break;
					if ( alpha::AminoAcid[(it->second)+1] == '*' ) {
						consensus[i] = alpha::AminoAcid[(it->second)+1];
						break;
					}
				}
			}
		}
		return consensus;
	}

};

#endif
