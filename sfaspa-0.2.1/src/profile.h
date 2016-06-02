/**
 * \file      profile.h
 * \brief     Profile matrix
 * \details   Profile matrix of multiple DNA or amino acids sequences.
 * \author    Youngik Yang
 * \version   0.001
 * \date      Tue 2013-12-17 06:29:58 PM
 * \copyright J. Craig Venter Institute.
 */

#ifndef __PROFILE_H__
#define __PROFILE_H__

#include "msa.h"


//-------------------------------------------------------------------------------------------------------------
// Amino Acids
// 'A'=1, 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'=20.
// 'B'=21 (Aspartic Acid, Asparagine), 'Z'=22 (Glutamic Acid, Glutamine), 'X'=23 (unknown), 
// 'U'=24, 'O'=25, '*'=26 (terminator)
// '-'=27 for gap
//-------------------------------------------------------------------------------------------------------------
const int naas = 27;
const int ndna = 6; // A,C,G,T,N,-

/**
 * \brief Profile
 */
class Profile
{
 private:
	unsigned nrow;
	unsigned ncol;
	unsigned **matrix;
	int type;
	bool     replaced;

	void printMatrix(std::ostream &out);
	void printAlphabets(std::ostream &out);

 public:
	Profile();
	Profile(MSA*, int);
	~Profile();
	void make(Sequences &padded_reads);
	void print(std::ostream &out);

   	void init(size_t c, int t);
	
	unsigned getRowCount() { return nrow; }
	unsigned getColCount() { return ncol; }
	
	void increment(size_t i, size_t j) { matrix[i][j]++; }
	unsigned get(size_t i, size_t j) { return matrix[i][j]; }

	void clean();

	/**
	 * Consensus module.
	 * Based on profile, report single most frequet AA in each column.
	 * If there is a tie involing stop codon, use stop as a consensus.
	 */
	std::string makeConsensus();

	bool replacedStopCodon() { return replaced; }
};

#endif
