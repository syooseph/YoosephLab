/** 
 * \file      codon.h
 * \brief     Codon table
 * \details   Codon table creation and translation of DNA sequence.
 *            
 * \author    Youngik Yang
 * \version   0.2
 * \date      2013
 * \bug       None.
 * \warning   None.
 * \date      Modified on Fri 2013-12-13 02:10:26 PM
 * \copyright J. Craig Venter Institute.
 */

#ifndef __CODON_H__
#define __CODON_H__

#include <tr1/unordered_map>
#include "dna.h"


typedef std::tr1::unordered_map<std::string, char> CodonMap;

/**
 * \class Codon
 * \brief Condon table and translation from DNA to AA
 */
class Codon
{
 private:
	CodonMap table; ///< codon table

 private:
	/**
	 * Fill condon table
	 */
	void init();

 public:
	/**
	 * Constructor 
	 */
	Codon();

	/** 
	 * Translate DNA sequence to amino acids sequence
	 */
	std::string translate(std::string &dna_str);
};

#endif
