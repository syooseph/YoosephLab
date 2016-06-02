/** 
 * \file      dna.h
 * \brief     DNA sequence handler
 * \author    Youngik Yang
 * \version   0.2
 * \date      2013
 * \bug       None.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __DNA_H__
#define __DNA_H__

#include <string>
#include <algorithm>

/**
 * \brief DNA definition
 */
namespace dna
{
	/**
	 * Reverse complement of given DNA sequences.
	 */
	inline std::string reverseComplement(std::string &str)
	{
		std::reverse(str.begin(), str.end());
		for ( size_t i = 0; i < str.length(); i++ ) {
			char base = str[i];
			char ndna = base;
			switch(base) {
			case 'A':
				ndna = 'T'; break;
			case 'C':
				ndna = 'G'; break;
			case 'G':
				ndna = 'C'; break;
			case 'T':
				ndna = 'A'; break;
			case 'a':
				ndna = 't'; break;
			case 'c':
				ndna = 'g'; break;
			case 'g':
				ndna = 'c'; break;
			case 't':
				ndna = 'a'; break;
			default:
				break;
			}
			str[i] = ndna;
		}
		return str;
	}
}

#endif
