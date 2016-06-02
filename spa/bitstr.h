/** 
 * \file bitstr.h
 * Sequence to bit conversion.
 * This class converts a sequence of amino acids to an array of short integers.
 * Each short integer packs three amino acids.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \pre       Sequence consists of valid amino acids.
 * \bug       Not known.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __BITSTR_H__
#define __BITSTR_H__

#include <cstring>
#include <cmath>
#include "alpha.h"

//-----------------------------------------
// 3 Amino Acids per unsigned short integer
// 5 bits for each amino acid
// 1 bit left-over
//-----------------------------------------
typedef unsigned short Trimer; ///< 3 Amino acids
typedef unsigned short Length; ///< sequence length

const int NAA = 3; ///< Amino acid count

/**
 *
 */
class BitString
{
 private:
	Length size; ///< sequence length
	Trimer *str; ///< bit string
	
 private:
	/** 
	 * This calculates the number of integers needed to store amino acids.
	 * \return Integer.
	 */
	int getTrimerSize() 
	{
		return (int)ceil( (float)size/NAA );
	}

 public:
	/** Constructor 
	 */
	BitString() 
	{ 
		size = 0;
	}
	
	/** Constructor.
	 * Construct trimers for given string.
	 * \param Cstring	amino acid sequence.
	 */
	BitString(char *cstr) 
	{
		size = (Length)strlen(cstr);
		if ( size == 0 ) return;		

		int ntri = getTrimerSize();
		str = new Trimer[ntri];

		std::string sstr(cstr);
		int pos = 0;
		for ( int i = 0; i < size; i+=NAA ) {
			int naa = NAA;
			if ( i + NAA >= size ) naa = size-i; // last iteration
			std::string sub = sstr.substr(i, naa);
			str[pos++] = alpha::AminoAcidToInteger<Trimer>(sub);
		}
	}


	/**
	 * Copy constructor
	 */
	BitString( const BitString &source ) 
	{
		size = source.size;
		int ntri = getTrimerSize(); 
		str = new Trimer[ntri];
		for ( int i = 0; i < ntri; i++ )
			str[i] = source.str[i];
	}

	/**
	 * Assignment operator overloading
	 */
	BitString& operator= (const BitString &source)
	{
		if (this == &source) return *this;

		size = source.size;
		int ntri = getTrimerSize(); 
		str = new Trimer[ntri];
		for ( int i = 0; i < ntri; i++ )
			str[i] = source.str[i];

		return *this;		
	}

	/** 
	 * Destructor
	 * Release bit strings.
	 */
	~BitString() 
	{ 
		if ( size ) delete[] str; 
	}
	
	/**
	 * Restore sequence from bit string.
	 * \return String, amino acids.
	 */
	std::string toString() 
	{
		std::string sstr = "";
		int ntri = getTrimerSize();
		for ( int i = 0; i < ntri; i++ ) {
			int naa = NAA;
			if ( (i+1)*NAA >= size ) naa = NAA - ((i+1)*NAA-size); // last iteration
			sstr += alpha::IntegerToAminoAcid<Trimer>(str[i], naa);
		}
		return sstr;
	}

	/**
	 * Return sequence length.
	 * \return Integer, sequence size
	 */
	int getSize() 
	{
		return size;
	}
};

#endif
