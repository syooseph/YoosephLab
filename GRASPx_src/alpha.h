/** 
 * \file      alpha.h
 * \brief     Conversion between amino acid and alphabet. 
 * \details   This includes character representation of amino acids,
 *            conversion between integer and amino acid.
 *            
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \bug       None.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __ALPHA_H__
#define __ALPHA_H__

#include <iostream>
#include <locale> // toupper
#include <math.h>


/**
 * Conversion between amino acid and alphabet. 
 */
namespace alpha 
{

	/** 
	 * 5 bits for amino acid
	 */
	const static int DIGIT = 5;

	/**
	 * Amino acids count
	 */
	const static int COUNT_AA = 26;

	/// Amino Acids.
	/// 'A'=1, 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'=20. \n
	/// 'B'=21 (Aspartic Acid, Asparagine), 'Z'=22 (Glutamic Acid, Glutamine), 'X'=23 (unknown). \n
	/// 'U'=24, UGA can code for selenocysteine (Sec or U), viewed as the 21st amino acid. \n
	/// 'O'=25, UAG can code for pyrrolysine (Pyl or O), which is viewed as the 22nd amino acid. \n
	/// '*'=26 (terminator). \n
	/// '-' is used for a gap. \n
	/// '0' is used for non-sense character.
	const static char AminoAcid[] = 
	{ 
		'0', 'A', 'R', 'N', 'D', 
		'C', 'Q', 'E', 'G', 'H', 
		'I', 'L', 'K', 'M', 'F', 
		'P', 'S', 'T', 'W', 'Y', 
		'V', 'B', 'Z', 'X', 'U',
		'O','*', '-' 
	};  // 0 - non-sense, '-' - gap

	//===========================================
	// Use four chars for start codon translation
	//===========================================

	/**
	 * Ascii value of each amino acid
	 */
	const static int AsciiToAA[] = 
	{ 
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 0 - 9
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 10 - 19
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 20 - 29
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 30 - 39
		-1, -1, 26, -1, -1, 27, -1, -1, -1, -1, // 40 - 49
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 50 - 59
		-1, -1, -1, -1, -1,  1, 21,  5,  4,  7, // 60 - 69
		14,  8,  9, 10, -1, 12, 11, 13,  3, 25, // 70 - 79
		15,  6,  2, 16, 17, 24, 20, 18, 23, 19, // 80 - 89
		22 
	};


	/// Reduced Amino Acids.
	/// {{C},{G},{A},{T,S},{N},{D,E},{Q,K,R},{V,I,L,M},{W,F,Y},{H},{P}}
	/// 'B' (Aspartic Acid:D, Asparagine:N) => 'D', 
	/// 'Z' (Glutamic Acid:E, Glutamine:Q) => 'Q'
	/// 'X'=23 (unknown) => '-'
	/// 'U'=24, UGA can code for selenocysteine (Sec or U), viewed as the 21st amino acid. => '-'
	/// 'O'=25, UAG can code for pyrrolysine (Pyl or O), which is viewed as the 22nd amino acid. '-'
	/// '*'=26 (terminator). 
	/// '-' is used for a gap.
	const static char ReducedAA[] = 
	{ 
		'C', 'G', 'A', 'T', 'N', 
		'D', 'Q', 'V', 'W', 'H', 
		'P', '*', '-' 
	}; 

	const static int AsciiToReducedAA[] = 
	{ 
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 0 - 9
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 10 - 19
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 20 - 29
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 30 - 39
		-1, -1, 26, -1, -1, 27, -1, -1, -1, -1, // 40 - 49
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 50 - 59
		-1, -1, -1, -1, -1,  1, 21,  5,  4,  7, // 60 - 69
		14,  8,  9, 10, -1, 12, 11, 13,  3, 25, // 70 - 79
		15,  6,  2, 16, 17, 24, 20, 18, 23, 19, // 80 - 89
		22 
	};


	
	
	/**
	 * Character to Amino Acid index conversion.
	 * This converts a single character to amino acid index in AminoAcid[] array.
	 * \param  Char		valid amino acid character.
	 * \return Integer	an index of the amino acid.
	 */
	inline int getAAIndex(char c) 
	{
		int k = int(toupper(c));
		if ( k != 42 && k!= 45 && ( k < 65 || k > 90 ) ) {
			return 0;
		}
		return AsciiToAA[k];
	}

	/**
	 * <em>i</em>th amino acid.
	 * \param  Integer	index in AminoAcid[] array.
	 * \return Amino acid.
	 */
	inline char getAminoAcid(int i) 
	{
		return AminoAcid[i];
	}

	/**
	 * Amino acid value of a given character.
	 * \param  Char	amino acid.
	 * \return Amino acid value in AminoAcid[] array.
	 */
	inline char getAminoAcid(char c) 
	{
		return AminoAcid[getAAIndex(c)];
	}

	/**
	 * String to integer conversion.
	 * This converts a string to an integer, which save memoris.
	 * \param	String	amino acids
	 * \return	Integer, bit string of amino acids.
	 * \warning	String must consists of valid amino acid, only.
	 */
	template<typename T>
		inline T AminoAcidToInteger(std::string &str) 
	{
		T num = 0;
		for ( size_t i = 0; i < str.size(); i++ ) {
			num += getAAIndex(str[i]);
			if ( i < str.size()-1 )
				num = num << DIGIT;
		}
		return num;
	}
	
	/* /\** */
	/*  * Char array to integer conversion. */
	/*  * This converts a char array to an integer, which save memoris. */
	/*  * \param	Char-array	amino acids. */
	/*  * \param   Integer		array size. */
	/*  * \return	Integer, bit string of amino acids. */
	/*  * \warning	String must consists of valid amino acid, only. */
	/*  *\/ */
	/* template<typename T> */
	/* 	inline T AminoAcidToInteger(char str[], int size)  */
	/* { */
	/* 	T num; */
	/* 	for ( int i = 0; i < size; i++ ) { */
	/* 		num = alpha::getAAIndex(str[i]); */
	/* 		if ( i != size-1 ) num = num << DIGIT; */
	/* 	} */
	/* 	return num; */
	/* } */

	/**
	 * Integer to string conversion.
	 * This converts an integer to correspoding amino acids.
	 * \param  Template	number representing amino acids.
	 * \param  Integer	size of amino acids packed into the number.
	 * \return String, amino acids.
	 */
	template<typename T>
		inline std::string IntegerToAminoAcid(T num, int size) 
	{
		std::string aa;

		for ( int i = 1; i <= size; i++ ) {
			int k = num % (int)pow(2.0, DIGIT);
			aa = getAminoAcid(k) + aa;
			num = num >> DIGIT;
		}
		return aa;
	}


	/**
	 * Last amino acid.
	 * \param  Template	number representing amino acids.
	 * \param  Integer	size of amino acids packed into the number.
	 * \return Last amino acid.
	 */
	template<typename T>
		inline char getLastAminoAcid(T num, int size) 
	{
		std::string aa = IntegerToAminoAcid(num, size);
		return aa[ aa.size()-1 ];
	}

	/**
	 * First amino acid.
	 * \param  Template	number representing amino acids.
	 * \param  Integer	size of amino acids packed into the number.
	 * \return First amino acid.
	 */
	template<typename T>
		inline char getFirstAminoAcid(T num, int size) 
	{
		std::string aa = IntegerToAminoAcid(num, size);
		return aa[0];
	}
}

#endif
