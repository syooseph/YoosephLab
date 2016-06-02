#ifndef __TRANSLATION_H__
#define __TRANSLATION_H__

#include <tr1/unordered_map>
#include "dna.h"

typedef std::tr1::unordered_map<std::string, char> CodonMap;

struct Codon
{
	CodonMap table;

	Codon()
	{
		init();
	}

	void init()
	{
		table.insert( std::pair<std::string, char>("TTT", 'F')); // Phenyl-alanine
		table.insert( std::pair<std::string, char>("TTC", 'F')); 
		table.insert( std::pair<std::string, char>("TTA", 'L')); // Leucine
		table.insert( std::pair<std::string, char>("TTG", 'L')); // Leucine, Start codon

		table.insert( std::pair<std::string, char>("TCT", 'S')); // Serine
		table.insert( std::pair<std::string, char>("TCC", 'S'));
		table.insert( std::pair<std::string, char>("TCA", 'S'));
		table.insert( std::pair<std::string, char>("TCG", 'S'));

		table.insert( std::pair<std::string, char>("TAT", 'Y')); // Tyrosine
		table.insert( std::pair<std::string, char>("TAC", 'Y'));
		table.insert( std::pair<std::string, char>("TAA", '*')); // Stop codon
		table.insert( std::pair<std::string, char>("TAG", '*')); // Stop codon

		table.insert( std::pair<std::string, char>("TGT", 'C')); // Cysteine
		table.insert( std::pair<std::string, char>("TGC", 'C'));
		table.insert( std::pair<std::string, char>("TGA", '*')); // Stop codon
		table.insert( std::pair<std::string, char>("TGG", 'W')); // Tryptophan


		table.insert( std::pair<std::string, char>("CTT", 'L')); // Leucine
		table.insert( std::pair<std::string, char>("CTC", 'L'));
		table.insert( std::pair<std::string, char>("CTA", 'L'));
		table.insert( std::pair<std::string, char>("CTG", 'L')); // Leucine, Start codon

		table.insert( std::pair<std::string, char>("CCT", 'P')); // Proline
		table.insert( std::pair<std::string, char>("CCC", 'P'));
		table.insert( std::pair<std::string, char>("CCA", 'P'));
		table.insert( std::pair<std::string, char>("CCG", 'P'));

		table.insert( std::pair<std::string, char>("CAT", 'H')); // Histidine
		table.insert( std::pair<std::string, char>("CAC", 'H'));
		table.insert( std::pair<std::string, char>("CAA", 'Q')); // Glutamine
		table.insert( std::pair<std::string, char>("CAG", 'Q'));

		table.insert( std::pair<std::string, char>("CGT", 'R')); // Arginine
		table.insert( std::pair<std::string, char>("CGC", 'R'));
		table.insert( std::pair<std::string, char>("CGA", 'R'));
		table.insert( std::pair<std::string, char>("CGG", 'R'));


		table.insert( std::pair<std::string, char>("ATT", 'I')); // Isoleucine, Start codon
		table.insert( std::pair<std::string, char>("ATC", 'I'));
		table.insert( std::pair<std::string, char>("ATA", 'I'));
		table.insert( std::pair<std::string, char>("ATG", 'M')); // Methionine, Start codon

		table.insert( std::pair<std::string, char>("ACT", 'T')); // Threonine
		table.insert( std::pair<std::string, char>("ACC", 'T'));
		table.insert( std::pair<std::string, char>("ACA", 'T'));
		table.insert( std::pair<std::string, char>("ACG", 'T'));

		table.insert( std::pair<std::string, char>("AAT", 'N')); // Asparagine
		table.insert( std::pair<std::string, char>("AAC", 'N'));
		table.insert( std::pair<std::string, char>("AAA", 'K')); // Lysine
		table.insert( std::pair<std::string, char>("AAG", 'K'));

		table.insert( std::pair<std::string, char>("AGT", 'S')); // Serine
		table.insert( std::pair<std::string, char>("AGC", 'S'));
		table.insert( std::pair<std::string, char>("AGA", 'R')); // Arginine
		table.insert( std::pair<std::string, char>("AGG", 'R'));



		table.insert( std::pair<std::string, char>("GTT", 'V')); // Valine
		table.insert( std::pair<std::string, char>("GTC", 'V'));
		table.insert( std::pair<std::string, char>("GTA", 'V'));
		table.insert( std::pair<std::string, char>("GTG", 'V')); // Valine, Start codon 

		table.insert( std::pair<std::string, char>("GCT", 'A')); // Alanine
		table.insert( std::pair<std::string, char>("GCC", 'A'));
		table.insert( std::pair<std::string, char>("GCA", 'A'));
		table.insert( std::pair<std::string, char>("GCG", 'A'));

		table.insert( std::pair<std::string, char>("GAT", 'D')); // Aspartic acid
		table.insert( std::pair<std::string, char>("GAC", 'D'));
		table.insert( std::pair<std::string, char>("GAA", 'E')); // Glutamic acid
		table.insert( std::pair<std::string, char>("GAG", 'E'));

		table.insert( std::pair<std::string, char>("GGT", 'G')); // Glycine
		table.insert( std::pair<std::string, char>("GGC", 'G'));
		table.insert( std::pair<std::string, char>("GGA", 'G'));
		table.insert( std::pair<std::string, char>("GGG", 'G'));
	}
};

namespace tran
{
	Codon codon;

	inline std::string translate(std::string &dna_str)
	{
/* 		if ( dna_str == "" ) return ""; */
		if ( dna_str.size() < 3 ) return "";

/* 		if ( dna_str.length() % 3 != 0 ) { */
/* 			std::cerr << "[Error] sequence length is not multiple of 3\n"; */
/* 			return ""; */
/* 		} */

		std::string aa_str;
		for ( size_t i = 0; i < dna_str.length(); i += 3 ) {
			if ( i+3 > dna_str.size() ) break;
			if ( codon.table.find(dna_str.substr(i, 3)) == codon.table.end() ) 
				aa_str += 'X';
			else
				aa_str += codon.table[dna_str.substr(i, 3)];
		}
		return aa_str;
	}
/* 	inline bool isStartCodon( char aa ) */
/* 	{ */
		
/* 	} */
}

#endif
