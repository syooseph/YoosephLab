/** 
 * \file      sequence.h
 * \brief     FASTA sequence handler
 * \details   This includes character representation of amino acids,
 *            conversion between integer and amino acid.
 *            
 * \author    Youngik Yang
 * \version   0.2
 * \date      2010-2011
 * \bug       None.
 * \warning   None.
 * \date      Modified on Tue 2013-12-17 02:06:04 PM
 * \copyright J. Craig Venter Institute.
 */

#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include <list>
#include <sstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "kmer.h"
#include "alpha.h"
#include "utils.h"
#include "timer.h"
#include "file.h"

typedef std::string TagType;
typedef std::string SeqType;
typedef char* CTagType;
typedef char* CSeqType;

enum LOADMODE { TAGONLY, SEQONLY, TAGSEQ };

/**
 * \brief Sequence loader/checker
 */
namespace seq
{
	/** 
	 * Verify a sequence is not too short.
	 * \param tag sequence
	 * \param seq sequence tag
	 * \param k   minimum length
	 */
	inline bool checkLength(TagType &tag, SeqType &seq, int k) 
	{
		if ( seq.size() < (unsigned)k ) {
			std::cout << "Too short sequence (" << seq.size() << "): " <<  tag << "\n";
			return 0;
		}
		return 1;
	}
	
	/**
	 * Check validity of given sequence.
	 * 1. Valid amino acid letters.
	 * 2. No stop codon in first base.
	 * \param seq   sequence
	 * \param debug debug option
	 */
	inline bool validSequence(SeqType &seq, bool debug = false) 
	{
		std::vector<char> chars;
		std::vector<int>  where;
		for ( unsigned i = 0; i < seq.size()-1; i++ ) {
			int curr = alpha::getAAIndex(seq[i]);
			if (  curr  < 1 || curr > 25 ) {
				chars.push_back(seq[i]);
				where.push_back(i);
			}
		}

		// allow stop codon *
		int last = alpha::getAAIndex(seq[seq.size()-1]);
		if ( last < 1 || last > 26 ) {
			chars.push_back(seq[seq.size()-1]);
			where.push_back(unsigned(seq.size()-1));
		}
		
		if (chars.size() > 0) {
			if ( debug ) {
				std::cout << "Invalid letters:";
				for ( unsigned i = 0; i < chars.size(); i++ )
					std::cout << chars[i] << "(" << where[i] << ") ";
				std::cout << "\n";
			}
			return 0;
		}
		return 1;
	}
	
	/**
	 * Save sequence information to stroages.
	 * \param tags       2D array of sequence reads
	 * \param seqs       2D sequence tags
	 * \param fasta_tag  a sequence
	 * \param fasta_seq  a seqence tag
	 * \param index      array index to save record
	 * \param mode       load mode (sequence only, tag only, or both)
	 */
	inline void update( char **tags, 
				 char **seqs, 
				 TagType &fasta_tag, 
				 SeqType &fasta_seq, 
				 int &index, 
				 int mode)
	{
		if ( mode == TAGONLY || mode == TAGSEQ ) {
			int len = fasta_tag.size();
			tags[index] = new char[len];
			/* drop '>' */
			strncpy( tags[index], fasta_tag.substr(1, len-1).c_str(), len );
		}
		if ( mode == SEQONLY || mode == TAGSEQ ) {
			seqs[index] = new char[ fasta_seq.size()+1 ];
			strncpy( seqs[index], fasta_seq.c_str(), fasta_seq.size()+1 );
		}
	}
	

	/**
	 * Read a sequence file save all reads.
	 * \param filename read file name
	 * \param tags     2D array of sequence tags
	 * \param seqs     2D array of sequences
	 * \param count    sequence count
	 * \param mode     sequence load mode
	 * \param debug    debug option
	 */ 
	inline void readSequences( const char* filename, 
							   char **tags, 
							   char **seqs,
							   int &count,
							   int mode,
							   bool debug = false)
	{
		std::ifstream fstrm(filename, std::ios_base::in | std::ios_base::binary);
		
		if ( !fstrm ) {
			std::cerr << "Can't open " << filename << "\n";
			exit (1);
		}
		
		std::string line;
		TagType fasta_tag;
		SeqType fasta_seq;
		
		boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
		if ( fio::getFileExtension( std::string(filename) ) == "gz" ) 
			in.push(boost::iostreams::gzip_decompressor());
		in.push(fstrm);
		std::istream istrm(&in);

		while ( std::getline(istrm, line) ) {
			if ( line[0] == '>' ) {
				if ( fasta_tag != "" && fasta_seq != "" ) {
					update(tags, seqs, fasta_tag, fasta_seq, count, mode);
					count++;
				}
				fasta_tag = line; fasta_seq = ""; 
				if ( debug && count && count % 1000000 == 0 ) std::cout << "\t" << count << " sequences\n";
			}
			else fasta_seq += line;
		}
		if ( fasta_tag != "" && fasta_seq != "") {
			update(tags, seqs, fasta_tag, fasta_seq, count, mode);
			count++;
		}
		fstrm.close();
	}

	/**
	 * Get the number of sequences in a FASTA file
	 */
	inline int getSequenceCount( const char *filename )
	{
		std::fstream fstrm(filename, std::ios_base::in | std::ios_base::binary);
	
		if ( !fstrm ) {
			std::cerr << "Can't open " << filename << "\n";
			exit (1);
		}

		boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
		if ( fio::getFileExtension( std::string(filename) ) == "gz" ) 
			in.push(boost::iostreams::gzip_decompressor());
		in.push(fstrm);
		std::istream istrm(&in);

		int count = 0;
		std::string line;
		while ( std::getline(istrm, line) ) {
			if ( line[0] == '>' ) count++;
		}
	
		fstrm.close();
		return count;
	}

	/**
	 * Get a sequence count of multiple FASTA files
	 */
	inline int totalSequenceCount( const std::vector<std::string> &input_files )
	{
		int total = 0;
		for ( size_t i = 0; i < input_files.size(); i++ ) {
			total += getSequenceCount( input_files[i].c_str() );//, zflag );
		}
		return total;
	}


	/** 
	 * Load a sequence file
	 */
	inline void loadSequenceFile( std::string input_file, 
								  char **tags,
								  char **seqs,
								  int &count,
								  int mode,
								  bool debug = false) 
	{
		std::string file = input_file;
		readSequences(file.c_str(), tags, seqs, count, mode, debug);//, zflag );
	}

	/**
	 * Load multiple sequence files
	 */
	inline void loadSequences( const std::vector<std::string> &input_files, 
							   char **tags,
							   char **seqs,
							   int mode,
							   bool debug = false)
	{
		double time1 = mytime();
	
		if (debug) std::cout << "\nLoading sequences ...\n";
		int count = 0;
		for ( unsigned i = 0; i < input_files.size(); i++ ) {
			//std::cout << "File:" << (i+1) << "\n";
			loadSequenceFile( input_files[i], tags, seqs, count, mode );//, zflag );
		}
		if ( debug ) 
			std::cout << count << " Sequences loaded: ("
					  << mytime()-time1 << " sec)\n";
	}

	/** 
	 * Release memory
	 */
	inline void purge( char **mems, int len ) 
	{
		for ( int i = 0; i < len; i++ )
			delete[] mems[i];
		delete[] mems;
		mems = NULL;
	}
}

#endif

