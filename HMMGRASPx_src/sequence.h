#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include <ctype.h>
#include <list>
#include <tr1/unordered_map>
//#include <seqan/file.h>
#include <math.h>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

//#include "datatype.h"
//#include "kmer.h"
#include "alpha.h"
#include "utils.h"
#include "timer.h"
#include "file.h"

typedef std::string TagType;
typedef std::string SeqType;
typedef char* CTagType;
typedef char* CSeqType;

enum LOADMODE { TAGONLY, SEQONLY, TAGSEQ };

namespace seq
{
	inline bool checkLength(TagType &tag, SeqType &seq, int k) 
	{
		if ( seq.size() < (unsigned)k ) {
			std::cout << "Too short sequence (" << seq.size() << "): " <<  tag << "\n";
			return 0;
		}
		return 1;
	}
	
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
		//index++;
	}
	
/* 	//============================================================================== */
/* 	// Tue 2011-05-24 10:13:05 AM */
/* 	// O.K. */
/* 	//============================================================================== */
/* 	void readSequences( const char* filename,  */
/* 						char **tags,  */
/* 						char **seqs, */
/* 						int &count, */
/* 						int mode ) */
/* 	{ */
/* 		std::fstream fstrm(filename, std::ios_base::in | std::ios_base::binary); */
		
/* 		if ( !fstrm ) { */
/* 			std::cerr << "Can't open " << filename << "\n"; */
/* 			exit (1); */
/* 		} */
		
/* 		//int nread = 0; */
/* 		std::string line; */
/* 		TagType fasta_tag; */
/* 		SeqType fasta_seq; */
		
/* 		while ( !fstrm.eof() ) { */
/* 			getline(fstrm, line); */
/* 			if ( line[0] == '>' ) { */
/* 				if ( fasta_tag != "" && fasta_seq != "" ) { */
/* 					update(tags, seqs, fasta_tag, fasta_seq, count, mode); */
/* 					count++; */
/* 				} */
/* 				fasta_tag = line; fasta_seq = "";  */
/* 				if ( count && count % 1000000 == 0 ) std::cerr << "\t" << count << " reads\n"; */
/* 			} */
/* 			else fasta_seq += line; */
/* 		} */
/* 		if ( fasta_tag != "" && fasta_seq != "") { */
/* 			update(tags, seqs, fasta_tag, fasta_seq, count, mode); */
/* 			count++; */
/* 		} */
/* 		fstrm.close(); */
/* 	} */
/* 	void readGzipSequences( const char* filename,  */
/* 							char **tags,  */
/* 							char **seqs, */
/* 							int &count, */
/* 							int mode ) */
/* 	{ */
/* 		std::fstream fstrm(filename, std::ios_base::in | std::ios_base::binary); */
	
/* 		if ( !fstrm ) { */
/* 			std::cerr << "Can't open " << filename << "\n"; */
/* 			exit (1); */
/* 		} */
 
/* 		boost::iostreams::filtering_streambuf<boost::iostreams::input> in; */
/* 		in.push(boost::iostreams::gzip_decompressor()); */
/* 		in.push(fstrm); */
   
/* 		//int count = 0; */
/* 		std::string line; */
/* 		TagType fasta_tag; */
/* 		SeqType fasta_seq; */
    
/* 		std::istream istrm(&in); */
/* 		while ( std::getline(istrm, line) ) { */
/* 			if ( line[0] == '>' ) { */
/* 				if ( fasta_seq != "" ) { */
/* 					update(tags, seqs, fasta_tag, fasta_seq, count, mode); */
/* 					count++; */
/* 				} */
/* 				fasta_tag = line; fasta_seq = "";  */
/* 				if ( count && count % 1000000 == 0 ) std::cerr << "\t" << count << " reads\n"; */
/* 			} */
/* 			else fasta_seq += line; */
/* 		} */
/* 		if ( fasta_tag != "" && fasta_seq != "") { */
/* 			update(tags, seqs, fasta_tag, fasta_seq, count, mode); */
/* 			count++; */
/* 		} */
/* 		fstrm.close(); */
/* 	} */
  inline void reverseSequences(unsigned int num_seqs, char** seqs, char** rev_seqs)  {
    unsigned int i, j;
    for(i = 0; i < num_seqs; ++ i) {
      unsigned int len = strlen(seqs[i]);
      rev_seqs[i] = new char[len + 1];
      for(j = 0; j < len; ++ j) {
        rev_seqs[i][len - j - 1] = seqs[i][j];
      }
      rev_seqs[i][len] = '\0';
    } 
    return;
  }
  
  inline void reAssignSpecialChar(SeqType& seq) {
    unsigned int i;
    for(i = 0; i < seq.length(); ++ i) {
      if(seq[i] != 'A' && seq[i] != 'R' && seq[i] != 'N' && seq[i] != 'D' &&
         seq[i] != 'C' && seq[i] != 'Q' && seq[i] != 'E' && seq[i] != 'G' &&
         seq[i] != 'H' && seq[i] != 'I' && seq[i] != 'L' && seq[i] != 'K' &&
         seq[i] != 'M' && seq[i] != 'F' && seq[i] != 'P' && seq[i] != 'S' &&
         seq[i] != 'S' && seq[i] != 'T' && seq[i] != 'W' && seq[i] != 'Y' &&
         seq[i] != 'V' && seq[i] != 'B' && seq[i] != 'J' && seq[i] != 'Z' &&
         seq[i] != 'X' && seq[i] != '*')  {
        seq[i] = 'X';
      }
    }
    return;
  }
  
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
		//if ( compressed ) 
		if ( fio::getFileExtension( std::string(filename) ) == "gz" ) 
			in.push(boost::iostreams::gzip_decompressor());
		in.push(fstrm);
		std::istream istrm(&in);


/* 		std::istream* istrm = &fstrm; */
/* 		if ( compressed ) { */
/* 			boost::iostreams::filtering_streambuf<boost::iostreams::input> in; */
/* 			in.push(boost::iostreams::gzip_decompressor()); */
/* 			in.push(fstrm); */
/* 			std::istream bis(&in); */
/* 			istrm = &bis; */
/* 		}  */

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
		  reAssignSpecialChar(fasta_seq);
			update(tags, seqs, fasta_tag, fasta_seq, count, mode);
			count++;
		}
		fstrm.close();
	}
  

	inline std::set<int> loadIndices( const char *filename ) 
	{
		std::list<int> qlist;
		
		std::fstream fin;
		fio::openFile(fin, filename, std::ios::in);
	
		int item;
		fin >> item;
		while ( !fin.eof() ) {
			qlist.push_back(item);
			fin >> item;
		}
		fin.close();
    
		std::list<int>::iterator it = qlist.begin();
		std::list<int>::iterator jt = qlist.end();
		return std::set<int>(it, jt);
	}

	inline std::set<std::string> loadTags( const char *filename ) 
    {
		std::list<std::string> qlist;
    
		std::fstream fin;
		fio::openFile(fin, filename, std::ios::in);
	
		std::string item;
		getline(fin, item);
		while ( !fin.eof() ) {
			qlist.push_back(item);
			getline(fin, item);
		}
		fin.close();

		std::list<std::string>::iterator it = qlist.begin();
		std::list<std::string>::iterator jt = qlist.end();
		return std::set<std::string>(it, jt);
	}

						   

/* 	// FASTA file sequence count */
/* 	int getSequenceCount( std::string input_file )  */
/* 	{ */
/* 		std::string file = input_file; */
/* 		const char *filename = input_file.c_str(); */
/* 		std::fstream fstrm(filename, std::ios_base::in | std::ios_base::binary); */
	
/* 		if ( !fstrm ) { */
/* 			std::cerr << "Can't open " << filename << "\n"; */
/* 			exit (1); */
/* 		} */
    
/* 		int count = 0; */
    
/* 		std::string line; */
/* 		getline(fstrm, line);		 */
/* 		while ( !fstrm.eof() ) { */
/* 			if ( line[0] == '>' ) count++; */
/* 			getline(fstrm, line);		 */
/* 		} */
	
/* 		return count; */
/* 	} */


/* 	int getGzipSequenceCount( std::string input_file ) */
/* 	{ */
/* 		std::string file = input_file; */
/* 		const char *filename = input_file.c_str(); */
/* 		std::fstream fstrm(filename, std::ios_base::in | std::ios_base::binary); */
	
/* 		if ( !fstrm ) { */
/* 			std::cerr << "Can't open " << filename << "\n"; */
/* 			exit (1); */
/* 		} */
    
/* 		boost::iostreams::filtering_streambuf<boost::iostreams::input> in; */
/* 		in.push(boost::iostreams::gzip_decompressor()); */
/* 		in.push(fstrm); */

/* 		int count = 0; */
    
/* 		std::istream istrm(&in); */
/* 		std::string line; */
/* 		while ( std::getline(istrm, line) ) { */
/* 			if ( line[0] == '>' ) count++; */
/* 		} */
/* 		return count; */
/* 	} */

	inline int getSequenceCount( const char *filename ) //, bool compressed ) 
	{
		std::fstream fstrm(filename, std::ios_base::in | std::ios_base::binary);
	
		if ( !fstrm ) {
			std::cerr << "Can't open " << filename << "\n";
			exit (1);
		}

		boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
		//if ( compressed ) 
		//if ( io::isGzipFile( filename) ) 
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

	// Sequence count of multiple FASTA files
	inline int totalSequenceCount( const std::vector<std::string> &input_files )//, bool zflag ) 
	{
		int total = 0;
		for ( size_t i = 0; i < input_files.size(); i++ ) {
			//if ( !zflag ) total += getSequenceCount( input_files[i] );
			//else total += getGzipSequenceCount( input_files[i] );
			total += getSequenceCount( input_files[i].c_str() );//, zflag );
			//std::cerr << i << "\t" << total << "\n";
		}
		return total;
	}


	inline void loadSequenceFile( std::string input_file, 
								  char **tags,
								  char **seqs,
								  int &count,
								  int mode,
								  bool debug = false) //, 
	//bool zflag )
	{
		//double time1 = mytime();
		std::string file = input_file;
		//if ( !zflag ) readSequences(file.c_str(), tags, seqs, count, mode );
		//else readGzipSequences(file.c_str(), tags, seqs, count, mode );
		readSequences(file.c_str(), tags, seqs, count, mode, debug);//, zflag );
		//std::cerr << "\telapsed: " << mytime()-time1 << " sec\n";
	}


	inline void loadSequences( const std::vector<std::string> &input_files, 
							   char **tags,
							   char **seqs,
							   int mode,
							   bool debug = false)//,
	//bool zflag )
		
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

	inline void purge( char **mems, int len ) 
	{
		for ( int i = 0; i < len; i++ )
			delete[] mems[i];
		delete[] mems;
		mems = NULL;
	}
	
	inline void extractByIndex( const char *qfile, char **seqs, char **tags, int size  )
	{
		std::set<int> query = loadIndices(qfile);
		for ( int i = 0; i < size; i++ )
			if ( query.find(i) != query.end() )
				std::cout << ">" << tags[i] << "\n" << seqs[i] << "\n";
	}

	inline void extractByTag( const char *qfile, char **seqs, char **tags, int size )
	{
		std::set<std::string> query = loadTags(qfile);
		for ( int i = 0; i < size; i++ )
			if ( query.find(std::string(tags[i])) != query.end() )
				std::cout << ">" << tags[i] << "\n" << seqs[i] << "\n";
	}

	/* 	bool checkLength(TagType &tag, SeqType &seq, int k); */
	/* 	bool validSequence(int index , SeqType &seq); */
	/* 	void update( char **tags, char **seqs, TagType &fasta_tag, SeqType &fasta_seq, int &index, int mode); */
	/* 	void readSequences(  const char* filename, char **tags, char **seqs, int &index, int mode); */
	/* 	void readGzipSequences(  const char* filename, char **tags, char **seqs, int &index, int mode); */
	/* 	std::set<int> loadIndices( const char *filename ); */
	/* 	std::set<std::string> loadTags( const char *filename ); */
	/* 	int getSequenceCount( std::string input_file ); */
	/* 	int getGzipSequenceCount( std::string input_file  ); */
	/* 	int totalSequenceCount( const std::vector<string> &input_files, bool zflag ); */
	/* 	void loadSequenceFile( std::string input_file, char **tags, char **seqs, int &count, int mode, bool zflag); */
	/* 	void loadSequences( const std::vector<string> &input_files, char **tags, char **seqs, int mode, bool zflag); */
	/* 	void purge( char **mems, int len );		 */
	/* 	void extractByIndex( const char *qfile, char **seqs, char **tags, int size  ); */
	/* 	void extractByTag( const char *qfile, char **seqs, char **tags, int size ); */
	/* 	void extractSequences( const char *qfile, std::vector<std::string> &seq_files, bool byid, int field); */
}

#endif

