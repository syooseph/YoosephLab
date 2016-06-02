/** 
 * \file      InvertedIndex.h
 * \brief     Inverted index.
 * \details   This class creates an inverted index between kmer and reads.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \bug       Not known.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __INVERTED_INDEX_H__
#define __INVERTED_INDEX_H__


#include <list>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "file.h"
#include "kmer.h"

typedef std::tr1::unordered_map<KmerId, Read*> InvertedIndexMap;
typedef std::set<KmerId> KmerSet;

/** Data entry for inverted index */
class Record 
{
 private:
	ReadId rid;
	KmerSet kmers;
 public:
	Record( ReadId &read, KmerSet &kids ) { rid = read; kmers = kids; }
	ReadId getRid() { return rid; }
	KmerSet getKmers() { return kmers; }
};

/** 
 * Inverted Index.
 * Inverted index between a kmer to a set of read IDs.
 */
class InvertedIndex 
{
 private:
	InvertedIndexMap index;

 public:

	/** Getter: inverted index */
	InvertedIndexMap getInvertedIndex()
	{
		return index;
	}

	/** Reclaim memory */
	void clear() 
	{
		for ( InvertedIndexMap::iterator it = index.begin(); it != index.end(); ++it ) 
			delete it->second;
		index.clear();
	}

	/** 
	 * Add entry.
	 * For each kmer in a given read, add the read to index.
	 * \param Record        Read ID and kmers of the sequence.
	 * \param CoverageMap	Count of reads for each kmer.
	 * \pre   CoverageMap must be valid. For each kmer entry, the 
	 *        count of reads must correctly cacluated beforehand.
	 */
	void add(Record &r, CoverageMap &coverages) 
	{
		ReadId rid    = r.getRid();
		KmerSet kmers = r.getKmers();

		std::set<KmerId>::iterator it;
		for ( it = kmers.begin(); it != kmers.end(); ++it ) {
			//============================================
			// if an index for a kmer does not exist, 
			// reserve the size of read count of the kmer.
			//============================================
			if ( index.find(*it) == index.end() ) 
				index[*it] = new Read(coverages[*it]);
			index[*it]->add(rid); // add new read ID.
		}
	}

	/** Does exist the key ? */
	bool has( KmerId &key ) 
	{
		if ( index.find(key) == index.end() ) return 0;
		return 1;
	}

	/** Erase the key/value from index */
	void erase( KmerId &key ) 
	{
		if ( ! has(key) ) {
			std::cerr << "Key does not exist\n";  
			return;
		}
		delete index[key];
		index.erase(key);
	}

	/** Set new reads to an existing index */
	void update(KmerId &key, ReadId* rids, size_t size )
	{
		erase(key);

		index[key] = new Read(size);
		for ( size_t i = 0; i < size; i++ ) 
			index[key]->add(rids[i]);
	}
	
	/** Return values by key */
	Read* getValue( KmerId &key ) 
	{
		if ( ! has(key) ) 
			return NULL;
		return index[key];
	}

	/** Index size */
	int getSize() 
	{
		return index.size();
	}

	/** Dump index to a binary file */
	void write( const char *file ) 
	{
		std::fstream out;
		io::openFile( out, file, std::ios::out | std::ios::binary );

		InvertedIndexMap::iterator it;
		CoverageType size;
		Read *reads;
		for ( it = index.begin(); it != index.end(); ++it ) {
			out.write((char*)&(*it).first, sizeof(KmerId));

			reads = it->second;
			size = reads->size;
			out.write((char*)&size, sizeof(CoverageType));
			for ( int i = 0; i < (int)size; i++ ) 
				out.write((char*)&(reads->rid[i]), sizeof(ReadId));
		}
		out.close();
	}

	/** Read an index entry */
	//void __readRecord( std::fstream &in, KmerId &kid  )
	void __readRecord( std::istream &in, KmerId &kid  )
	{
		CoverageType count = 0;
		ReadId rid;
		in.read((char*)&kid, sizeof(KmerId));
		if ( in.eof() ) return; // enforce this to avoid appending junk
		in.read((char*)&count, sizeof(CoverageType));
		if ( in.eof() ) return; // enforce this to avoid appending junk
		index[kid] = new Read(count);
		for ( int i = 0; i < (int)count; i++ ) { 
			in.read((char*)&rid, sizeof(ReadId));   
			index[kid]->add(rid);
		}
	}

	//===========================
	// load saved index from disk
	//=========================== 
	/** Load index from a binary file */
	void load( const char *file ) {
		std::fstream fstrm;
		io::openFile( fstrm, file, std::ios::in | std::ios::binary );
		fstrm.seekg(0);

		std::string fname = file;
		boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
		if ( fname.substr( fname.find_last_of(".") + 1 ) == "gz" )
			in.push(boost::iostreams::gzip_decompressor());
		in.push(fstrm);
		std::istream iin(&in);
		
		KmerId kid;
		while ( !iin.eof() ) {
			__readRecord( iin, kid );
		}			
	}
};

#endif
