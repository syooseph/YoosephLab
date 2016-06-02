/** 
 * \file      biostr.h
 * \brief     Protein sequence handler
 * \details   This namespace includes conversion between kmers and 
 *            protein seuqence and gap removals.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \bug       None.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __BIOSTR_H__
#define __BIOSTR_H__

//#include "relation.h"
#include "graph.h"

/**
 * Conversion between kmers and sequences
 */
namespace biostr 
{
	/**
	 * \brief   Convert kmers to a protein sequence.
	 * \details 
	 * \param   kmers An array of kmers
	 * \param   nkmer The size of kmer array
	 * \param   k     A size of kmer
	 * \return  A string of protein sequence
	 */
	inline std::string getSequenceString( const KmerId* kmers, int nkmer, int k )
	{
		std::string seq;
		for ( int i = 0; i < nkmer; ++i ) {
			KmerId kmer = kmers[i];
			if ( i == 0 ) 
				seq = alpha::IntegerToAminoAcid(kmer, k);
			else 
				seq += alpha::getLastAminoAcid(kmer, k);
		}
		return seq;
	}
	

	/**
	 * \brief   Convert a node to a kmer length string.
	 * \param   vertex      Graph node
	 * \param   vertex_map  Map between graph node to kmer
	 * \param   k           A size of kmer
	 * \return  A kmer length string
	 */
	inline std::string getSequenceString( Vertex vertex, VertexToKmerMap &vertex_map, int k )
	{
		std::string seq;
		KmerId kmer = vertex_map[vertex];
		return alpha::IntegerToAminoAcid(kmer, k);
	}
	
	
	/**
	 * \brief   Convert nodes to amino acids.
	 * \param   NodeArray	A set of nodes
	 * \param   vertex_map  Map between graph node to kmer
	 * \param   k           A size of kmer
	 * \return  A string of amino acids.
	 */
	inline std::string getSequenceString( const NodeArray &vertices, VertexToKmerMap &vertex_map, int k )
	{
		std::string seq;
		for ( size_t i = 0; i < vertices.size(); ++i ) {
			KmerId kmer = vertex_map[vertices[i]];
			if ( i == 0 ) 
				seq = alpha::IntegerToAminoAcid(kmer, k);
			else 
				seq += alpha::getLastAminoAcid(kmer, k);
		}
		return seq;
	}
	
	
	/**
	 * \brief   Convert path to amino acids.
	 * \param   PathType	A path that consists of nodes.
	 * \param   vertex_map  Map between graph node to kmer
	 * \param   k           A size of kmer
	 * \return  A string of amino acids.
	 */
	inline std::string getSequenceString( const PathType &path, VertexToKmerMap &vertex_map, int k )
	{
		NodeArray vertices = NodeArray( path.begin(), path.end() );
		return getSequenceString( vertices, vertex_map, k );
	}
	
	/**
	 * \brief   Get a sequence length.
	 * \param   PathType	A path that consists of nodes.
	 * \param   vertex_map  Map between graph node to kmer
	 * \param   k           A size of kmer
	 * \return  A length of a string of amino acids.
	 */
	inline size_t getSequenceLength( const PathType &path, VertexToKmerMap &vertex_map, int k ) 
	{
		return path.size()+k-1;
	}
	
	/**
	 * \brief   Get a sequence length.
	 * \param   NodeArray	A set of nodes.
	 * \param   vertex_map  Map between graph node to kmer
	 * \param   k           A size of kmer
	 * \return  A length of a string of amino acids.
	 */
	inline size_t getSequenceLength( const NodeArray &vertices, VertexToKmerMap &vertex_map, int k ) 
	{
		return vertices.size()+k-1;
	}


	/**
	 * \brief   Remove gaps.
	 * \param   String a sequence string
	 * \return  A string with no gaps.
	 */
	inline std::string stripGap( std::string seq )
	{
		for ( int i = (int)seq.size()-1; i >= 0; i-- )
			if ( seq[i] == '-' ) seq.erase(i, 1);
		return seq;
	}

	/**
	 * \brief   Convert a protein sequence to kmers.
	 * \details 
	 * \param   str  A protein sequence
	 * \param   k    A size of kmer
	 * \return  A vector of kmers
	 */
	inline KmerArray getKmers( std::string str, int k )
	{
		KmerList kids;
		for (int i = 0; i <= (int)str.size() - k; i++) {
			KmerType kmer = str.substr(i,k);
			kids.push_back( alpha::AminoAcidToInteger<KmerId>(kmer) );
		}
		return KmerArray( kids.begin(), kids.end() );
	}
}

#endif
