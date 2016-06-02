/**
 * \file      seed.h
 * \brief     Path seeds
 * \details   Path seeds are maintained by two hash tables. 
 *            First hash table stores vertex and bin ID.
 *            Second hash table stroes bin ID and a set of vertices.
 * \author    Youngik Yang
 * \version   0.2
 * \date    2013
 * \copyright J. Craig Venter Institute.
 */

#ifndef __SEED_H__
#define __SEED_H__


#include <boost/unordered_map.hpp>
#include <unordered_map>
#include <vector>
#include "graph.h"
#include "param.h"

typedef double BinType;
typedef boost::unordered_map<Vertex, BinType> VertexToBinMap;
// Do not use tr1/unordered_map for this
typedef boost::unordered_map<Vertex, bool> NodeFlagMap;
typedef std::map<BinType, NodeFlagMap> BinMemberMap;
typedef std::pair<BinType, bool> BinFlag;

/**
 * \class Seed
 * \brief Path Seeds for path extraction,
 */
class Seed
{
 private:
	VertexToBinMap node2bins; ///< Hash (key:vertex, value:bin)
	BinMemberMap   bin2nodes; ///< Hash (key:bin, value: a set of vertices)

 public:
	/** Insert a vertex to seed bin */
	void insert(Vertex n, BinType v);

	/** Update vertex to given bin value */ 
	void update(Vertex n, BinType v);

	/** Erase a vertex from seed pool */
	void erase(Vertex n);

	/** Get one seed from pool */
	Vertex getSeed();

	/** Pop one seed from seed pool */
	Vertex pop();

	/** Check existence of seed/bin pair */
	bool has(Vertex, BinType&);

	/** Get size of seed pool */
	size_t getSize();

	/** Clear seed pool */
	void clear();

	/** Get bin of given vertex */
	BinType getBin(Vertex node);

	/** Get a number of bins */
	size_t getBinSize();
};

#endif
