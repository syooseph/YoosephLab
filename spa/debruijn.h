/** 
 * \file debruijn.h
 * Debruijn graph handler.
 * This file includes definition of debruijn graph and declaration of 
 * utility function of the graph object.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \pre       Boost graph library must be installed.
 * \bug       Not known.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef DEBRUJIN_H
#define	DEBRUJIN_H

#include <iostream>
#include <tr1/unordered_map>
#include <boost/graph/adjacency_list.hpp>
#include "timer.h"
#include "smath.h"

/*====================================*/
/* definition of node/edge properties */
/*====================================*/
/** Edge weight */
typedef float EdgeWeight; 

/** Empty vertex property */
struct VertexProperties 
{
};

/** Edge property */
struct EdgeProperties 
{
	EdgeWeight weight; ///< Count of supporting reads.
};


/** Definition of graph */
typedef boost::adjacency_list 
	< 
	boost::setS,     // disallow parallel edges
	boost::listS,    // vertex container
	boost::directedS,
	VertexProperties,
	EdgeProperties
	> DeBruijnGraph;


#endif	/* DEBRUJIN_H */
