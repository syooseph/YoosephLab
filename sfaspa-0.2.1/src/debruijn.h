/** 
 * \file      debruijn.h
 * \brief     Debruijn graph handler.
 * \details   This file includes definition of debruijn graph and declaration of 
 *            utility function of the graph object.
 * \author    Youngik Yang
 * \version   0.2
 * \date      2010-2011
 * \date      modifed on Tue 2013-12-17 06:21:01 PM
 * \pre       Boost graph library must be installed.
 * \bug       Not known.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef DEBRUJIN_H
#define	DEBRUJIN_H

#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include "timer.h"
#include "smath.h"

/*====================================*/
/* definition of node/edge properties */
/*====================================*/

/** Edge weight */
typedef float EdgeWeight; 

/** 
 * \brief vertex property for deBrujin Graph
*/
struct VertexProperties 
{
};

/** 
 * \brief Edge property deBrujin Graph
*/
struct EdgeProperties 
{
	EdgeWeight weight; ///< Count of supporting reads.
};


/*====================================*/
/* definition of node/edge properties */
/*====================================*/
typedef boost::adjacency_list 
	< 
	boost::setS,      // disallow parallel edges
	boost::listS,     // vertex container
	boost::directedS, 
	VertexProperties,
	EdgeProperties
	> DeBruijnGraph;


#endif	/* DEBRUJIN_H */
