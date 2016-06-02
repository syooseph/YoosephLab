#include "graph.h"

void graph::build( DeBruijnGraph &graph,
				   CoverageMap &kmer_coverage,
				   VertexToKmerMap &vertex_map,
                   int kmer_size,
                   std::string &graph_input,
                   bool verbose )
{
    KmerToVertexMap kmer_map;   // kmer-id to vertex (address) mapping
	
	graph::buildGraph( graph, kmer_coverage, kmer_map, kmer_size, graph_input, verbose );
    vertex_map = revertMap( kmer_map );
}


void graph::buildGraph( DeBruijnGraph &graph, 
						CoverageMap &kmer_coverage,
						KmerToVertexMap &vertex_map,
                        int kmer_size, 
                        std::string &graph_input,
                        bool verbose )
{
    std::fstream fstrm;
    fio::openFile( fstrm, graph_input.c_str(), std::ios::in | std::ios::binary );

    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    if ( graph_input.substr( graph_input.find_last_of(".") + 1 ) == "gz" )
        in.push(boost::iostreams::gzip_decompressor());
    in.push(fstrm);
    std::istream gin(&in);
    
    KmerId kmer_id;
    CoverageType coverage;
    CoverageType left[alpha::COUNT_AA];

    gin.read((char*)&kmer_id, sizeof(KmerId));
    gin.read((char*)&coverage, sizeof(CoverageType));
    for ( int i = 0; i < alpha::COUNT_AA; i++ ) 
        gin.read((char*)&(left[i]), sizeof(CoverageType));

    double t1 = mytime();
    int count = 0;    
    while ( !gin.eof() ) {
        updateTopology( graph, vertex_map, kmer_id, left, kmer_size );
        kmer_coverage[kmer_id] = coverage;

        gin.read((char*)&kmer_id, sizeof(KmerId));
        gin.read((char*)&coverage, sizeof(CoverageType));
        for ( int i = 0; i < alpha::COUNT_AA; i++ ) 
            gin.read((char*)&(left[i]), sizeof(CoverageType));

        ++count;
        if ( count % 10000000 == 0 ) {
            if ( verbose ) std::cout << "\t" << count << " lines processed (" 
                                     << num_vertices(graph) << " vertices, " 
                                     << num_edges(graph) << " edges, " << mytime()-t1 << " sec)\n";
            t1 = mytime();
        }
    }
    fstrm.close();
}

/**
 * \pre Kmer ID must be valid.
 */
void graph::updateTopology( DeBruijnGraph &graph, 
                     KmerToVertexMap &vertex_map,
                     KmerId kmer_id,
                     CoverageType left[],
                            int kmer_size )
{
    Vertex target, source;

    target = vertex_map[kmer_id];
    if ( ! target ) {
        target = add_vertex(graph);
        vertex_map[kmer_id] = target;
    } 

    for ( int i = 0; i < alpha::COUNT_AA; i++ ) {
        if ( left[i] == 0 ) continue;
        KmerType left_kmer = alpha::AminoAcid[i+1] + alpha::IntegerToAminoAcid(kmer_id, kmer_size).substr(0,kmer_size-1);
        KmerId left_id = alpha::AminoAcidToInteger<KmerId>(left_kmer);

        source = vertex_map[left_id];
        if ( ! source ) {
            source = add_vertex(graph);
            vertex_map[left_id] = source;
        } 
        std::pair<Edge, bool> edge = add_edge(source, target, graph);
        graph[edge.first].weight = left[i];
    }
}


bool graph::formCycle( NodeArray &path, 
                           NodeArray &nodes )
{
    for ( size_t i = 0; i < nodes.size(); ++i ) 
        for ( size_t j = 0; j < path.size(); ++j ) 
            if ( nodes[i] == path[j] ) 
                return 1;
    return 0;
}

bool graph::formCycle( NodeArray &path, 
                           NodeType  target )
{
    for ( size_t i = 0; i < path.size(); ++i ) 
        if ( target == path[i] ) 
            return 1;
    return 0;
}

bool graph::formCycle( NodeList &path, 
                           NodeType target )
{
    for ( NodeList::iterator it = path.begin(); it != path.end(); ++it )
        if ( target == *it ) 
            return 1;
    return 0;
}

bool graph::hasNode( PathType &path,
                     NodeType target )
{
    return formCycle(path, target);
}


void graph::dropAncestor( AdjacencyMap &ancestors, 
                              Vertex key,
                              Vertex val)
{
    NodeArray preds = ancestors[key];
    NodeArray::iterator it;
    for ( it = preds.begin(); it != preds.end(); ) {
        if ( *it == val ) it = preds.erase(it);
        else ++it;
    }
    ancestors[key] = preds;
}

void graph::dropAncestors( AdjacencyMap &ancestors, 
                               PathType &minor_path,
                               PathType &major_path )
{
    NodeArray major = std::vector<NodeType>(major_path.begin(), major_path.end());
    NodeArray minor = std::vector<NodeType>(minor_path.begin(), minor_path.end());

    for ( size_t i = 0; i < minor.size(); i++ ) {
        if ( minor[i] == major[i] ) continue;
        ancestors.erase( ancestors.find(minor[i]) );
    }
}


NodeFlags graph::trimGraph( DeBruijnGraph &graph,
                            CoverageMap &kmer_coverage,
                            VertexToKmerMap &vertex_map,
                            int min_depth,
                            bool verbose )
{
    if ( verbose ) {
        std::cout << "... Trimming Vertices\n";
    }
    NodeList trimmed = graph::trimVertices(graph, kmer_coverage, vertex_map, min_depth);
    
    if ( verbose ) std::cout << "... Trimming Edges\n";
    graph::trimEdges(graph, min_depth);

    if ( verbose ) std::cout << "... Trimming Islands\n";
    NodeList dropped = graph::dropIslands(graph, vertex_map);

    NodeFlags flags;
    for ( auto node : trimmed )
        flags[node] = true;
    for ( auto node : dropped )
        flags[node] = true;

    return flags;

    // for ( NodeList::iterator it = dropped.begin(); it != dropped.end(); ++it )
    //     trimmed.push_back(*it);

    // return NodeSet( trimmed.begin(), trimmed.end() );
}


NodeList graph::trimVertices( DeBruijnGraph &graph,
                              CoverageMap &kmer_coverage,
                              VertexToKmerMap &vertex_map,
                              size_t min_coverage )
{
    AdjacencyMap ancestors = predecessorMap(graph);


    std::tr1::unordered_map<Vertex, bool> drop_map;

    Vertex_iter vc, ve, next;
    boost::tie(vc, ve) = vertices(graph);
    for ( next = vc; vc != ve; vc = next ) {
        ++next;

        if ( kmer_coverage[vertex_map[*vc]] >= min_coverage ) continue;

        NodeArray preds = ancestors[*vc];
        for ( size_t i = 0; i < preds.size(); i++ ) {
            if ( drop_map.find(preds[i]) != drop_map.end() ) continue;
            if ( edge(preds[i], *vc, graph).second)
                remove_edge(preds[i], *vc, graph);
        }

        NodeArray succs = successors(graph, *vc);
        for ( size_t i = 0; i < succs.size(); i++ )
            if ( edge(*vc, succs[i], graph).second )
                remove_edge(*vc, succs[i], graph);

        drop_map.insert( std::pair<Vertex, bool>(*vc, true) );
        remove_vertex(*vc, graph);
    }

    NodeList deleted;
    std::tr1::unordered_map<Vertex, bool>::iterator it;
    for ( it = drop_map.begin(); it != drop_map.end(); ++it )
        deleted.push_back(it->first);

    return deleted;
}

void graph::trimEdges( DeBruijnGraph &graph,
                           size_t min_coverage )
{
	Edge_iter ec, ee, next;
    boost::tie(ec, ee) = edges(graph);
    for ( next = ec; ec != ee; ec = next ) {
        ++next;
        if ( graph[*ec].weight < min_coverage )
            remove_edge( *ec, graph );
    }
}

NodeList graph::dropIslands( DeBruijnGraph &graph, 
                                 VertexToKmerMap &vertex_map
                                 ) 
{
    NodeList dropped;
    AdjacencyMap ancestors = predecessorMap(graph);

    Vertex_iter vc, ve, next;
    boost::tie(vc, ve) = vertices(graph);
    for ( next = vc; vc != ve; vc = next ) {
        ++next;
        if ( ancestors[*vc].size() == 0 && successors(graph, *vc).size() == 0 ) {
            dropped.push_back(*vc);
            remove_vertex(*vc, graph);
        }
    }
    return dropped;
}


//==============================================================================
// Fri 2011-04-01 03:52:10 PM
// BFS graph of given depth
// 2. Path update
// 3. Find a bubble
//==============================================================================
void graph::BFS( PathQueue &queue,
                     DeBruijnGraph &graph,
                     int max_depth )
{
    //-----------------------------------
    // Get the top sequence in the queue
    // Get the first node in the sequence
    //-----------------------------------
    
    PathType top = queue.front();     
    if ( (int)top.size() >= max_depth ) return;

    queue.pop();
    Vertex last  = top.back();

    if ( last == NULL ) {
        PathType ntop = top;
        ntop.push_back(NULL);
        queue.push(ntop);
    }
    
    else {
        NodeArray children = successors( graph, last );
        if ( children.size() == 0 ) {
            PathType ntop = top;
            ntop.push_back(NULL);
            queue.push(ntop);
        }

        for ( size_t i = 0; i < children.size(); ++i ) {
            PathType ntop = top;
            ntop.push_back(children[i]);
            queue.push(ntop);
        }
    }
    BFS(queue, graph, max_depth );
}

void graph::reverseBFS( PathQueue &queue,
                            DeBruijnGraph &graph,
                            AdjacencyMap &ancestors,
                            int max_depth )
{
    //-----------------------------------
    // Get the top sequence in the queue
    // Get the first node in the sequence
    //-----------------------------------
    
    PathType top = queue.front();     
    if ( (int)top.size() >= max_depth ) return;

    queue.pop();
    Vertex front  = top.front();

    if ( front == NULL ) {
        PathType ntop = top;
        ntop.insert(ntop.begin(), NULL);
        queue.push(ntop);
    }
    
    else {
        NodeArray predecessors = ancestors[front];
        if ( predecessors.size() == 0 ) {
            PathType ntop = top;
            ntop.insert(ntop.begin(), NULL);
            queue.push(ntop);
        }

        for ( size_t i = 0; i < predecessors.size(); ++i ) {
            PathType ntop = top;
            ntop.insert(ntop.begin(), predecessors[i]);
            queue.push(ntop);
        }
    }
    reverseBFS(queue, graph, ancestors, max_depth );
}


/*===================*/
/* Utility functions */
/*===================*/

/**
 * In-degree map.
 * Boost graph library does not supports indegree map for a directed 
 * graph. This function creates indegree map.
 * \param g DeBruijnGraph object.
 * \return IndegreeMap, count of indegrees of all the node in the graph.
 */
InDegreeMap graph::indegrees(DeBruijnGraph &g) 
{
	InDegreeMap in_degree;
	Vertex_iter i, iend;
	OutEdge_iter j, jend;
	for (boost::tie(i, iend) = vertices(g); i != iend; ++i)
		for (boost::tie(j, jend) = out_edges(*i, g); j != jend; ++j)
			in_degree[target(*j, g)] += 1;
	
	return in_degree;
}

/**
 * A set of predecessors.
 * \param g	DeBruijn Graph
 * \param v	vertex of interest 
 * \return predecessors of a vertex
 */ 
NodeArray graph::predecessors(DeBruijnGraph &g, Vertex v) 
{	
	NodeArray parents;

	Edge_iter ec, ee;
	for ( boost::tie(ec,ee) = edges(g); ec != ee; ++ec ) 
		if ( target(*ec, g) == v ) 
			parents.push_back(source(*ec,g));
	
	return parents;
}

/**
 * Predecessor map.
 * \param  g		DeBruijn Graph
 * \return Mapping for predecessor nodes of all nodes.
 */ 
AdjacencyMap graph::predecessorMap(DeBruijnGraph &g)  
{
	AdjacencyMap ancestors;
	Edge_iter ec, ee;
	for ( boost::tie(ec,ee) = edges(g); ec != ee; ++ec ) 
        ancestors[target(*ec,g)].push_back(source(*ec,g));

	return ancestors;
}


/**
 * Child nodes.
 * \param g	DeBruijn Graph
 * \param v	vertex of interest 
 * \return successors of a vertex.
 */ 
NodeArray graph::successors(DeBruijnGraph &g, Vertex v) 
{
	NodeArray children;

	OutEdge_iter ec, ee;
	for ( boost::tie(ec,ee) = out_edges(v, g); ec != ee; ++ec ) 
		children.push_back(target(*ec,g));
	
	return children;
}

/**
 * Successor map.
 * \param  g		DeBruijn Graph
 * \return Mapping for successor nodes of all nodes.
 */ 
AdjacencyMap graph::successorMap(DeBruijnGraph &g) 
{
	AdjacencyMap children;
	Edge_iter ec, ee;
	for ( boost::tie(ec,ee) = edges(g); ec != ee; ++ec ) {
        children[source(*ec,g)].push_back(target(*ec,g));
	}	
	return children;
}

/**
 * Source nodes
 * \param  g		DeBruijn Graph
 * \return A list of all source nodes
 */ 
NodeList graph::getSources(DeBruijnGraph &g) 
{
	Vertex_iter vc, ve;
	InDegreeMap in_degree = indegrees(g);

	NodeList sources;
	for ( boost::tie(vc, ve) = vertices(g); vc != ve; ++vc )
		if ( out_degree(*vc, g) >  0 && in_degree[*vc] == 0 ) 
			sources.push_back(*vc);

	return sources;
}

/**
 * Sink nodes
 * \param  g		DeBruijn Graph
 * \return A list of all sink nodes
 */ 
NodeList graph::getSinks(DeBruijnGraph &g) 
{
	Vertex_iter vc, ve;
	InDegreeMap in_degree = indegrees(g);

	NodeList sinks;
	for ( boost::tie(vc, ve) = vertices(g); vc != ve; ++vc )
		if ( out_degree(*vc, g) == 0 && in_degree[*vc] >  0 ) 
			sinks.push_back(*vc);

	return sinks;
}

/**
 * Junction nodes.
 * This function extracts junction nodes from the graph, where 
 * the junction node is defined as indegree(v)*outdegree(v) > 1.
 * \param  g		DeBruijn Graph
 * \return A list of all junction nodes
 */ 
NodeList graph::getJunctions(DeBruijnGraph &g) 
{
	Vertex_iter vc, ve;
	InDegreeMap in_degree = indegrees(g);

	NodeList junctions;
	for ( boost::tie(vc, ve) = vertices(g); vc != ve; ++vc )
		if ( out_degree(*vc, g) * in_degree[*vc] >  1 ) 
			junctions.push_back(*vc);

	return junctions;
}


/**
 * Find a vertex.
 * \param graph  DeBruijn Graph
 * \param vertex vertex of interest 
 * \return True/False.
 */ 
bool graph::hasVertex(DeBruijnGraph &graph, Vertex vertex) 
{
	Vertex_iter vc, ve;
	for ( boost::tie(vc,ve) = vertices(graph); vc != ve; ++vc ) {
		if ( *vc == vertex ) 
			return 1;
	}
	return 0;
}


bool graph::isSource(Vertex v,
					 //DeBruijnGraph &g, 
					 AdjacencyMap &ancestors) 
{
    if ( ancestors[v].size() == 0 )
        return true;
    return false;
}

bool graph::isSink(Vertex v, 
				   DeBruijnGraph &g )
{
    if ( out_degree(v,g) == 0 )
        return true;
    return false;
}


/*===================*/
/* Summary functions */
/*===================*/
/**
 * Display graph summary
 */
void graph::graphSummary( DeBruijnGraph & graph, std::ostream &out ) 
{
    int n = num_vertices(graph);
    
    std::vector<double> degrees;
    degrees.reserve(n);

    int source, sink, one2one, one2many, many2one, many2many, island;
    source = sink = one2one = one2many = many2one = many2many = island = 0;

    AdjacencyMap ancestors = predecessorMap(graph);
    Vertex_iter vc, ve;
    for ( boost::tie(vc, ve) = vertices(graph); vc != ve; ++vc ) {
        NodeArray preds = ancestors[*vc];
        NodeArray succs = successors(graph, *vc);
        
        degrees.push_back( double( preds.size()+succs.size() ) );
        
        if ( preds.size() == 0 && succs.size() == 0 ) island++;
        else if ( preds.size() == 0 && succs.size() >  0 ) source++;
        else if ( preds.size() >  0 && succs.size() == 0 ) sink++;
        else if ( preds.size() == 1 && succs.size() == 1 ) one2one++;
        else if ( preds.size() == 1 && succs.size() >  1 ) one2many++;
        else if ( preds.size() >  1 && succs.size() == 1 ) many2one++;
        else if ( preds.size() >  1 && succs.size() >  1 ) many2many++;
    }

 	out << "\nGraph Summary\n";
	out << "# nodes  :" << num_vertices(graph) << std::endl;
	out << "# edges  :" << num_edges(graph) << std::endl;
    //out << "\t# islands:" << island << "\n";
    out << "# sources:" << source << "\n";
    out << "# sinks  :" << sink << "\n";
    // out << "# 1 to 1 :" << one2one << "\n";
    // out << "# 1 to M :" << one2many << "\n";
    // out << "# M to 1 :" << many2one << "\n";
    // out << "# M to M :" << many2many << "\n";

    out << "Degrees: ( ";
    out << "max:" << math::max( &degrees[0], n ) << " ";
    out << "min:" << math::min( &degrees[0], n ) << " ";
    out << "avg:" << math::mean( &degrees[0], n ) << " ";
    out << "med:" << math::median( &degrees[0], n, false ) << " )\n";
}

void graph::coverageSummary( DeBruijnGraph &graph, 
                             CoverageMap &kmer_coverage,
                             VertexToKmerMap &vertex_map,
                             std::ostream &out ) 
{
    int n = num_vertices(graph);
    
    std::vector<double> coverages;
    coverages.reserve(n);

    Vertex_iter vc, ve;
    for ( boost::tie(vc, ve) = vertices(graph); vc != ve; ++vc ) {
        KmerId kid = vertex_map[*vc];
        coverages.push_back( (double)kmer_coverage[kid] );
    }

    out << "\nKmer coverage summary\n";
    out << "max:" << math::max( &coverages[0], n ) << " ";
    out << "min:" << math::min( &coverages[0], n ) << " ";
    out << "avg:" << math::mean( &coverages[0], n ) << " ";
    out << "med:" << math::median( &coverages[0], n, false ) << "\n";
}


/** Exract indegrees and save to double array */
void graph::__extractIndegrees(DeBruijnGraph &graph, InDegreeMap &in_degree, double indeg[])
{
	int i = 0;
	Vertex_iter vc, ve;
	for ( boost::tie(vc,ve) = vertices(graph); vc != ve; ++vc ) {
		indeg[i] = in_degree[*vc];
		i++;
	}
}

/** Exract outdegrees and save to double array */
void graph::__extractOutdegrees(DeBruijnGraph &graph, double outdeg[])
{
	int i = 0;
	Vertex_iter vc, ve;
	for ( boost::tie(vc,ve) = vertices(graph); vc != ve; ++vc ) {
		outdeg[i] = out_degree(*vc, graph);
		i++;
	}
}

/** Dispaly degree statistics */
void graph::degreeStats( DeBruijnGraph &graph ) 
{

    int n = num_vertices(graph);
    double indeg[n], outdeg[n];

    InDegreeMap in_degree = indegrees(graph);
    __extractIndegrees(graph, in_degree, indeg);
    __extractOutdegrees(graph, outdeg);

    math::sort(indeg, n);
    math::sort(outdeg,n);

    std::cout << "\nGraph Degree Statistics\n";
    std::cout << "\tin-degree\tout-degree\n";
    std::cout << "avg\t" << math::mean(indeg, n) << "\t" << math::mean(outdeg, n) << "\n";
    std::cout << "max\t"  << math::max(indeg, n) << "\t" << math::max(outdeg, n) << "\n";
    std::cout << "min\t"  << math::min(indeg, n) << "\t" << math::min(outdeg, n) << "\n";
    std::cout << "med\t"  << math::median(indeg, n, 1) << "\t" << math::median(outdeg, n, 1) << "\n";

    for ( int i = 10; i <= 90; i+= 10 ) {
        std::cout << i << " quantile\t" << math::quantile(indeg, n, i/100.0, 1) << "\t" << math::quantile(outdeg, n, i/100.0, 1) << "\n";
    }
    for ( int i = 91; i <= 99; i++ ) {
        std::cout << i << " quantile\t" << math::quantile(indeg, n, i/100.0, 1) << "\t" << math::quantile(outdeg, n, i/100.0, 1) << "\n";
    }
    for ( double i = 99.1; i <= 100; i+=0.1 ) {
        std::cout << i << " quantile\t" << math::quantile(indeg, n, i/100.0, 1) << "\t" << math::quantile(outdeg, n, i/100.0, 1) << "\n";
    }
    std::cout << "\n";
}




VertexToKmerMap graph::revertMap( KmerToVertexMap &kmer_map ) 
{
    VertexToKmerMap vertex_map;
    KmerToVertexMap::iterator it;
    for ( it = kmer_map.begin(); it != kmer_map.end(); ) {
        vertex_map[ it->second ] = it->first;
        it = kmer_map.erase(it);
    }
    return vertex_map;
}

