/*
 * assembly.cpp
 *
 *  Created on: Feb 15, 2012
 *      Author: Youngik Yang
 */

#include "assembly.h"


void assem::proceed( DeBruijnGraph &graph,
					 InvertedIndex &iindex, 
					 PathToAlnMap &path2aln_map,
					 VertexToKmerMap &vertex_map,
					 CoverageMap &kmer_coverage,
					 BitString *bstrs,
					 char *strands,
					 ReadId *pairs,
					 int &nreads, 
					 PathId *used_reads,
					 std::set<KmerType>  &debug_kmers,
					 Param &param )
{
    //------------------------------------------------------------
    // Trim weakly supported node and edges.
    // Then, remove corresponding entries from the inverted index.
    //------------------------------------------------------------
    trim(graph, kmer_coverage, vertex_map, iindex, param);
	
    //--------------------
    // Discover kmer paths
    //--------------------
    extractPaths(graph, bstrs, strands, pairs, nreads, vertex_map, iindex, kmer_coverage, used_reads, debug_kmers, path2aln_map, param );

    //--------------------------------
    // Merge or latch discovered paths
    //--------------------------------
    combinePaths(path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param );

    //---------------------
    // Recruit unused reads
    //---------------------
    recruitReads( path2aln_map, iindex, bstrs, strands, pairs, used_reads, nreads, param );

    //-----------------------------------------------------------------------
    // Latch paths with paried end information with no supporting latch reads
    //-----------------------------------------------------------------------
    connectPaths(path2aln_map, bstrs, strands, pairs, used_reads, nreads, iindex, param );
    //extendShortOverlapPaths(path2aln_map, bstrs, strands, pairs, used_reads, nreads, iindex, param );
      

    //-------------------------------
    // Correct misplaced paired-reads
    //-------------------------------
    mergePairedPaths(path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param );

    //tunePaths( path2aln_map, bstrs, strands, pairs, used_reads, nreads, iindex, param );
    //scaffoldPaths(path2aln_map, bstrs, strands, pairs, used_reads, param );
}


void assem::trim( DeBruijnGraph &graph, CoverageMap& kmer_coverage, VertexToKmerMap &vertex_map, InvertedIndex &iindex, Param &param )
{
    if ( ! param.trim_flag ) return;

    std::cout << "\nTrimming graph ...\n";
    double time1 = mytime();
    NodeSet trimmed = graph::trimGraph(graph, kmer_coverage, vertex_map, param.min_depth, param.verbose );
    std::cout << "Graph trimmed:" << mytime()-time1 << " sec\n";
	graph::graphSummary(graph, std::cout);    

    //---------------------------------------------
    // Now trim away singleton list from vertex map
    //---------------------------------------------
    trimMap(kmer_coverage, vertex_map, trimmed, param);
    trimReads(iindex, trimmed, vertex_map, param);

    cov::coverageSummary(iindex, std::cout);

//     UnDeBruijn ungraph;
//     std::tr1::unordered_map<Vertex, UnVertex> v2v_map;
//     convert(graph, ungraph, v2v_map);

//     computeCC(ungraph);

}

void assem::trimMap( CoverageMap &kmer_coverage,
              VertexToKmerMap &vertex_map,
              NodeSet &dropped,
			  Param &param )
{
    NodeSet::iterator it;
    for ( it = dropped.begin(); it != dropped.end(); ++it ) {
        if ( vertex_map.find(*it) == vertex_map.end() ) {
            if ( param.verbose ) std::cout << "\t** " << *it << ": vertex map does not exist\n"; 
            continue;
        }
        kmer_coverage.erase( kmer_coverage.find(vertex_map[*it]) );
    }
}

void assem::trimReads( InvertedIndex &iindex, 
                NodeSet &trimmed, 
                VertexToKmerMap &vertex_map,
				Param &param )
{
    double time1 = mytime();
    if ( param.verbose ) std::cout << "\nTrimming reads ...\n";
    NodeSet::iterator it;
    for ( it = trimmed.begin(); it != trimmed.end(); ++it ) {
        if ( vertex_map.find(*it) == vertex_map.end() ) {
            std::cerr << "\t** " << *it << " vertex map does not exist\n";
            exit(0);
        }
        iindex.erase( vertex_map[*it] );
    }
    if ( param.verbose ) std::cout << "... Read trimmed:" << mytime()-time1 << " sec\n";
}


std::pair<Vertex, double> assem::findMaxPriority( std::tr1::unordered_map<Vertex, double> &cov_map )
{
    Vertex v = NULL;
    double m = -1.0;

    std::tr1::unordered_map<Vertex, double>::iterator it = cov_map.begin();
    for ( ; it != cov_map.end(); ++it ) {
        if ( it->second > m ) {
            v = it->first;
            m = it->second;
        }
    }
    return std::pair<Vertex, double>(v,m);
}

void assem::dropBadPriorities(std::tr1::unordered_map<Vertex, double> &cov_map,
                              NodeSet &deleted_nodes)
{
    for (NodeSet::iterator it = deleted_nodes.begin(); it != deleted_nodes.end(); ++it )
        if ( cov_map.find(*it) != cov_map.end() ) 
            cov_map.erase(*it);

    std::tr1::unordered_map<Vertex, double>::iterator it = cov_map.begin();
    for ( ; it != cov_map.end(); ) {
        if ( it->second < 0 ) 
            it = cov_map.erase(it++);
        else ++it;
    }
}

// void assem::updatePriorities(PathType &max_path, 
//                       std::tr1::unordered_map<Vertex, double> &cov_map,
//                       DeBruijnGraph &graph,
//                       InvertedIndex &iindex, 
//                       VertexToKmerMap &vertex_map,
//                       AdjacencyMap &ancestors,
//                              NodeSet &deleted_nodes,                              
//                       Param &param )
// {
//     for ( PathType::iterator it = max_path.begin(); it != max_path.end(); ++it ) {
//         if ( deleted_nodes.find(*it) == deleted_nodes.end() ) {
//             updatePathSearchPriority( *it, cov_map, graph, iindex, vertex_map, ancestors, param );
//         } 
//         //else cov_map[*it] = -1;
//     }
//     dropBadPriorities(cov_map, deleted_nodes);
// }

void assem::extractPaths( DeBruijnGraph &graph,
                    BitString *bstrs,
                    char *strands,
                    ReadId *pairs,
                    int &nreads, 
                    VertexToKmerMap &vertex_map,
					InvertedIndex &iindex, 
                    CoverageMap &kmer_coverage,
                    PathId *used_reads,
                    std::set<KmerType>  &debug_kmers,
                    PathToAlnMap &path2aln_map,
                    Param &param )
{
    if ( !param.path_flag ) return;

//     /* debugging */
//     std::string consensus_file = param.out_dir + "/debug.fasta";
//     std::string statistic_file = param.out_dir + "/debug.cover";
//     std::string placement_file = param.out_dir + "/debug.place";
//     std::string alignment_file = param.out_dir + "/debug.align";
//     std::string profile_file   = param.out_dir + "/debug.profile";
//     std::fstream cons, stat, plac, alig, prof;

//     io::openFile(cons, consensus_file.c_str(), std::ios::out);
//     io::openFile(stat, statistic_file.c_str(), std::ios::out);
//     io::openFile(plac, placement_file.c_str(), std::ios::out);
//     io::openFile(alig, alignment_file.c_str(), std::ios::out);
//     io::openFile(prof, profile_file.c_str(), std::ios::out);

//     std::fstream patb;
//     std::string path_file   = param.out_dir + "/path." + param.dump_suffix;
//     io::openFile(patb, path_file.c_str(), std::ios::out | std::ios::binary );

//     cons.rdbuf()->pubsetbuf(0,0);
//     stat.rdbuf()->pubsetbuf(0,0);
//     plac.rdbuf()->pubsetbuf(0,0);
//     alig.rdbuf()->pubsetbuf(0,0);
//     prof.rdbuf()->pubsetbuf(0,0);
//     patb.rdbuf()->pubsetbuf(0,0);



    double time1 = mytime();

    setbuf(stdout, NULL); // no buffering
    //std::cout << "\nSTAGE 1:\nSearching paths ...\n";
    std::cout << "\nSearching paths ...\n";

    //--------------------------------
    // Set minimum depth of seed kemrs
    //--------------------------------
    setStartMinCoverage(iindex, param);

    AdjacencyMap ancestors = graph::predecessorMap(graph);
    std::multimap<double, Vertex> cov_map;
    //std::tr1::unordered_map<Vertex, double> cov_map;
    initPathSearchPriority( cov_map, graph, iindex, vertex_map, ancestors, param );

    NodeSet deleted_nodes;
    int nstart = 0;
    int nseeds = cov_map.size();
    double pratio = 99.0;
    std::multimap<double, Vertex>::reverse_iterator it;
    std::cout << "Seed kmer count:" << cov_map.size() << "\n";
    double lt0 = mytime();

    size_t c_good_rids, c_short_path, c_short_con, c_small_reads, c_msa_error, c_zero_depth, c_med_depth, c_zero_prof;
    c_good_rids = c_short_path = c_short_con = c_small_reads = c_msa_error = c_zero_depth = c_med_depth = c_zero_prof = 0;
    std::tr1::unordered_map<KmerId, bool> succ_kmers, fail_kmers;
    
    while (cov_map.size()) {
        //double st = mytime();

        //
        //double tt = mytime();
        //std::pair<Vertex, double> vp = findMaxPriority(cov_map);
        //if ( param.verbose ) std::cout << "Max priority search:" << mytime()-tt << " sec\n";
        //Vertex v = vp.first;
        //double p = vp.second;
        
        //if ( v == NULL ) break; 
        //if ( p == -1.0 ) break;
        it = cov_map.rbegin();
        Vertex v = it->second;

//         //---------------------------------
//         // Remove 1st element from the last
//         //---------------------------------
        cov_map.erase(--it.base());
        //cov_map[v] = -1; 
        if ( cov_map.size()/(double)nseeds*100 <= pratio ) {
            fprintf( stdout, "%d/%d: %.2f sec\n", int(nseeds-cov_map.size()), nseeds, mytime()-lt0);
            pratio -= 1;
        }
        
        //--------------------------------------
        // Check whether it is valid search node
        //--------------------------------------
        if ( ! validStartNode(v, iindex, ancestors, vertex_map, deleted_nodes, debug_kmers, param) ) continue;

        nstart++;

        //--------------------
        // Current pivot k-mer
        //--------------------
        KmerId kid = vertex_map[v];
        
        //--------------------
        // Depth of start kmer
        //--------------------
        size_t start_depth = iindex.getValue(kid)->size;

//         std::cout << "Seed:" << alpha::IntegerToAminoAcid(kid,param.kmer_size ) << "\tcoverage:" << start_depth << "\t";
//         std::cout << "#indeg:" << ancestors[v].size() << "\t" << "#outdeg:" << graph::successors(graph, v).size() << "\t";
        
        if ( param.verbose ) {
            std::cout << "\n*Start:" << alpha::IntegerToAminoAcid(kid,param.kmer_size ) << "\tcoverage:" << start_depth << "\t";
            std::cout << "#indeg:" << ancestors[v].size() << "\t" << "#outdeg:" << graph::successors(graph, v).size() << "\n";
        }


        /*===========================================*/
        /* Discover path with the current start node */
        /*===========================================*/
        PathType max_path;
        ReadIdArray path_rids;
        int lstop, rstop;
        //double et0=mytime();
        extractGreedyBestPath(max_path, path_rids, v, graph, iindex, vertex_map, ancestors, lstop, rstop, param);
        //double et1=mytime();
        //std::cout << "plength:" << max_path.size()+param.kmer_size-1 << "\tlstop:" << conditions[lstop] << "\trstop:" << conditions[rstop] << "\tstime:" << et1-et0 << "\t";

        if ( param.verbose ) 
            std::cout << "\nBest Path:" << biostr::getSequenceString(max_path, vertex_map, param.kmer_size) << "\n\n";

        /*=================*/
        /* Path Validation */
        /*=================*/
        ReadIdArray good_rids = trimUsedReads(path_rids, used_reads);
        //std::cout << "creads:" << good_rids.size() << "\t";
        if ( good_rids.size() < (size_t)param.min_share ) {
            //std::cout << "small_reads(" << good_rids.size() << ")\n";
            c_good_rids++; 
            fail_kmers.insert( std::pair<KmerId, bool>( kid, true ) );
            continue; 
        }

        if ( shortPath(max_path, vertex_map, param) ) {
            //std::cout << "short_paths\n";
            c_short_path++; 
            fail_kmers.insert( std::pair<KmerId, bool>( kid, true ) );
            continue; // Path is short
        }

        //---------------------------------------------
        // Pile of reads belongs to the discovered path
        //---------------------------------------------
        //double pt0 = mytime();
        SpaPath *spath = getPathReadPile(max_path, vertex_map, iindex, bstrs, good_rids, used_reads, param);
        //double pt1 = mytime();

        std::string consensus = spath->getConsensusString();
        if ( (int)consensus.size() <= param.kmer_size ) {
            //std::cout << "short_consensus(" << consensus.size() << ")\tatime:" << pt1-pt0 << "\n";
            c_short_con++;
            fail_kmers.insert( std::pair<KmerId, bool>( kid, true ) );
            if ( param.verbose ) std::cout << "\n** Too short path:" << consensus << "\n"; 
            spath->resetReads(used_reads, NOT_PATH);
            delete spath; continue;
        }
        
        //------------------------------------------
        // Bad path because of small number of reads
        //------------------------------------------
        if ( spath->getReadCount() < (unsigned)param.min_share ) {
            //std::cout << "small_aligned_reads(" << spath->getReadCount() << ")\tatime:" << pt1-pt0 << "\n";
            c_small_reads++;
            fail_kmers.insert( std::pair<KmerId, bool>( kid, true ) );
            if ( param.verbose ) std::cout << "\n** Weak supported path (#total reads:" << spath->getReadCount() 
                                     << "\t#reads in path:" << spath->getReadCount() << ")\n";
            
            spath->resetReads(used_reads, NOT_PATH);
            delete spath; continue;
        }
                           
        //--------------------------
        // Mapping reads to the path
        //--------------------------
        int error = 0;
        //double mt0 = mytime();
		MSA msa(*spath, (PathId)path2aln_map.size(), used_reads, bstrs, strands, pairs, PILE|FLAT|TRIM|CODON, error, param);
        //double mt1 = mytime();
        if ( error ) {
            //std::cout << "msa_error\tatime:" << pt1-pt0 << "\tmtime:" << mt1-mt0 << "\n";
            c_msa_error++;
            fail_kmers.insert( std::pair<KmerId, bool>( kid, true ) );
            if ( param.verbose ) std::cout << "MSA fail:\tcode:" << error << "\n";
            spath->resetReads(used_reads, NOT_PATH);
            delete spath; continue;
        }
        if ( param.verbose ) {
            std::cout << "\nAlignment\n"; msa.printAlignment(std::cout, *spath, 100);
            std::cout << "\nProfile\n";   msa.printProfile(std::cout);
        }

        //-------------------------------------------
        // Check existence of no read covered regions
        //------------------------------------------- 
        if ( msa.getZeroCoverageRegions().size() > 0 ) {
            //std::cout << "zero_depth_region\tatime:" << pt1-pt0 << "\tmtime:" << mt1-mt0 << "\n";
            c_zero_depth++;
            fail_kmers.insert( std::pair<KmerId, bool>( kid, true ) );
            if ( param.verbose ) 
                std::cout << "\n** Zero coverage region exists\n";
            spath->resetReads(used_reads, NOT_PATH);
            delete spath; continue;
        }

        //----------------------------------
        // Check coverage of each Amino Acid
        //----------------------------------
        ScoreSummary path_score = getPathDepth( msa, param );
        if ( path_score.med < param.med_depth ) {
            //std::cout << "poor_med_depth(" << path_score.med << ")\tatime:" << pt1-pt0 << "\tmtime:" << mt1-mt0 << "\n";
            c_med_depth++;
            fail_kmers.insert( std::pair<KmerId, bool>( kid, true ) );
            if ( param.verbose ) {
                std::cout << "\n** Weak supported (# reads:" << spath->getReadCount() << ")\n";
                std::cout << "Median support:" << path_score.med << "\n";
            }
            spath->resetReads(used_reads, NOT_PATH);
            delete spath; continue;
        }

        ProfileVector *vprof = spath->getProfileVector();
        if ( vprof->getSize() == 0 ) {
            //std::cout << "zero_length_profile\tatime:" << pt1-pt0 << "\tmtime:" << mt1-mt0 << "\n";
            c_zero_depth++;
            fail_kmers.insert( std::pair<KmerId, bool>( kid, true ) );
            if ( param.verbose ) std::cout << "Zero length profile\n";
            spath->resetReads(used_reads, NOT_PATH);
            delete spath; continue;
        }


        //---------------------------------
        // Check entropy of each Amino Acid
        //---------------------------------
        Profile profile = msa.getProfile();
        ScoreSummary entropy = computeEntropy( profile );

        //---------------------------------------------------------------
        // Save the current read pile mapping
        //-----------------------------------
        savePath( spath, path2aln_map );
        //PathId pid = savePath( spath, path2aln_map );

//         /* debugging */
//         std::cout << "PID:" << pid << "\t";
//         patb.write((char*)&(pid), sizeof(PathId));
//         spath->dump(patb);
//         std::cout << "... path dumped\n";

//         writeConsensus( cons, spath, pid );
//         std::cout << "... consensus written\n";
//         writePlacement( plac, spath, pid );
//         std::cout << "... placement written\n";
//         writeStatistic( stat, spath, pid );
//         std::cout << "... statistic written\n";
//         writeAlignment( alig, spath, pid, bstrs );
//         std::cout << "... alignment written\n";
//         writeProfile( prof, spath, pid );
//         std::cout << "... profile written\n";
        

        
        //-------------------------
        // Display the current path
        //-------------------------
        if ( param.verbose ) 
            printPathSummary(msa.getConsensus(), start_depth, kid, ancestors[v].size(), graph::successors(graph, v).size(), lstop, rstop, path_score, entropy, spath->getReadCount(), param);
        

        //-------------------------
        // drop reads found in path
        //-------------------------
        good_rids = ReadIdArray( spath->getReads(), spath->getReads() + spath->getReadCount() );
        adjustReadIds ( max_path, graph, good_rids, iindex, vertex_map, param.kmer_size );

        //----------------------------------------
        // remove weak kmers after read adjustment
        //----------------------------------------
        if ( param.trim_flag ) dropWeakNodes( deleted_nodes, max_path, graph, ancestors, iindex, vertex_map, param );

        //--------------------------------------------
        // Put vertex if it is still a good start kmer
        //--------------------------------------------
//         double tp = mytime();
//         updatePriorities(max_path, cov_map, graph, iindex, vertex_map, ancestors, deleted_nodes, param );
//         if ( param.verbose ) std::cout << "Update search priority:" << mytime()-tp << " sec\n";

        updatePathSearchPriority( v, cov_map, graph, iindex, vertex_map, ancestors, param );

        succ_kmers.insert( std::pair<KmerId, bool>( kid, true ) );

        //std::cout << "good_path\tatime:" << pt1-pt0 << "\tmtime:" << mt1-mt0 << "\tttime:" << mytime()-st << "\n";
    }

    size_t c_fail_kmers = 0;
    size_t c_both_kmers = 0;
    for ( std::tr1::unordered_map<KmerId, bool>::iterator fit = fail_kmers.begin(); fit != fail_kmers.end(); ++fit ) {
        if ( succ_kmers.find(fit->first) == succ_kmers.end() ) 
            c_fail_kmers++;
        else
            c_both_kmers++;
    }
    
//     std::cout << "\nFailure summary\n";
//     std::cout << "Success seed kmers:" << succ_kmers.size() << "\n"; 
//     std::cout << "Failure seed kmers:" << c_fail_kmers << "\n"; 
//     std::cout << "Good/bad seed kmers:" << c_both_kmers << "\n";
//     std::cout << "Small reads in path:" << c_good_rids << "\n";
//     std::cout << "Short path:" << c_short_path << "\n";
//     std::cout << "Short consensus:" << c_short_con << "\n";
//     std::cout << "Small reads in alignment:" << c_small_reads << "\n";
//     std::cout << "MSA error:" << c_msa_error << "\n";
//     std::cout << "Zero depth:" << c_zero_depth << "\n";
//     std::cout << "Low med depth:" << c_med_depth << "\n";
//     std::cout << "Zero length profile:" << c_zero_prof << "\n\n";

    if ( param.verbose ) displayPathSearchSummary(nstart, deleted_nodes.size(), path2aln_map.size(), nreads-getUnusedReadCount(used_reads, nreads));

    //------------------------
    // Trim away deleted nodes
    //------------------------
    for ( NodeSet::iterator it = deleted_nodes.begin(); it != deleted_nodes.end(); ++it ) {
        iindex.erase(vertex_map[*it]);
        remove_vertex(*it, graph);
    }

    std::cout << "Path count:" << path2aln_map.size() << "\n";
    std::cout << "Path search:" << mytime()-time1 << " sec\n";

    if ( param.verbose ) 
        graph::graphSummary(graph, std::cout);
}

ScoreSummary assem::getPathDepth(MSA &msa, Param &param)
{
    unsigned nc  = msa.getProfile().ncol;
    unsigned **m = msa.getProfile().matrix;
    double *csums = new double[nc];
    for ( unsigned i = 0; i < nc; i++ ) {
        csums[i] = 0;
    }
    for ( unsigned i = 0; i < nc; i++ ) {
        for ( int j = 0; j < naas; j++ ) {
            csums[i] += m[i][j];
        }
    }
    math::sort(csums,nc);
    double psum = math::sum( csums, nc );
    double pmin = math::min( csums, nc );
    double pmax = math::max( csums, nc );
    double pavg = math::mean( csums, nc );
    double pmed = math::median( csums, nc, true );
    if ( param.verbose ) std::cout << "\tRead support\tMax:" << pmax << "\tMin:" << pmin << "\tAvg:" << pavg << "\tMed:" << pmed << "\n";
    delete[] csums;
    
    return ScoreSummary(psum, pmin, pmax, pavg, pmed);
}

ScoreSummary assem::computeEntropy( Profile &profile )
{
    unsigned nc  = profile.ncol;
    unsigned **m = profile.matrix;

    double *entropies = new double[nc];
    for ( unsigned i = 0; i < nc; i++ ) {
        std::map<int, int> cmap;
        int sum = 0;
        for ( int j = 0; j < naas; j++ ) {
            if ( m[i][j] > 0 ) {
                cmap.insert( std::pair<int, int>(j, m[i][j]) );
                sum += m[i][j];
            }
        }
        double e = eval::entropy<int>(cmap, sum);
        entropies[i] = e;
    }
    math::sort(entropies, nc);
    double esum = math::sum( entropies, nc );
    double emin = math::min( entropies, nc );
    double emax = math::max( entropies, nc );
    double eavg = math::mean( entropies, nc );
    double emed = math::median( entropies, nc, true );
    delete[] entropies;
 
    return ScoreSummary(esum, emin, emax, eavg, emed);
}

ScoreSummary assem::computeDepth( Profile &profile )
{
    unsigned nc  = profile.ncol;
    unsigned **m = profile.matrix;

    double *depths = new double[nc];
    for ( unsigned i = 0; i < nc; i++ ) {
        std::map<int, int> cmap;
        int sum = 0;
        for ( int j = 0; j < naas; j++ ) {
            if ( m[i][j] > 0 ) {
                cmap.insert( std::pair<int, int>(j, m[i][j]) );
                sum += m[i][j];
            }
        }
        depths[i] = sum;
    }
    math::sort(depths, nc);
    double esum = math::sum( depths, nc );
    double emin = math::min( depths, nc );
    double emax = math::max( depths, nc );
    double eavg = math::mean( depths, nc );
    double emed = math::median( depths, nc, true );
    delete[] depths;
 
    return ScoreSummary(esum, emin, emax, eavg, emed);
}



void assem::setStartMinCoverage( InvertedIndex &iindex, Param &param )
{
    InvertedIndexMap imap = iindex.getInvertedIndex();
    int n = imap.size();
    double *coverage = new double[n];

    int i = 0;
    for (InvertedIndexMap::iterator it = imap.begin(); it != imap.end(); ++it )
        coverage[i++] = it->second->size;

    math::sort(coverage, n);        
    double value = math::quantile(coverage, n, (100-param.kmer_percentile)/100.0, 1);
    if ( param.min_seed > (int)value ) param.min_seed = (int) value;
    std::cout << "Min. coverage of seed kemrs:" << param.min_seed << "\n";
    delete[] coverage;
}


void assem::initPathSearchPriority( std::multimap<double, Vertex> &cov_map,
                                    //std::tr1::unordered_map<Vertex, double> &cov_map,
                                   DeBruijnGraph &graph,
                                   InvertedIndex &iindex,
                                   VertexToKmerMap &vertex_map,
                                   AdjacencyMap &ancestors,
                                   Param &param )
{
    Vertex_iter vc, ve;
    for ( boost::tie(vc,ve) = vertices(graph); vc != ve; ++vc ) {
        updatePathSearchPriority(*vc, cov_map, graph, iindex, vertex_map, ancestors, param);
    }
    //NodeSet empty;
    //dropBadPriorities(cov_map, empty);
    if ( param.verbose ) std::cout << "\n# Search nodes:" << cov_map.size() << "\n";
}

bool assem::validStartNode( Vertex &v, 
                     InvertedIndex &iindex, 
                     AdjacencyMap &ancestors, 
                     VertexToKmerMap &vertex_map, 
                     NodeSet &deleted_nodes, 
                     std::set<KmerType> &debug_kmers,
					 Param &param  )
{
        
    //------------------------------------------
    // Shouldn't occur (No vertex mapping exist)
    //------------------------------------------
    if ( vertex_map.find(v) == vertex_map.end() ) { 
        std::cerr << "[Error] Vertex map for :" << v << " does not exit\n"; 
        exit(-1); 
    }
    
    //--------------------
    // Current pivot k-mer
    //--------------------
    KmerId kid = vertex_map[v];
    
    //-------------------
    // Invalid start kmer
    //-------------------
    if ( alpha::getFirstAminoAcid<KmerId>(kid, param.kmer_size) == '*') {
        if ( param.verbose )  std::cout << "** Invalid Start Kmer:" << alpha::IntegerToAminoAcid(kid, param.kmer_size ) << "\tskip\n"; 
        return false;
    }
    
    //---------------------------
    // Skip if deleted previously
    //---------------------------
    if ( deleted_nodes.find(v) != deleted_nodes.end() )  {
        if ( param.verbose )  std::cout << "** Deleted:" << alpha::IntegerToAminoAcid(kid,param.kmer_size ) << "\tskip\n"; 
        return false;
    }
    
    //---------------------------------------------
    // Skip if not a debuggin kmer in debuggin mode
    //---------------------------------------------
    if ( debug_kmers.size() > 0  && debug_kmers.find(alpha::IntegerToAminoAcid(kid,param.kmer_size)) == debug_kmers.end() ) {
        return false;
    }
    
    //-----------------------------------        
    // Shouldn't occur (No index mapping)
    //-----------------------------------        
    if ( ! iindex.has(kid) ) { 
        if ( param.verbose ) std::cerr << "** Inverted index does not exist: " << alpha::IntegerToAminoAcid(kid,param.kmer_size ) << "\n"; 
        exit(1);
    }
    
    int size = iindex.getValue(vertex_map[v])->size;
    //-----------------------------
    // Skip if either low supported
    //-----------------------------
    if ( size < param.min_seed ) {
        if ( param.verbose ) std::cout << "** Low read coverage:" << alpha::IntegerToAminoAcid(kid,param.kmer_size ) << "\tcoverage:" << size << "\tskip\n";
        return false;
    }
    
    //-------------------------------------
    // Shouldn't occur (No predecessor map)
    //-------------------------------------
    if ( ancestors.find(v) == ancestors.end() ) {
        if ( param.verbose ) std::cerr << "** Predecessor map not defined: " << alpha::IntegerToAminoAcid(kid,param.kmer_size ) << "\n"; 
        exit(1);
        //continue;
    }
    return true;
}


void assem::extractGreedyBestPath( PathType &best_path, 
							ReadIdArray &good_rids, 
							Vertex curr, 
							DeBruijnGraph &graph, 
							InvertedIndex &iindex, 
							VertexToKmerMap &vertex_map, 
							AdjacencyMap &ancestors,
							int &lstop,
							int &rstop,
							Param &param )
{
    double t0 = mytime();

    NodeArray preds = ancestors[curr];
    NodeArray succs = graph::successors(graph, curr);
    
    best_path.push_back(curr);
   
    std::queue<NodeList> rid_queue;
    ReadIdArray poor_rids, ncyc_rids;

    KmerId kid = vertex_map[curr];    
    good_rids = ReadIdArray( iindex.getValue(kid)->rid, iindex.getValue(kid)->rid + iindex.getValue(kid)->size );
    std::sort(good_rids.begin(), good_rids.end());

    PathType cyclic_path;
    Vertex cycle_start = NULL;
    Vertex cycle_trigger = NULL;

    int reason;
    int direction = LEFT;
    if ( param.verbose )  {
        direction == LEFT ? std::cout << "(-) " : std::cout << "(+) " ;
        std::cout << alpha::IntegerToAminoAcid(kid, param.kmer_size) << "\n";
    }

    extendGreedyBestNode( best_path, cyclic_path, cycle_start, cycle_trigger, rid_queue, good_rids, poor_rids, ncyc_rids, graph, iindex, vertex_map, ancestors, direction, reason, param );
    direction == LEFT ? lstop = reason : rstop = reason;
    

    //---------------------------------------------------
    // reset and prepare to extend the opposite direction
    //---------------------------------------------------
    cyclic_path.clear();
    clear<NodeList>(rid_queue);
    direction == LEFT ? kid  = vertex_map[best_path.back()] : kid  = vertex_map[best_path.front()] ;
    ReadIdArray curr_rids = ReadIdArray( iindex.getValue(kid)->rid, iindex.getValue(kid)->rid + iindex.getValue(kid)->size );
    std::sort(curr_rids.begin(), curr_rids.end());
    poor_rids = Set::Difference<ReadId>( curr_rids, good_rids, true );
    if ( param.verbose ) std::cout << "\nGood size:" << good_rids.size() << "\tPoor size:" << poor_rids.size() << "\n";

    ncyc_rids = good_rids;

    cycle_start = cycle_trigger = NULL;
    //---------------------------------
    // extension in the other direction
    //---------------------------------
    direction = RIGHT;
    if ( param.verbose ) {
        direction == LEFT ? std::cout << "(-) " : std::cout << "(+) " ;
        std::cout << alpha::IntegerToAminoAcid(kid, param.kmer_size) << "\n";
    }
    extendGreedyBestNode( best_path, cyclic_path, cycle_start, cycle_trigger, rid_queue, good_rids, poor_rids, ncyc_rids, graph, iindex, vertex_map, ancestors, direction, reason, param );
    direction == LEFT ? lstop = reason : rstop = reason;

    if ( param.verbose ) std::cout << "\n** Best path search (length:" << best_path.size()+param.kmer_size-1 << "\t#reads:" << good_rids.size() << "\ttime:" << mytime() - t0 << " sec)\n";
}

bool assem::shortPath( PathType &max_path,
                VertexToKmerMap &vertex_map, 
				Param &param )
{
    std::string path_str = biostr::getSequenceString(max_path, vertex_map, param.kmer_size);
    if ( (int)path_str.length() < param.min_length ) {
        if (param.verbose) std::cout << "\nShort Path\n";
        return true;
    }
    return false;
}

void assem::extendGreedyBestNode( PathType &best_path,
                   PathType &cyclic_path,
                   Vertex &cycle_start,
                   Vertex &cycle_trigger,
                   std::queue<NodeList> &rid_queue,
                   ReadIdArray &good_rids, 
                   ReadIdArray &poor_rids,
                   ReadIdArray &ncyc_rids,
                   DeBruijnGraph &graph, 
                   InvertedIndex &iindex,
                   VertexToKmerMap &vertex_map, 
                   AdjacencyMap &ancestors,
                   int direction,
                   int &reason,
                   Param &param )
{
    double t0 = mytime();

    double tic, toc;

    if ( param.verbose ) {
        std::cout << "--------------------------------------------------------------------------------\n";
        std::cout << "Best:" << biostr::getSequenceString(best_path, vertex_map, param.kmer_size) << "\n";
    }
    Vertex curr;
    direction == LEFT ? curr = best_path.front() : curr = best_path.back();
    KmerId cur_kid = vertex_map[curr];

    if ( direction == RIGHT && alpha::getLastAminoAcid<KmerId>( cur_kid, param.kmer_size ) == '*' ) {
        if ( param.verbose ) std::cout << "Stop extendsion\n";
        reason = PATHEND;
        return;
    }
    
    /* Path reads so far
       If path extension fails, use this values.
     */
    ReadIdArray good_reads = good_rids;

    ReadIdArray curr_reads = ReadIdArray( iindex.getValue(cur_kid)->rid,  iindex.getValue(cur_kid)->rid + iindex.getValue(cur_kid)->size );
    std::sort(curr_reads.begin(), curr_reads.end());

    ReadIdArray comm_reads = Set::Intersect<ReadId>( good_rids, curr_reads, true );

    if ( param.verbose ) 
        std::cout << "Curr:" 
                  << alpha::IntegerToAminoAcid(vertex_map[curr], param.kmer_size) << " ("
                  << ancestors[curr].size() << "|"
                  << graph::successors(graph, curr).size() << "|" 
                  << iindex.getValue(vertex_map[curr])->size << "|" 
                  << comm_reads.size() << ")\n";
    
    
    std::pair<Vertex, int> nd = getTraceNode( best_path, direction, param );
    Vertex back = nd.first;
    int dist = nd.second;
    
    ReadIdArray back_comm;
    if ( back != NULL ) {
        KmerId back_kid = vertex_map[back];
        ReadIdArray back_reads = ReadIdArray( iindex.getValue(back_kid)->rid, iindex.getValue(back_kid)->rid + iindex.getValue(back_kid)->size );
        std::sort(back_reads.begin(), back_reads.end());
        back_comm = Set::Intersect<ReadId>( good_rids, back_reads, true );
        
        ReadIdArray comm  = Set::Intersect<ReadId>(back_comm, curr_reads, true);
        if ( comm.size() == 0 ) {
            if ( param.verbose ) std::cout << "** No common reads in distance " << dist << "\n"; 
            reason = NONCOMMON;
            return;
        }
        if ( param.verbose ) {
            std::cout << "\ndist\tkmer1\tkmer2\tdepth1\tdepth2\tcommon\n";
            std::cout << dist << "\t" 
                      << alpha::IntegerToAminoAcid(cur_kid, param.kmer_size) << "\t"
                      << alpha::IntegerToAminoAcid(back_kid, param.kmer_size) << "\t"
                      << iindex.getValue(cur_kid)->size << "\t"
                      << iindex.getValue(back_kid)->size << "\t"
                      << comm.size() << "\n";
        }
    }
    


    size_t init_depth = comm_reads.size();
    size_t last_depth;

    PathQueue queue; queue.push(best_path);
    int depth = best_path.size()+1;

    direction == RIGHT ? 
        graph::BFS(queue, graph, depth) :
        graph::reverseBFS(queue, graph, ancestors, depth);

    Vertex max_node = NULL;
    ReadIdArray max_comm;
    size_t sum_reads;
    std::tr1::unordered_map<Vertex, ReadIdArray> back_rids;
    
    NodeSet vset = findGoodNeighbors(graph, curr, queue, vertex_map, direction, param); // determine path ends
    if ( vset.size() == 0 ) {
        reason = PATHEND;
        if ( param.verbose ) std::cout << "** Path ends\n"; return; 
    }

    findGreedyBestNeighbor(vset, graph, ancestors, curr, back_comm, curr_reads, good_rids, max_node, max_comm, sum_reads, back_rids, direction, vertex_map, iindex, param );
    
    if ( max_node == NULL ) { 
        if ( param.verbose ) std::cout << "** Extension fail\n"; 
        reason = EXTENDFAIL;
        return; 
    }

    if ( alpha::getFirstAminoAcid<KmerId>(vertex_map[max_node], param.kmer_size) == '*') {
        if ( param.verbose ) std::cout << "Bad extension\n"; 
        reason = EXTENDFAIL;
        return;
    }
    if ( max_node == curr ) { 
        if ( param.verbose ) std::cout << "** Self loop\n"; 
        reason = SELFLOOP;
        return; 
    }

    if ( param.verbose ) {
        std::cout << "Max k-mer:" << alpha::IntegerToAminoAcid(vertex_map[max_node], param.kmer_size)  << "\n";
        std::cout << "# reads: (curr:" << comm_reads.size() << "\t"
                  << "sum:" << sum_reads << "\tmax:" << max_comm.size() << ")\n"; 
    }
//     if ( sum_reads * 2 < comm_reads.size() ) {
//         if ( param.verbose ) std::cout << "\n** Weak read support - sum:" << sum_reads << "\n";  //return; 
//     }

    if ( (int)max_comm.size() < param.min_share ) {
        if ( param.verbose ) std::cout << "\n** Weak read support - max:" << max_comm.size() << "\n";
        reason = WEAKMAX;
        return; 
    }



    if ( graph::formCycle(best_path, max_node) ) {
        if ( graph::formCycle(cyclic_path, max_node) ) {
            if ( param.verbose ) std::cout << "\n** Cyclic path:" << biostr::getSequenceString(cyclic_path, vertex_map, param.kmer_size) << "\n";
            NodeArray all_nodes = NodeArray(best_path.begin(), best_path.end());
            NodeArray sub_nodes;
            direction == LEFT ? sub_nodes = NodeArray(all_nodes.begin()+cyclic_path.size(),  all_nodes.end()) :
                sub_nodes = NodeArray(all_nodes.begin(), all_nodes.begin() + all_nodes.size()-cyclic_path.size());
            best_path = PathType(sub_nodes.begin(), sub_nodes.end());

            good_rids = ncyc_rids;
            reason = CYCLICPATH;
            return;
        }
        if ( cyclic_path.size() == 0 ) { 
            if ( param.verbose ) std::cout << "\n** Cycle start\n";
            cycle_trigger = curr;
            cycle_start = max_node;
            cyclic_path.push_back(max_node);
        } else if ( cyclic_path.size() > 0 ) {
            if ( param.verbose ) std::cout << "\n** Cycle extend\n";
            if ( direction == LEFT ) 
                cyclic_path.push_front(max_node);
            else
                cyclic_path.push_back(max_node);
        }
    } else {
        if ( cyclic_path.size() ) {
            //            std::cout << "Cycle end\n";
//             std::cout << "Curr:" << alpha::IntegerToAminoAcid(vertex_map[curr], k) << "\n";
//             std::cout << "Max Node:" << alpha::IntegerToAminoAcid(vertex_map[max_node], k) << "\n";
            //std::cout << "Trigger:" << cycle_trigger << "\t" << alpha::IntegerToAminoAcid(vertex_map[cycle_trigger],k) << "\n";
            //std::cout << "Start:" << cycle_start << "\t" << alpha::IntegerToAminoAcid(vertex_map[cycle_start],k) << "\n";
        }
        if ( cycle_trigger == curr ) {
            if ( param.verbose ) {
                std::cout << "\nRedundant path\n";
                std::cout << "Old:" << biostr::getSequenceString(best_path, vertex_map, param.kmer_size) << "\n";
                std::cout << "Cyc:" << biostr::getSequenceString(cyclic_path, vertex_map, param.kmer_size) << "\n";
            }
            NodeArray nodes = NodeArray(best_path.begin(), best_path.end());
            if ( direction == LEFT )
                nodes.erase( nodes.begin(), nodes.begin()+cyclic_path.size() );
            else
                nodes = NodeArray(nodes.begin(), nodes.end()-cyclic_path.size() );

            best_path = PathType(nodes.begin(), nodes.end());
            if ( param.verbose ) std::cout << "New:" << biostr::getSequenceString(best_path, vertex_map, param.kmer_size) << "\n";
        }
        ncyc_rids = good_rids;
        cyclic_path.clear();
        cycle_trigger = cycle_start = NULL;
    }


    bool stop = false;

    KmerId max_kid = vertex_map[max_node];

    if ( param.verbose ) 
        direction == LEFT ?
            std::cout << "\nDiscarding reads (other -> curr) ...\n" :
            std::cout << "\nDiscarding reads (other -> max)...\n";
    
    NodeList nodes;
    ReadIdList del_list;

    tic = mytime();
    NodeArray preds;
    direction == LEFT ? preds = ancestors[curr] : preds = ancestors[max_node];
    for ( NodeArray::iterator it = preds.begin(); it != preds.end(); ++it ) {
        if ( direction == LEFT )  {
            if ( *it == max_node ) continue;
        } else {
            if ( *it == curr ) continue;
            if ( *it == max_node ) continue; //self
        }
        KmerId kid = vertex_map[*it];
        if ( param.verbose ) {
            std::cout << "\t" 
                      << alpha::IntegerToAminoAcid(kid, param.kmer_size) << " -> ";
            direction == LEFT ? std::cout << alpha::IntegerToAminoAcid(vertex_map[curr], param.kmer_size) :
                std::cout << alpha::IntegerToAminoAcid(max_kid, param.kmer_size);
            std::cout << "\t" << iindex.getValue(kid)->size << " (";
            if ( direction == LEFT ) { 
                if ( param.verbose ) std::cout << back_rids[*it].size() << ")\n";
            }
            else {
                ReadIdArray rids = ReadIdArray(iindex.getValue(kid)->rid, iindex.getValue(kid)->rid + iindex.getValue(kid)->size);
                std::sort(rids.begin(), rids.end());
                ReadIdArray comm2 = Set::Intersect<ReadId>(good_rids, rids, true);
                if ( param.verbose ) std::cout << comm2.size() << ")\n";
            }
        }

        if ( !graph::hasNode( best_path, *it )  ) 
             nodes.push_back(*it);
        else {
            if ( param.verbose ) std::cout << "\t! Node in path\t skip\n";
            continue; // do not delete reads discovered before
        }
        append( del_list, iindex.getValue(kid)->rid, iindex.getValue(kid)->size );
    }
        
    toc = mytime();
    //std::cout << "$ Cleaning in-edges:" << toc-tic << " sec\n";



    if ( param.verbose ) 
        direction == LEFT ?
            std::cout << "\nDiscarding reads (max -> other) ...\n" :
            std::cout << "\nDiscarding reads (curr -> other) ...\n";
    
    tic = mytime();
    NodeArray succs;
    direction == LEFT ? succs = graph::successors(graph, max_node) : succs = graph::successors(graph, curr);
    for ( NodeArray::iterator it = succs.begin(); it != succs.end(); ++it ) {
        if ( direction == LEFT ) {
            if ( *it == curr ) continue;
            if ( *it == max_node ) continue; // self
        } else {
            if ( *it == max_node ) continue;
        }
        KmerId kid = vertex_map[*it];
        if ( param.verbose ) {
            std::cout << "\t";
            direction == LEFT ? std::cout << alpha::IntegerToAminoAcid(max_kid, param.kmer_size) :
                std::cout << alpha::IntegerToAminoAcid(vertex_map[curr],param.kmer_size) ;
            std::cout << " -> " << alpha::IntegerToAminoAcid(kid, param.kmer_size)
                      << "\t" << iindex.getValue(kid)->size << " (";
            if ( direction == LEFT ) {
                ReadIdArray rids = ReadIdArray(iindex.getValue(kid)->rid, iindex.getValue(kid)->rid + iindex.getValue(kid)->size);
                std::sort(rids.begin(), rids.end());
                ReadIdArray comm2 = Set::Intersect<ReadId>(good_rids, rids, true);
                std::cout << comm2.size() << ")\n"; 
            } else 
                std::cout << back_rids[*it].size() << ")\n";
        }
        if ( !graph::hasNode( best_path, *it )  ) 
             nodes.push_back(*it);
        else {
            if ( param.verbose ) std::cout << "\t! Node in path\tskip\n";
            continue; // do not delete reads discovered before
        }
        append( del_list, iindex.getValue(kid)->rid, iindex.getValue(kid)->size );
    }

        
    toc = mytime();
    //std::cout << "$ Cleaning out-eges:" << toc-tic << " sec\n";

    tic = mytime();
    if ( del_list.size() > 0 ) {
        ReadIdArray del_rids(del_list.begin(), del_list.end());
        std::sort(del_rids.begin(), del_rids.end());
        invalidate( good_rids, poor_rids, del_rids, true );
    }
    toc = mytime();
    //std::cout << "$ Updating supports:" << toc-tic << " sec\n";

    tic = mytime();
    //std::cout << "@ Updating queue:\t#old:" << poor_rids.size() << "\t";
    //std::cout << "\n# bad reads:" << poor_rids.size() << " -> ";
    updateQueue( rid_queue, nodes, poor_rids, iindex, vertex_map, 35 );
    toc = mytime();
    //std::cout << poor_rids.size() << "\n";
    //std::cout << "#new:" << poor_rids.size() << "\t" << toc-tic << " sec\n";


    tic = mytime();

    curr_reads = Set::Intersect<ReadId>( good_rids, curr_reads, true );
    ReadIdArray max_rids = ReadIdArray( iindex.getValue(max_kid)->rid, iindex.getValue(max_kid)->rid + iindex.getValue(max_kid)->size );
    std::sort(max_rids.begin(), max_rids.end());
    ReadIdArray new_rids = getNewSupports( good_rids, poor_rids, max_rids );//iindex.getValue(max_kid)->rid, iindex.getValue(max_kid)->size );
        
    good_rids = Set::Union<ReadId>( good_rids, new_rids, true );
        
    toc = mytime();
    //std::cout << "$ New supports:" << toc-tic << " sec\n";
        
    tic = mytime();
    ReadIdArray next_rids = Set::Intersect<ReadId>( good_rids, max_rids, true );
    if ( param.verbose ) 
        std::cout << "\n# seqlen:" << (best_path.size()+param.kmer_size-1)
                  << "\t# curr:" << curr_reads.size() 
                  << "\t# new:" << new_rids.size() 
                  << "\t# next:" << next_rids.size()
                  << "\t# good:" << good_rids.size() 
                  << "\t# bad:" << poor_rids.size() 
                  << "\t# queue:" << rid_queue.size() << "\n";
    toc = mytime();
    //std::cout << "$ Display stat:" << toc-tic << " sec\n";

    last_depth = curr_reads.size();

    if ( (int)curr_reads.size() < param.min_share ) {
        if ( param.verbose ) std::cout << "\n** Poor read supports (current)\n"; 
        reason = WEAKCURRENT;
        stop = true;
    }

    if ( new_rids.size() > curr_reads.size() ) {
        if ( param.verbose ) std::cout << "\n** Sudden increase of new reads\n"; //stop = true;
    }

    if ( init_depth > 2*last_depth ) {
        if ( param.verbose )  std::cout << "\n** Sudden decrease of supporting reads\n"; //stop = true;
    }

    //    tracePath( best_path, good_rids, graph, ancestors, iindex, vertex_map, param.kmer_size, direction, stop, 20);



    if ( param.verbose ) std::cout << "\n$ Path extension:" << mytime()-t0 << " sec\n";
    if ( stop ) { good_rids = good_reads; return; }

    if ( direction == LEFT ) {
        best_path.push_front(max_node);
    } else {
        best_path.push_back(max_node);
    }
    //extendGreedyBestNode( best_path, cyclic_path, cycle_start, cycle_trigger, rid_queue, good_rids, poor_rids, ncyc_rids, graph, iindex, vertex_map, ancestors,  k, min_support, max_trace, direction, reason );
    extendGreedyBestNode( best_path, cyclic_path, cycle_start, cycle_trigger, rid_queue, good_rids, poor_rids, ncyc_rids, graph, iindex, vertex_map, ancestors, direction, reason, param );
}



void assem::append( ReadIdList &reads, 
             ReadId rids[],
             size_t count ) 
{
    for ( size_t i = 0; i < count; i++ ) 
        reads.push_back(rids[i]);
}

void assem::invalidate( ReadIdArray &good_rids,
                 ReadIdArray &poor_rids,
                 ReadIdArray &del_rids,
                 bool sorted )
{
    good_rids = Set::Difference<ReadId>( good_rids, del_rids, true );
    poor_rids = Set::Union<ReadId>( poor_rids, del_rids, sorted );
}


void assem::updateQueue( std::queue<NodeList> &nq,
                  NodeList &nodes,                                      
                  ReadIdArray &poor_rids, 
                  InvertedIndex &iindex, 
                  VertexToKmerMap &vertex_map, 
                  size_t q_size )
{
    if ( nq.size() == q_size ) {
        ReadIdList list;
        for ( NodeList::iterator nt = nq.front().begin(); nt != nq.front().end(); ++nt ) {
            KmerId kid = vertex_map[*nt];
            for ( size_t i = 0; i < iindex.getValue(kid)->size; i++ ) 
                list.push_back( iindex.getValue(kid)->rid[i] );
        }
        
        ReadIdArray old_rids = ReadIdArray( list.begin(), list.end() );        

        std::sort(old_rids.begin(), old_rids.end());
        poor_rids = Set::Difference<ReadId>( poor_rids, old_rids, true);

        nq.pop();
    } 
    nq.push(nodes);
}


std::pair<Vertex, int> assem::getTraceNode( PathType &best_path,
                                     int direction,
									 Param &param )
{
    if ( best_path.size() <= 1 ) return std::pair<Vertex, int>(best_path.front(), 0);
    
    PathType apath = best_path;
    // drop curr node
    direction == LEFT ? apath.pop_front() : apath.pop_back();
    NodeArray nodes;
    direction == LEFT ?
        nodes = NodeArray(apath.begin(), apath.end()) :
        nodes = NodeArray(apath.rbegin(), apath.rend());

    int i = 0; 
    Vertex vv; 
    if ( (int)nodes.size() < param.back_trace ) {
        vv = nodes.back();
        i = nodes.size();
    }
    else {
        for ( NodeArray::iterator it = nodes.begin(); it != nodes.end(); ++it ) {
            ++i;
            if ( i == param.back_trace ) {
                vv = *it; break;
            }
        }
    }
    return std::pair<Vertex, int>(vv, i);
}

NodeSet assem::findGoodNeighbors( DeBruijnGraph &graph,
                           Vertex curr_node, 
                           PathQueue queue,  // pass by copy
                           VertexToKmerMap &vertex_map,
                           int direction,
                           Param &param )
{
    Paths paths = QueueToList<PathType>(queue);

    NodeSet vset;
    for ( Paths::iterator it = paths.begin(); it != paths.end(); ++it )  {
        PathType path = *it;
        Vertex v;
        direction == LEFT ? v = (*it).front() : v = (*it).back();
        if ( v == NULL ) continue;
        
        if ( vertex_map.find(v) == vertex_map.end() ) {
            std::cerr << "\n[ERROR] " << v << "\tVertex map does not exist\n";
            exit(1);
        }

        std::pair<Edge, bool> ep;
        if ( direction == LEFT ) 
            ep = edge(v, curr_node, graph) ;
        else
            ep = edge(curr_node, v, graph) ;

        if ( ep.second == false ) {
            std::cerr << "\n[ERROR] Edge to " << alpha::IntegerToAminoAcid(vertex_map[v], param.kmer_size) << " not exist\n";
            exit(1);
        }

        EdgeWeight weight;
        direction == LEFT ? weight = graph[edge(v, curr_node, graph).first].weight : weight = graph[edge(curr_node, v, graph).first].weight;
        if ( weight < param.min_share ) {
            if ( param.verbose ) std::cout << "\t" << alpha::IntegerToAminoAcid(vertex_map[v], param.kmer_size) << " : weak edge support (" << weight << ")\n";
            continue;
        }
        vset.insert(v);
    }
    return vset;
}


/**
 * Find the best neighboring extension node.
 */
void assem::findGreedyBestNeighbor( NodeSet &vset,
                         DeBruijnGraph &graph,
                         AdjacencyMap &ancestors,
                         Vertex curr_node,
                         ReadIdArray &back_comm,
                         ReadIdArray &curr_rids,
                         ReadIdArray &good_rids,
                         Vertex &max_node,
                         ReadIdArray &max_comm,
                         size_t &sum_reads,
                         std::tr1::unordered_map<Vertex, ReadIdArray> &back_rids, 
                         int direction,
                         VertexToKmerMap &vertex_map,
                         InvertedIndex &iindex,
                         Param &param
                         )
{
    if ( vset.size() == 0 ) return;

    if ( param.verbose ) 
        direction == LEFT ? std::cout << "\nPrev k-mers:\n" : std::cout << "\nNext k-mers:\n";

    sum_reads = 0;
    ReadIdArray comm_cp, comm_cb, comm_cg;
    std::multimap<int, Vertex> size2node_map;
    std::map<Vertex, ReadIdArray> cur2next_map, cur2back_map;
    
    Vertex      back_node = NULL;;
    ReadIdArray comm_rids;

 
    //------------------------------
    // Neighbor node search priority
    //------------------------------
    std::multimap<size_t, Vertex> depths;
    for ( NodeSet::iterator it = vset.begin(); it != vset.end(); ++it )  {
        depths.insert( std::pair<size_t, Vertex>( iindex.getValue(vertex_map[*it])->size, *it ) );
    }
    
    std::multimap<size_t, Vertex>::reverse_iterator it;
    for ( it = depths.rbegin(); it != depths.rend(); ++it ) {
        Vertex v = it->second;
        KmerId kid = vertex_map[v];
        if ( back_node == NULL ) back_node = v;

        if ( iindex.getValue(kid)->size < comm_rids.size() ) {
            if ( param.verbose ) {
                std::cout << "\t" 
                          << alpha::IntegerToAminoAcid(vertex_map[v], param.kmer_size) 
                          << "\t"
                          << "i:" << ancestors[v].size()
                          << " o:"  << graph::successors(graph,v).size() 
                          << " c:"  << iindex.getValue(vertex_map[v])->size
                          << " - skip\n";
            }
            continue;
        }

        ReadIdArray crids = ReadIdArray( iindex.getValue(kid)->rid, iindex.getValue(kid)->rid + iindex.getValue(kid)->size);
        std::sort(crids.begin(), crids.end());

        comm_cp = Set::Intersect<ReadId>( curr_rids, crids, true );
        if ( comm_cp.size() == 0 || comm_cp.size() < comm_rids.size() ) {
            if ( param.verbose ) {
                std::cout << "\t" 
                          << alpha::IntegerToAminoAcid(vertex_map[v], param.kmer_size) 
                          << "\t"
                          << "i:" << ancestors[v].size()
                          << " o:"  << graph::successors(graph,v).size() 
                          << " c:"  << iindex.getValue(vertex_map[v])->size
                          << " n:" << comm_cp.size() 
                          << " - weak sharing\n";
            }
            continue;
        }

        comm_cb = Set::Intersect<ReadId>( crids, back_comm, true );        
        if ( comm_cb.size() > comm_rids.size() ) {
            back_node = v; comm_rids = comm_cb;
        }

        comm_cg = Set::Intersect<ReadId>( good_rids, crids, true );

        EdgeWeight weight;
        direction == LEFT ? weight = graph[edge(v, curr_node, graph).first].weight : weight = graph[edge(curr_node, v, graph).first].weight;
        if ( param.verbose ) std::cout << "\t" 
                                 << alpha::IntegerToAminoAcid(vertex_map[v], param.kmer_size) 
                                 << "\t"
                                 << "i:" << ancestors[v].size()
                                 << " o:"  << graph::successors(graph,v).size() 
                                 << " c:"  << iindex.getValue(vertex_map[v])->size
                                 << " n:" << comm_cp.size() 
                                 << " g:" << comm_cg.size() 
                                 << " b:" << comm_cb.size()
                                 << "\n";
        
        back_rids[v] = comm_cb;
        sum_reads += comm_cp.size();        

        if ( comm_cb.size() > 0 ) {
            size2node_map.insert(std::pair<int, Vertex>(comm_cb.size(), v));
            cur2next_map.insert(std::pair<Vertex, ReadIdArray>(v, comm_cg));
            cur2back_map.insert(std::pair<Vertex, ReadIdArray>(v, comm_cb));
        }
    }

    if ( size2node_map.size() == 0 ) return;

    std::multimap<int, Vertex>::reverse_iterator rt = size2node_map.rbegin();
    Vertex v = rt->second;
    max_node = v;
    max_comm = cur2next_map[v];
}

ReadIdArray assem::getNewSupports ( ReadIdArray &good_rids,
                             ReadIdArray &poor_rids, 
                             ReadIdArray &curr_rids )
{
    ReadIdArray curr_good = Set::Difference<ReadId>( curr_rids, poor_rids, true );
    return Set::Difference<ReadId>( curr_good, good_rids, true );
}



ReadIdArray assem::trimUsedReads( ReadIdArray &path_rids,
                           PathId *used_reads )
{
    ReadIdList rlist;
    for ( ReadIdArray::iterator it = path_rids.begin(); it != path_rids.end(); ++it )
        if ( used_reads[*it] == NOT_PATH ) rlist.push_back(*it);

    return ReadIdArray( rlist.begin(), rlist.end() );
}

SpaPath* assem::getPathReadPile( PathType &max_path,
                         VertexToKmerMap &vertex_map,
                         InvertedIndex &iindex,
                         BitString *bstrs,
                         ReadIdArray &good_rids, 
                         PathId *used_reads,
						 Param &param )
{
    //------------------------
    // align reads in the path
    //------------------------
    NodeArray nodes = NodeArray(max_path.begin(), max_path.end());
    size_t nreads = good_rids.size();
    size_t nnodes = max_path.size();
    KmerArray kmers;
    for ( size_t i = 0; i < nnodes; i++ )
        kmers.push_back(vertex_map[nodes[i]]);
    
    double t0 = mytime();
    std::string seq = biostr::getSequenceString( max_path, vertex_map, param.kmer_size );
    SpaPath *spath = new SpaPath( seq.c_str(), &kmers[0], &good_rids[0], nnodes, nreads );
    spath->align( iindex, bstrs, param );
    spath->validate(bstrs, param);
    double t1 = mytime();
    if ( param.verbose ) std::cout << "\n$ Read mapping:" << t1 - t0 << " secs\n";

    return spath;
}

PathId assem::savePath( SpaPath *spath,
               PathToAlnMap &path2aln_map )
{
    PathId pid = (PathId)path2aln_map.size();
    //ReadIdArray path_rids = ReadIdArray(spath->getReads(), spath->getReads() + spath->getReadCount() );
    path2aln_map.insert( std::pair<PathId, SpaPath*>( pid, spath ) );
    return pid;
}

void assem::printPathSummary( std::string consensus,
                       size_t kmer_depth,
                       KmerId start_kmer,
                       int nparents,
                       int nchildren,
                       int lstop,
                       int rstop,
                       ScoreSummary &coverage,
                       ScoreSummary &entropy,
                       int creads,
					   Param &param
                       )
{
    std::cout << "Start(" << kmer_depth << "|" << nparents << "|" << nchildren << "):" 
              << alpha::IntegerToAminoAcid(start_kmer, param.kmer_size ) << "\t"
              << "Stop(" << conditions[lstop] << "|" << conditions[rstop] << ")\t"
              << "Coverage(" << coverage.sum << "|" << coverage.max << "|" << coverage.min << "|" << coverage.avg << "|" << coverage.med << ")\t"
              << "Entropy(" << entropy.sum << "|" << entropy.max << "|" << entropy.min << "|" << entropy.avg << "|" << entropy.med << ")\t"
              << "Consensus" 
              << "("
              << consensus.size() 
              << "|"
              << creads
              << "):"
              << consensus
              << "\n";
    
}

void assem::adjustReadIds ( PathType &max_path, 
                     DeBruijnGraph &graph,
                     ReadIdArray &path_reads,
                     InvertedIndex &iindex, //std::tr1::unordered_map<KmerId, Read*> &kmer2read_map,
                     VertexToKmerMap &vertex_map, 
                     int k )
{
    if ( max_path.size() <= 1 ) return;

    std::sort( path_reads.begin(), path_reads.end() );

    Vertex source = NULL;
    Vertex target = NULL;

    ReadIdArray rids, diff, scomm, tcomm, stcomm;
    for ( PathType::iterator it = max_path.begin(); it != max_path.end(); ++it ) {
        target = *it;
        KmerId kid = vertex_map[*it];
        
        ReadIdArray rids = ReadIdArray( iindex.getValue(kid)->rid, iindex.getValue(kid)->rid + iindex.getValue(kid)->size);
        std::sort(rids.begin(), rids.end());
        tcomm = Set::Intersect<ReadId>( path_reads, rids, true );
//         std::cout << "\tKmer:" << alpha::IntegerToAminoAcid(vertex_map[*it], k) 
//                   << "\tRead:" << iindex.getValue(kid)->size
//                   << "\tShare:" << tcomm.size();

        diff = Set::Difference<ReadId>( rids, tcomm, true );
        updateReadIds( diff, kid, iindex );

//         std::cout << "\tDiff:" << iindex.getValue(kid)->size << "\n";

        if ( source != NULL ) {
            stcomm = Set::Intersect<ReadId>( scomm, tcomm, true );
            Edge e = edge(source, target, graph).first;
            graph[e].weight -= stcomm.size();
        }
        source = target;
        scomm = tcomm;
    }
}

void assem::updateReadIds ( ReadIdArray &nrids,
                     KmerId &kid,
                     InvertedIndex &iindex ) //std::tr1::unordered_map<KmerId, Read*> &kmer2read_map )
{
    iindex.update(kid, &nrids[0], nrids.size());
}

void assem::dropWeakNodes( NodeSet &deleted_nodes,
                    PathType &max_path,
                    DeBruijnGraph &graph,
                    AdjacencyMap &ancestors, 
                    InvertedIndex &iindex, //std::tr1::unordered_map<KmerId, Read*> &kmer2read_map,
                    VertexToKmerMap &vertex_map,
					Param &param )
//                    size_t min_support )
{
    if ( param.verbose ) std::cout << "\n";
    for ( PathType::iterator it = max_path.begin(); it != max_path.end(); ++it ) {
        //if ( iindex.getValue(vertex_map[*it])->size > min_support ) continue;
        if ( iindex.getValue(vertex_map[*it])->size > (size_t)param.min_depth ) continue;

        NodeArray succs = graph::successors(graph, *it);
        for ( NodeArray::iterator jt = succs.begin(); jt != succs.end(); ++jt ) {
            if ( edge( *it, *jt, graph ).second ) 
                remove_edge( *it, *jt, graph );
            else {
                if ( param.verbose ) 
                    std::cout << "Edge from " 
                              << alpha::IntegerToAminoAcid(vertex_map[*jt], param.kmer_size) << " to "
                              << alpha::IntegerToAminoAcid(vertex_map[*it], param.kmer_size) << " not exist\n";
            }
            graph::dropAncestor(ancestors, *jt, *it);         
        }
        NodeArray preds = ancestors[*it];
        for ( NodeArray::iterator jt = preds.begin(); jt != preds.end(); ++jt ) 
            if ( edge(*jt, *it, graph).second )
                remove_edge( *jt, *it, graph ); 
            else {
                if ( param.verbose ) 
                    std::cout << "Edge from " 
                              << alpha::IntegerToAminoAcid(vertex_map[*jt], param.kmer_size) << " to "
                              << alpha::IntegerToAminoAcid(vertex_map[*it], param.kmer_size) << " not exist\n";
            }
        ancestors.erase(ancestors.find(*it));

        deleted_nodes.insert(*it);
//         if ( param.verbose ) 
//             std::cout << "\tDeleted:" 
//                       << alpha::IntegerToAminoAcid(vertex_map[*it], param.kmer_size) 
//                       << "\t"
//                       << iindex.getValue(vertex_map[*it])->size//kmer2read_map[vertex_map[*it]]->size
//                       << "\n";
    }
}

void assem::updatePathSearchPriority( Vertex v,
                                      std::multimap<double, Vertex> &cov_map,
                                      //std::tr1::unordered_map<Vertex, double> &cov_map,
                               DeBruijnGraph &graph,
                               InvertedIndex &iindex,
                               VertexToKmerMap &vertex_map,
                               AdjacencyMap &ancestors,
							   Param &param )
{    
    size_t sd = ancestors[v].size() +  graph::successors(graph, v).size();
    if ( sd == 0 ) return;
//     if ( sd == 0 ) {
//         cov_map[v] = -1; return;
//     }
    int nread = iindex.getValue(vertex_map[v])->size;
    if ( nread < param.min_seed ) return;
//     if ( nread < param.min_seed ) {
//         cov_map[v] = -1; return;
//     }
    double weight = nread / exp( (double)sd );
    cov_map.insert(std::pair<double, Vertex>(weight, v));
    //cov_map[v] = weight;
}

void assem::displayPathSearchSummary( int nstart, int delnodes, int npaths, int reads )
{
    std::cout << "# Start:" << nstart << "\n";
    std::cout << "# Paths:" << npaths << "\n";
    std::cout << "# Reads:" << reads << "\n";
}



//PathLengthMap assem::getPathLengths( PathToAlnMap &pmap, PathIdSet &path_ids )
PathLengthMap assem::getPathLengths( PathToAlnMap &pmap )
{
    PathLengthMap lmap;
    for ( PathToAlnMap::iterator it = pmap.begin(); it != pmap.end(); ++it ) {
        PathId pid  = it->first;
        SpaPath *spath = it->second;
        //size_t len  = spath->getKmerCount();
        size_t len  = (spath->getConsensusString()).size();
        //if ( path_ids.find(pid) == path_ids.end() ) continue;
        lmap.insert( std::pair<size_t, PathId>( len, pid ) );
    }
    return lmap;
}

// // current latch
// // start and end kmer must be identical
// // So, if there is any mismatch in a given length (default:10), it is ignored.
// int assem::getMinKmerCount( int seqlen, Param &param  )
// {
//     int latch_k = param.latch_length - param.kmer_size + 1;

//     return latch_k;
// }


// /**
//  * Allow mismatch in latch region
//  */
// int assem::getMinKmerCount( int seqlen, Param &param  )
// {
//     int latch_k = int (param.latch_length * param.filter_score) - param.kmer_size + 1;
//     assert(latch_k > 0);
//     return latch_k;
// }



// void assem::initPosFlags( std::vector<int> &pvec,
//                    int qnkmer, 
//                    int rnkmer,
//                    KmerId *qkmers,
//                    KmerId *rkmers,
//                    std::map<int, std::set<int> > &collide )
// {
//     std::set<int> inits;
//     for ( int i = 0; i < qnkmer; i++ )  {
//         for ( int j = 0; j < rnkmer; j++ )  {
//             if ( qkmers[i] == rkmers[j] ) {
//                 collide[j].insert(i);
//                 if ( pvec[j] == -1 ) pvec[j] = i;
//                 else { 
//                     if ( j == 0 ) inits.insert(i);
//                     //put same kmers in right place
//                     else if ( abs(int(pvec[j]-pvec[j-1])) > abs(int(i-pvec[j-1])) ) pvec[j] = i;
                    
//                 }
//             }
//         }
//     }            
// }

void assem::initPosFlags( std::vector<int> &pvec,
                   int qnkmer, 
                   int rnkmer,
                   KmerId *qkmers,
                   KmerId *rkmers,
                   std::map<int, std::set<int> > &collide )
{
    std::set<int> inits;
    for ( int i = 0; i < qnkmer; i++ )  {
        for ( int j = 0; j < rnkmer; j++ )  {
            if ( qkmers[i] == rkmers[j] ) {
                collide[j].insert(i);
                if ( pvec[j] == -1 ) pvec[j] = i;
                else {
                    if ( j == 0 ) continue;
                    // make tight ragne when same kmers appear multiple time
                    if ( abs(int(pvec[j]-pvec[j-1])) > abs(int(i-pvec[j-1])) ) pvec[j] = i;
                }
            }
        }
    }            
}

/*
bool assem::__latchableLeft( std::vector<int> &pvec, 
                             int &qs, 
                             int &qe,
                             int &ms, 
                             int &me,
                             int direction,
                             Param &pram )
{
    if ( ms >= qs ) return false;
    
}

bool assem::latchableLeft( std::vector<int> &pvec, 
                           IntPair &qrange, 
                           IntPair &rrange,
                           int direction,
                           Param &pram )
{
    if ( param.pair_flag && direction == NEGATIVE )
        __latchablePairEndLeft( pvec, qrange, rrange, param );
//     if ( param.pair_flag && direction == NEGATIVE )
//         __latchablePairEndRight( pvec, qrange, rrange, param );


    int i = 0;
    while ( pvec[i] == -1 ) i++;
    if ( i == pvec.size() ) return false;

    
}
*/

bool assem::latchable( std::map<int, std::set<int> > &collide,
                std::vector<int> &pvec,
                int qnkmer, 
                int rnkmer,
                int &count,
                int min_length,
                int direction,
				Param &param  )
{
    int start = 0;
    if ( direction == 1 ) start = rnkmer-1; // latch right
    int end = start;

    if ( collide.find(start) == collide.end() ) return false;
    
    if ( direction == -1 ) {
        while ( end < (int)pvec.size() && pvec[end] != -1 ) end++;
        int len = (end-start) + param.kmer_size-1;
        if ( len < min_length ) return false;

        for ( int i = end-1; i >= start; i-- ) {
            int srch = (qnkmer-1)-count;
            if ( collide[i].find(srch) == collide[i].end() ) return false;
            count++;
        }
    } else {
        while ( end >= 0 && pvec[end] != -1 ) end--;
        int len = (start-end) + param.kmer_size-1;
        if ( len < min_length ) return false;

        for ( int i = end+1; i <= start; i++ ) {
            int srch = count;
            if ( collide[i].find(srch) == collide[i].end() ) return false;
            count++;
        }
    }

    return true;
}



std::vector<int> assem::updatePosFlags(std::vector<int> &pvec, int count, int rstart, int qstart )
{
    std::vector<int> nvec = pvec;
    for ( int i = 0; i < count; i++ ) 
        nvec[rstart+i] = qstart+i;
    return nvec;
}

bool assem::goodLatch(std::vector<int> &pvec, KmerId *qkmers, KmerId *rkmers, int qnkmer, int rnkmer, int direction, Param &param  )
{
    int s,e;
    if ( direction == -1 ) {
        s = 0;
        e = s+1;
        while ( e < (int)pvec.size() && pvec[e] != -1 ) { e++; }
        e--;
    } else {
        e = (int)pvec.size()-1;
        s = e-1;
        while ( s >= 0 && pvec[s] != -1 ) { s--; }
        s++;
    }

    if ( pvec[s] + (e-s+1) > qnkmer ) return 0;
    if ( s + (e-s+1) > rnkmer ) return 0;


    std::string query = biostr::getSequenceString( qkmers+pvec[s], e-s+1, param.kmer_size );
    std::string sbjct = biostr::getSequenceString( rkmers+s, e-s+1, param.kmer_size );

    if ( query == sbjct ) return 1;
    return 0;
}

                
std::vector<int> assem::makePosFlags( int qnkmer, 
                               int rnkmer, 
                               KmerId *qkmers,
                               KmerId *rkmers,
                               int min_length, Param &param)
{
    std::vector<int> pvec(rnkmer, -1);
    std::map<int, std::set<int> > collide;
    initPosFlags( pvec, qnkmer, rnkmer, qkmers, rkmers, collide );
    
    int count = 0;
    if ( latchable(collide, pvec, qnkmer, rnkmer, count, min_length, -1, param) ) { // latch left
        if ( qnkmer-count>=0) {
            std::vector<int> nvec = updatePosFlags(pvec, count, 0, qnkmer-count);
            if ( goodLatch(nvec, qkmers, rkmers, qnkmer, rnkmer, -1, param) ) pvec = nvec;
        }
    } else if ( latchable(collide, pvec, qnkmer, rnkmer, count, min_length, 1, param) ) { // latch right
        if ( rnkmer-count>=0) {
            std::vector<int> nvec = updatePosFlags(pvec, count, rnkmer-count, 0);
            if ( goodLatch(nvec, qkmers, rkmers, qnkmer, rnkmer, 1, param) ) pvec = nvec;
        }
    }

    return pvec;
}

void assem::setKmerMap( std::tr1::unordered_map<KmerId, int> & count_map, 
                        KmerArray &kmers ) 
{
    count_map.clear();
    for ( size_t i = 0; i < kmers.size(); i++ ) {
        if (count_map.find(kmers[i]) == count_map.end() )
            count_map[kmers[i]] = 0;
        count_map[kmers[i]]++;
    }
}

size_t assem::countSameKmers( std::tr1::unordered_map<KmerId, int> &count_map,
                              KmerArray &query_kmers )
{
    //double t0 = mytime();
    std::tr1::unordered_map<KmerId, int>::iterator it;

    size_t count = 0;
    for ( size_t i = 0; i < query_kmers.size(); i++ ) {
        it = count_map.find(query_kmers[i]);
        if ( it == count_map.end() ) continue;
        if ( it->second > 0 ) {
            count++;
            it->second--;
        }
    }
    //std::cout << "map time w/ count map:" << mytime()-t0 << "\n";
    return count;
}

size_t assem::countSameKmers( std::string &sbjct,
                              std::string &query,
                              size_t small_k )
{
    //double t0 = mytime();
    KmerArray query_kmers = biostr::getKmers( query, small_k );
    KmerArray sbjct_kmers = biostr::getKmers( sbjct, small_k );
    
    std::tr1::unordered_map<KmerId, int> count_map;
    for ( size_t i = 0; i < query_kmers.size(); i++ ) {
        if (count_map.find(query_kmers[i]) == count_map.end() )
            count_map[query_kmers[i]] = 0;
        count_map[query_kmers[i]]++;
    }

    
    std::tr1::unordered_map<KmerId, int>::iterator it;

    size_t count = 0;
    for ( size_t i = 0; i < sbjct_kmers.size(); i++ ) {
        it = count_map.find(sbjct_kmers[i]);
        if ( it == count_map.end() ) continue;
        if ( it->second > 0 ) {
            count++;
            it->second--;
        }
    }
        
    //std::cout << "map time w/ strings:" << mytime()-t0 << "\n";
    return count;
    //KmerArray comm = Set::Intersect<KmerId>(query_kmers, sbjct_kmers, false );
    //return comm.size();
}


size_t assem::countSameKmers( KmerToPathMap &pathid_map,
                              PathToAlnMap &path2aln_map,
                              PathId query_pid,
                              PathId sbjct_pid,
                              size_t small_k )
{

    std::string query = path2aln_map[query_pid]->getConsensusString();
    std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();
    return countSameKmers( sbjct, query, small_k );

//     KmerArray query_kmers = biostr::getKmers( path2aln_map[query_pid]->getConsensusString(), small_k );
//     KmerArray sbjct_kmers = biostr::getKmers( path2aln_map[sbjct_pid]->getConsensusString(), small_k );

//     size_t count = 0;
//     KmerArray comm = Set::Intersect<KmerId>(query_kmers, sbjct_kmers, false );
//     return comm.size();
}


// size_t assem::countSameKmers( KmerToPathMap &pathid_map,
//                               PathToAlnMap &path2aln_map,
//                               PathId query_pid,
//                               PathId sbjct_pid )
// {
//     KmerId *kmers = path2aln_map[query_pid]->getKmers();
//     size_t  nkmer = path2aln_map[query_pid]->getKmerCount();    

//     size_t count = 0;
//     KmerSet kset = KmerSet( kmers, kmers+nkmer );
//     for ( KmerSet::iterator kit = kset.begin(); kit != kset.end(); ++kit ) {
//         if ( pathid_map.find(*kit) == pathid_map.end() ) continue;
//         std::set<PathId>::iterator pit;
//         for ( pit = pathid_map[*kit].begin(); pit != pathid_map[*kit].end(); ++pit ) 
//             if ( *pit == sbjct_pid ) count++;
//     }
//     return count;
// }

std::map<PathId, int>  assem::getKmerCountMap( PathIdSet skip_pids,
                                        KmerId *kmers,
                                        int nkmer,
                                        KmerToPathMap &pathid_map,
                                        PathIdSet &merged_paths )
{
    std::map<PathId, int> count_map;

    KmerSet kset = KmerSet( kmers, kmers+nkmer );
    for ( KmerSet::iterator kit = kset.begin(); kit != kset.end(); ++kit ) {
        if ( pathid_map.find(*kit) == pathid_map.end() ) continue;
        std::set<PathId>::iterator pit;
        for ( pit = pathid_map[*kit].begin(); pit != pathid_map[*kit].end(); ++pit ) {
            if ( skip_pids.find(*pit) != skip_pids.end() ) continue; 
            if ( merged_paths.find(*pit) != merged_paths.end() ) continue; 
            if ( count_map.find(*pit) == count_map.end() )
                count_map.insert( std::pair<PathId, int>(*pit, 0) );
            count_map[*pit]++;
        }
    }
    return count_map;
}

PosPathPairList assem::getKmerPosList( std::map<PathId, int> &count_map,
                                       KmerId *kmers,
                                       int nkmer,
                                       PathToAlnMap &path2aln_map,
                                       int mink, 
                                       bool is_query,
                                       Param &param )
{
    PosPathPairList pp_list;
    std::multimap<int, PathId> icount_map = util::invert_map<PathId, int>(count_map);
    std::multimap<int, PathId>::reverse_iterator rt = icount_map.rbegin();
    for ( ; rt != icount_map.rend(); ++rt ) {
        if ( rt->first < mink ) continue;
        SpaPath *spath = path2aln_map[rt->second];
        KmerId *mkmers = spath->getKmers();
        size_t nmkmer = spath->getKmerCount();    

        std::vector<int> pvec;
        if ( is_query ) pvec = makePosFlags(nkmer, nmkmer, kmers, mkmers, param.latch_length, param);
        else pvec = makePosFlags(nmkmer, nkmer, mkmers, kmers, param.latch_length, param);
        pp_list.push_back(PosPathPair(pvec, rt->second));
    }
    return pp_list;
}

PathIdSet assem::searchOverlapCandidates(  //PathIdSet pids,
                                          PathId sbjct_pid,
                                          KmerId *skmers,
                                          //int snkmer,
                                          int srch_nkmer,
                                          PathToAlnMap &path2aln_map,
                                          KmerToPathMap &pathid_map,
                                          PathIdSet &merged_paths,
                                          int mink,
                                          Param &param,
                                          int direction )
{
    PathIdSet avoids;
    avoids.insert(sbjct_pid);
    //PosPathPairList pos_paths = searchSimilarPaths(avoids, skmers, snkmer, path2aln_map, pathid_map, merged_paths, mink, false, param);
    PosPathPairList pos_paths = searchSimilarPaths(avoids, skmers, srch_nkmer, path2aln_map, pathid_map, merged_paths, mink, false, param);

    PathIdList candidates;
    for ( PosPathPairList::iterator it = pos_paths.begin(); it != pos_paths.end(); ++it ) {
        PathId query_pid = it->second;
        std::vector<int> pvec = it->first;
        if ( param.verbose ) {
            std::cout << "query:" << query_pid << "\n";
            printPositions(pvec);
        }
        
        if ( direction == LEFT ) {
            for ( int i = 0; i < srch_nkmer; i++ ) {
                if ( pvec[i] != -1 ) candidates.push_back(query_pid);
                break;
            }
        } else {
            int i;
            for ( i = srch_nkmer-1; i>=0; i-- ) {
                if ( pvec[i] != -1 ) candidates.push_back(query_pid);
                break;
            }
        }
    }

    return PathIdSet( candidates.begin(), candidates.end() );
}
    
//     std::map<PathId, int> count_map = getKmerCountMap( avoids, skmers, snkmer, pathid_map, merged_paths );
//     if ( count_map.size() == 0 ) return PosPathPairList();

//     //size_t srch_nkmer = param.latch_length - param.kmer_size + 1;                    
//     std::tr1::unordered_map<PathId, bool> overlap_paths;

//     for ( std::map<PathId, int>::iterator it = count_map.begin(); it != count_map.end(); ++it ) {
//         PathId query_pid = it->first;
//         KmerId *qkmers = path2aln_map[query_pid]->getKmers();
//         size_t  qnkmer = path2aln_map[query_pid]->getKmerCount();
//         //if ( qnkmer < srch_nkmer ) continue;
            
//         int beg, end;        
//         // sbjct:     oooooxxxxx
//         // query:xxxxxooooo
//         if ( direction == LEFT ) {
//             //beg = qnkmer-srch_nkmer;
//             beg = qnkmer/2;
//             end = qnkmer;
//         } 
//         // sbjct:xxxxxooooo
//         // query:     oooooxxxxx
//         else {
//             beg = 0;
//             //end = srch_nkmer;
//             end = qnkmer/2;
//         }
        
//         for ( int i = 0; i < snkmer; i++ ) {
//             for ( int j = beg; j < end; j++ ) {
//                 if ( skmers[i] == qkmers[j] ) {
//                     overlap_paths.insert( std::pair<PathId, bool>( query_pid, true ) );
//                     break;
//                 }
//             }
//         }
//     }
    
    
//     for ( std::map<PathId, int>::iterator it = count_map.begin(); it != count_map.end();  ) {
//         if ( overlap_paths.find( it->first ) != overlap_paths.end() &&
//              overlap_paths[it->first] == true ) {
//             ++it;
//         } else {
//             count_map.erase(it++);
//         }
//     }
//     if ( count_map.size() == 0 ) return PosPathPairList();


//     return getKmerPosList( count_map, skmers, snkmer, path2aln_map, mink, param );
//}

PosPathPairList assem::searchSimilarPaths( PathIdSet pids,
                                           KmerId *kmers,
                                           int nkmer,
                                           PathToAlnMap &path2aln_map,
                                           KmerToPathMap &pathid_map,
                                           PathIdSet &merged_paths,
                                           int mink,
                                           bool is_query,
                                           Param &param )
{
    if ( nkmer == 0 ) return PosPathPairList();
    
    std::map<PathId, int> count_map = getKmerCountMap( pids, kmers, nkmer, pathid_map, merged_paths );
    if ( count_map.size() == 0 ) return PosPathPairList();
    
    return getKmerPosList( count_map, kmers, nkmer, path2aln_map, mink, is_query, param );
}



void assem::addPathIds( PathId pid,
                    KmerId *kmers,
                    size_t nkmer, 
                    KmerToPathMap &pathid_map)
{
    if ( nkmer == 0 ) return;
    for ( size_t i = 0; i < nkmer; i++ )
        pathid_map[kmers[i]].insert(pid);
}

void assem::dropPathIds( PathId pid,
                  KmerId *kmers,
                  size_t nkmer, 
                  KmerToPathMap &pathid_map)
{
    if ( nkmer == 0 ) return;
    for ( size_t i = 0; i < nkmer; i++ )
        pathid_map[kmers[i]].erase(pid);
}


// bool assem::mergePathPair( PathId query_pid, 
//                     PathId sbjct_pid, 
//                     AlignSummary &summary,
//                     PathIdSet &merged_paths,
//                     PathToAlnMap &path2aln_map,
//                            //KmerToPathMap &pathid_map,
//                     BitString *bstrs,
//                     char *strands, 
//                     ReadId *pairs,
//                     PathId *used_reads,
//                     InvertedIndex &iindex,
// 					Param &param )
// {
//     SpaPath *sbjct_spath = path2aln_map[sbjct_pid];
//     SpaPath *query_spath = path2aln_map[query_pid];

//     if ( param.verbose ) {
//         std::cout << "Merging pairs:\tQuery:" << query_pid << "\tSbjct:" << sbjct_pid << "\n";
//         std::cout << "old str:" << biostr::getSequenceString(sbjct_spath->getKmers(), sbjct_spath->getKmerCount(), param.kmer_size) << "\n";
//         std::cout << "old nread:" << sbjct_spath->getReadCount() << "\n";
//     }


// //     dropPathIds( query_pid, query_spath->getKmers(), query_spath->getKmerCount(), pathid_map );
// //     dropPathIds( sbjct_pid, sbjct_spath->getKmers(), sbjct_spath->getKmerCount(), pathid_map );

//     if ( param.verbose ) 
//         std::cout << "Joining path\t"
//                   << "start:" << summary.range.first << "\t"
//                   << "lgap:" << summary.lgap.first << "\t"
//                   << "tgap:" << summary.egap.first << "\n";
    
//     if ( summary.ilist.size() > 0 || summary.dlist.size() > 0 ) {
//         if ( param.verbose ) std::cout << "Gappy joining\n";
//     }
   
//     sbjct_spath->join(query_spath, summary.range.first, summary.lgap.first, summary.egap.first, summary.ilist, summary.dlist, param.kmer_size, iindex, bstrs, param);
// //     addPathIds ( sbjct_pid, sbjct_spath->getKmers(), sbjct_spath->getKmerCount(), pathid_map );

// //     if ( ! success ) {
// //         if ( param.verbose ) std::cout << "Join failed\n";
// //         return false;
// //     }

//     if ( param.verbose ) {
//         std::cout << "new str:" << biostr::getSequenceString(sbjct_spath->getKmers(), sbjct_spath->getKmerCount(), param.kmer_size) << "\n";
//         std::cout << "new nread:" << sbjct_spath->getReadCount() << "\n";
//     }

//     int error = 0;

// 	MSA msa(*sbjct_spath, sbjct_pid, used_reads, bstrs, strands, pairs, PILE|FLAT|TRIM|CODON, error, param);
// //     if ( error ) {
// //         std::cout << "\n** MSA failed\n"; 
// //         delete spath; continue;
// //     }

//     if ( error > 0 ) {
//         if ( param.verbose ) std::cout << "MSA errors\n";
//         return false;
//     }

//     if ( param.verbose ) {
//         std::cout << "\nAlignment\n";
//         msa.printAlignment(std::cout, *sbjct_spath, 100);
//         std::cout << "\nProfile\n";
//         msa.printProfile(std::cout);
//     }



//     merged_paths.insert(query_pid);
//     return true;
// }

// IntPair assem::__getLatchRangeLeft(std::vector<int>& pos_vec, Param &param )
// {
//     int i = 0;
//     while ( pos_vec[i] == -1 ) {
//         if ( i == pos_vec.size() ) return IntPair(-1,-1);
//         i++;
//     }

//     int s = i;
//     int prev = pos_vec[s];
//     int curr = pos_vec[s];
//     for ( i = s+1; i < pos_vec.size(); i++ ) {
//         curr = pos_vec[i];
//         if ( curr != -1 && curr != prev+1 ) return IntPair(s, i-1);
//         if ( curr == -1 ) prev = prev + 1;
//         else prev = curr;
//     }
//     return IntPair(s, i);
// }

// IntPair assem::__getLatchRangeRight(std::vector<int>& pos_vec, Param &param )
// {
//     int i = pos_vec.size()-1;
//     int min = 1000000;
//     while ( pos_vec[i] != -1 && pos_vec[i] < min ) {
//         if ( i < 0 ) return IntPair(-1,-1);
//         if ( pos_vec[i] < min ) min = pos_vec[i];
//         else break;
//         i--;
//     }
//     int s = i;
//     int prev = pos_vec[s];
//     int curr = pos_vec[s];
//     for ( i = s+1; i < pos_vec.size(); i++ ) {
//         curr = pos_vec[i];
//         if ( curr != prev+1 ) return IntPair(s, i-1);
//         prev = curr;
//     }
//     return IntPair(s, i);
// }

std::vector<IntPair> assem::makeBlocks( std::vector<int> &pos_vec )
{
    std::vector<IntPair> blocks;

    int i = 0;
    while ( true ) {
        if ( i >= (int)pos_vec.size() ) break;
        
        int count = 0;
        int start = i;
        int prev  = -1;
        while ( i < (int)pos_vec.size() && pos_vec[i] != -1 ) { 
            if ( prev > pos_vec[i] ){
                if ( count > 0 ) blocks.push_back(IntPair(count, start));
                count = 0; start = i;
            }
            prev = pos_vec[i];
            count++; i++;
        }
        if ( count > 0 ) blocks.push_back(IntPair(count, start));
        while ( i < (int)pos_vec.size() && pos_vec[i] == -1 ) { i++; }
    }
    return blocks;
}

IntPair assem::expandBound( IntPair bound, 
                            std::vector<int> &pos_vec,
                            Param &param)
{
    int size = bound.first;
    int spos = bound.second;
    int epos = spos + size - 1;
    
    //if ( param.verbose ) std::cout << "os:" << spos << ":" << pos_vec[spos] << "\toe:" << epos << ":" << pos_vec[epos] << "\n";
    
    int nbad = 0;
    int nepos = epos;
    int dist = 1;
    int curr;
    int i = epos+1;
    while ( i < (int)pos_vec.size() ) {
            curr = pos_vec[i];
            if ( curr != -1 ) {
                if ( curr <= pos_vec[nepos] ) {
                    if ( nbad ) break;
                    nbad++;
                } else {
                    // tight bound 
                    // may need to give some allowing distances
                    if ( curr > dist+pos_vec[nepos] ) break;
                    else nepos = i;
                }
            } else nbad = 0;
            ++i; ++dist;
    }
    
    nbad = 0;
    int nspos = spos;
    dist = 1;
    i = spos-1;
    while ( i >= 0 ) {
        curr = pos_vec[i];
        if ( curr != -1 ) {
            if ( curr >= pos_vec[nspos] ) {
                if ( nbad ) break;
                nbad++;
            } else {
                if ( curr < pos_vec[nspos]-dist ) break;
                    else nspos = i;
            }
        } else nbad = 0;
        --i; ++dist;
    }
    
    //if ( param.verbose ) std::cout << "ns:" << nspos << ":" << pos_vec[nspos] << "\tne:" << nepos << ":" << pos_vec[nepos] << "\n";
    return IntPair(nspos, nepos);
}

bool assem::inorder(std::vector<int> &pos)
{
    int prev = pos[0];
    for ( size_t i = 1; i < pos.size(); i++ ) {
        int curr = pos[i];
        if ( curr == prev+1 ) 
            prev = curr;
        else 
            return false;
    }
    return true;
}

void assem::makeTightBound(IntPair &block, std::vector<int> &pos_vec)
{
    int start = block.second;
    int count = block.first;
    

    std::vector<int> sub = std::vector<int>( pos_vec.begin()+start, pos_vec.begin()+start+count );
    while ( sub.size() ) {
        if ( inorder(sub) ) return;
        std::vector<int> lsub = std::vector<int>( sub.begin(), sub.begin()+sub.size()-1 );
        if ( inorder(lsub) )  { 
            block.first--;
            break;
        }
        std::vector<int> rsub = std::vector<int>( sub.begin()+1, sub.end());
        if ( inorder(rsub) )  {
            block.second++;
            block.first--;
            break;
        }
        sub.pop_back();
        if ( sub.size() ) sub.erase(sub.begin());
    }
}


IntPair assem::getMatchRangeLeft(std::vector<int>& pos_vec, Param &param )
{
    std::vector<IntPair> blocks = makeBlocks( pos_vec );
    if ( blocks.size() == 0 ) return IntPair(-1,-1);

    if ( param.verbose ) 
        std::cout << "Before:" << blocks[0].first << "\t" << blocks[0].second << "\n";
    makeTightBound(blocks[0], pos_vec);
    if ( param.verbose ) 
        std::cout << "After:" << blocks[0].first << "\t" << blocks[0].second << "\n";
    return expandBound(blocks[0], pos_vec, param );
}

IntPair assem::getMatchRangeRight(std::vector<int>& pos_vec, Param &param )
{
    std::vector<IntPair> blocks = makeBlocks( pos_vec );
    if ( blocks.size() == 0 ) return IntPair(-1,-1);
    
    if ( param.verbose ) 
        std::cout << "Before:" << blocks.back().first << "\t" << blocks.back().second << "\n";
    makeTightBound(blocks.back(), pos_vec);
    if ( param.verbose ) 
        std::cout << "After:" << blocks.back().first << "\t" << blocks.back().second << "\n";
    return expandBound(blocks.back(), pos_vec, param);
}


IntPair assem::getMatchRangeMax(std::vector<int>& pos_vec, Param &param )
{
    std::vector<IntPair> blocks = makeBlocks( pos_vec );
    if ( blocks.size() == 0 ) return IntPair(-1,-1);

    std::map<int, int> block_size;
    for ( size_t i = 0; i < blocks.size(); i++ ) 
        block_size.insert(blocks[i]);
    std::map<int, int>::reverse_iterator it =  block_size.rbegin();
    
    return expandBound(IntPair(it->first, it->second), pos_vec, param);
}




// -1 0 1 2 3 4 5 6 7 8 9 10 11 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 -1 -1 -1 -1 -1 -1 63 64 

// to do 
// handle this
/*
Query:PEQQVDPRKAAVEAAIARAKARKREQQPANAEPEEQVDPRKAAV
Sbjct:DPRKAAVEAAIARAKARKLEQQQANAVPEEQVDPRKAAVEAAIARAKARKLEQQQANAEPEQQVDPRKAAVEAAIARAKARKLEQQQANAVP
pos_vec size:87
5 6 7 8 9 10 11 12 13 14 15 16 17 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 32 33 34 35 36 37 38 7 8 9 10 11 12 13 14 15 16 17 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 28 -1 -1 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 
os:59:0 oe:76:17
ns:59:0 ne:76:17
*/
IntPair assem::getMatchRange(std::vector<int>& pos_vec, Param &param )
{
    std::map<int, int> block_size;
    
    int i = 0;
    while ( true ) {
        if ( i >= (int)pos_vec.size() ) break;
        
        int count = 0;
        int start = i;
        int prev  = -1;
        while ( i < (int)pos_vec.size() && pos_vec[i] != -1 ) { 
            if ( prev > pos_vec[i] ){
                if ( count > 0 ) block_size.insert(IntPair(count, start));
                count = 0; start = i;
            }
            prev = pos_vec[i];
            count++; i++;
        }
        if ( count > 0 ) block_size.insert(IntPair(count, start));
        while ( i < (int)pos_vec.size() && pos_vec[i] == -1 ) { i++; }
    }
    i = 0;
    if ( param.verbose ) {
        std::cout << "Position vector:\n";
        while ( i < (int)pos_vec.size() ) 
            std::cout << pos_vec[i++] << " ";
        std::cout << "\n";
    }
    if ( block_size.size() == 0 ) return IntPair(-1,-1);

    std::map<int, int>::reverse_iterator it =  block_size.rbegin();
    int size = it->first;
    int spos = it->second;
    int epos = spos + size - 1;
    
    if ( param.verbose ) std::cout << "os:" << spos << ":" << pos_vec[spos] << "\toe:" << epos << ":" << pos_vec[epos] << "\n";
    
    int nbad = 0;
    int nepos = epos;
    int dist = 1;
    int curr;
    i = epos+1;
    while ( i < (int)pos_vec.size() ) {
            curr = pos_vec[i];
            if ( curr != -1 ) {
                if ( curr <= pos_vec[nepos] ) {
                    if ( nbad ) break;
                    nbad++;
                } else {
                    // tight bound 
                    // may need to give some allowing distances
                    if ( curr > dist+pos_vec[nepos] ) break;
                    else nepos = i;
                }
            } else nbad = 0;
            ++i; ++dist;
    }
    
    nbad = 0;
    int nspos = spos;
    dist = 1;
    i = spos-1;
    while ( i >= 0 ) {
        curr = pos_vec[i];
        if ( curr != -1 ) {
            if ( curr >= pos_vec[nspos] ) {
                if ( nbad ) break;
                nbad++;
            } else {
                if ( curr < pos_vec[nspos]-dist ) break;
                    else nspos = i;
            }
        } else nbad = 0;
        --i; ++dist;
    }
    
    if ( param.verbose ) std::cout << "ns:" << nspos << ":" << pos_vec[nspos] << "\tne:" << nepos << ":" << pos_vec[nepos] << "\n";
    return IntPair(nspos, nepos);
}

// int countMismatchKmers(std::vector<int>& pos_vec,
//                        size_t s,
//                        size_t e )
// {
//     int count = 0;
//     for ( size_t i = s; i <= e; i++ )
//         if ( pos_vec[i] == -1 ) count++;
//     return count;
// }

int assem::countMatchKmers(std::vector<int>& pos_vec,
                     size_t s,
                     size_t e )
{
    int count = 0;
    for ( size_t i = s; i <= e; i++ )
        if ( pos_vec[i] != -1 ) count++;
    return count;
}


bool assem::indeled( std::vector<int>& pos_vec,
             IntPair &range )
{
    bool indel_flag = true;

    size_t s = range.first;
    size_t e = range.second;
    size_t min = pos_vec[s];
    for ( size_t i = s+1; i <= e; i++ ) {
        if ( pos_vec[i] == -1 ) continue;
        if ( (int)min > pos_vec[i] ) min = pos_vec[i];
    }
    size_t max = pos_vec[e];
    for ( int i = e-1; i >= (int)s; i-- ) {
        if ( pos_vec[i] == -1 ) continue;
        if ( (int)max < pos_vec[i] ) max = pos_vec[i];
    }

    if ( (e-s) == (max-min) ) indel_flag = false;
    
    return indel_flag;
}

              
bool assem::ordered(std::vector<int>& pos_vec,
             IntPair &range )
{
    bool order_flag = true;

    size_t s = range.first;
    size_t e = range.second;

    int prev = -1;
    int curr = -1;
    for ( size_t i = s+1; i <= e; i++ ) {
        curr = pos_vec[i];
        if ( curr != -1 && curr < prev ) { order_flag = false; break; }
        prev = curr;
    }
    return order_flag;
}

// bool assem::latchableReadRight( size_t qs,
//                       size_t qe,
//                       size_t ms,
//                       size_t me,
//                       size_t qnkmer,
//                       size_t rnkmer,
//                       KmerId *qkmers,
//                       KmerId *rkmers,
//                       size_t spur_offset,
//                       Param &param)
// {
//     int query_len = qnkmer + param.kmer_size - 1;
//     int sbjct_len = me - ms + param.kmer_size;

//     // xxxxxxxxxxoooooooooo
//     // ----------oooooooooo
//     if ( qs == 0 ) {
//         // 1. end of reference sequence
//         // xxxxxxxxxxoooooooooo----------
//         // ----------ooooooooooxxxxxxxxxx
//         //if ( me == pos_vec.size()-1 ) {
//         if ( me == rnkmer-1 ) {
//             if ( param.verbose ) std::cout << "Type:Read Latch Right\n";       
//             if ( alpha::getLastAminoAcid(rkmers[rnkmer-1], param.kmer_size) == '*' ) {
//                 if (param.verbose) std::cout << "Invalid Latch - stop codon\n";
//                 return false;
//             }
//             else if ( sbjct_len < (int)latch_offset  ) {
//                 if (param.verbose) std::cout << "BAD LATCH - SHORT:" << sbjct_len << "\n";
//                 return false;
//             }
//             else { 
//                 if (param.verbose) std::cout << "TRUE LATCH\n"; 
//                 return true;
//             }
//         }
//     }
//     return false;
// }

// bool assem::latchableReadLeft( size_t qs,
//                       size_t qe,
//                       size_t ms,
//                       size_t me,
//                       size_t qnkmer,
//                       size_t rnkmer,
//                       KmerId *qkmers,
//                       KmerId *rkmers,
//                       size_t spur_offset,
//                       Param &param)
// {
//     int query_len = qnkmer + param.kmer_size - 1;
//     int sbjct_len = me - ms + param.kmer_size;

//     if ( qe == qnkmer-1 ) {
//         // 1. beginning of reference sequence
//         // ----------ooooooooooxxxxxxxxxx
//         // xxxxxxxxxxoooooooooo----------
//         if ( ms == 0 ) {
//             if ( param.verbose ) std::cout << "Type:Read Latch Left\n";
//             if ( alpha::getLastAminoAcid(qkmers[qnkmer-1], param.kmer_size) == '*' ) {
//                 if ( me == rnkmer-1 && alpha::getLastAminoAcid(rkmers[rnkmer-1], param.kmer_size) == '*' ) 
//                     return true;
//                 else {
//                     if (param.verbose) std::cout << "Invalid latch - stop codon\n";
//                     return false;
//                 }
//             }
//             else if ( sbjct_len < (int)latch_offset  ) {
//                 if ( param.verbose) std::cout << "BAD LATCH - SHORT:" << sbjct_len << "\n";
//                 return false;
//             }
//             else { 
//                 if (param.verbose) std::cout << "TRUE LATCH\n";
//                 return true;
//             }
//         }
//     }
//     return false;    
// }

// int assem::getPathReadType(IntPair &qrange,
//                        IntPair &rrange,
//                        KmerId *qkmers,
//                        KmerId *rkmers,
//                        size_t qnkmer,
//                        size_t rnkmer,
//                        size_t latch_offset, 
//                        size_t spur_offset,
//                        Param &param  )

// {    
//     if ( isBubble( qs, qe, ms, me, qnkmer, rnkmer, qkmers, rkmers, spur_offset, type, param ) ) return type;
//     if ( isSpurRight( qs, qe, ms, me, qnkmer, rnkmer, qkmers, rkmers, spur_offset, type, param ) ) return type;
//     if ( isSpurLeft( qs, qe, ms, me, qnkmer, rnkmer, qkmers, rkmers, spur_offset, type, param ) ) return type;
//     if ( isFrayedRope( qs, qe, ms, me, qnkmer, rnkmer, qkmers, rkmers, spur_offset, type, param  ) ) return type;
//     if ( latchableReadRight( qs, qe, ms, me, qnkmer, rnkmer, qkmers, rkmers, spur_offset, param  ) ) return LATCH;
//     if ( latchableReadLeft( qs, qe, ms, me, qnkmer, rnkmer, qkmers, rkmers, spur_offset, param ) ) return LATCH;
    
//     return NOTYPE;
// }

int assem::getLinkType(IntPair &qrange,
                       IntPair &rrange,
                       KmerId *qkmers,
                       KmerId *rkmers,
                       size_t qnkmer,
                       size_t rnkmer,
                       size_t latch_offset, 
                       size_t spur_offset,
                       Param &param  )
{
    int type = NOTYPE;

    /* reference start and end postions */
    size_t ms = rrange.first;
    size_t me = rrange.second;

    /* query start and end positions */
    size_t qs = qrange.first;
    size_t qe = qrange.second;


    int query_len = qnkmer + param.kmer_size - 1;
    int sbjct_len = me - ms + param.kmer_size;


    /* Bubble or substring */
    // xxxxxooooooooooxxxxx
    // -----oooooooooo-----
    if ( qs == 0 && qe == qnkmer-1 ) {
        /* Bubble */
        type = BUBBLE;
        if ( param.verbose ) std::cout << "Bubble\n";

        std::string sbjct = biostr::getSequenceString( &rkmers[ms], me-ms+1, param.kmer_size );
        std::string query = biostr::getSequenceString( &qkmers[qs], qe-qs+1, param.kmer_size );

        if ( sbjct == query ) {
            if ( param.verbose ) std::cout << "Substring\n";
            type = SUBSTRING;
        }
        else if ( alpha::getLastAminoAcid(qkmers[qnkmer-1], param.kmer_size) == '*' ) {
            int offset = rnkmer - me - 1;
            if ( offset > 2*(int)spur_offset ) {
                type = NOTYPE;
                if ( param.verbose ) {
                    std::cout << "Stop condon with long offset:" << offset << "\tNot bubble\n";
                }
            }
        } 
    }
    
    // xxxxxxxxxxoooooooooo
    // ----------oooooooooo
    else if ( qs == 0 ) {
        // 1. end of reference sequence
        // xxxxxxxxxxoooooooooo----------
        // ----------ooooooooooxxxxxxxxxx
        //if ( me == pos_vec.size()-1 ) {
        if ( me == rnkmer-1 ) {
            if ( param.verbose ) std::cout << "Type:Latch Right\n";       
            if ( alpha::getLastAminoAcid(rkmers[rnkmer-1], param.kmer_size) == '*' ) {
                if (param.verbose) std::cout << "Invalid Latch - stop codon\n";
            }
            else if ( sbjct_len < (int)latch_offset  ) {
                if (param.verbose) std::cout << "BAD LATCH - SHORT:" << sbjct_len << "\n";
            }
            else { 
                if (param.verbose) std::cout << "TRUE LATCH\n"; 
                type = LATCH;
            }
        }
        // 2. mismatch in query sequence
        // xxxxxxxxxxooooooooooxxxxxxxxxx
        // ----------ooooooooooxxxxx-----
        else {
            if ( param.verbose ) std::cout << "Type:Spur Right\n";
            if ( (query_len-qe) > 2*spur_offset ) { 
                if ( param.verbose ) std::cout << "BAD SPUR:\tLength:" << query_len << "\tEnd:" << qe << "\n";;
            }
            else if ( (query_len-qe) > spur_offset ) {
                if ( param.verbose ) std::cout << "WEAK SPUR\n";
                type = WEAK_SPUR;
            } else {
                if ( param.verbose ) std::cout << "TRUE SPUR\n";
                type = SPUR;
            }
            if ( type != NOTYPE && alpha::getLastAminoAcid(qkmers[qnkmer-1], param.kmer_size) == '*' ) {
                int offset = rnkmer - me - (qnkmer-qe);
                if ( offset > 2*(int)spur_offset ) {
                    if ( param.verbose ) std::cout << "Stop codon with long offset:" << offset << "\n";
                    type = NOTYPE;
                }
            }
        }
    }

    // ooooooooooxxxxxxxxxx
    // oooooooooo----------
    else if ( qe == qnkmer-1 ) {
        // 1. beginning of reference sequence
        // ----------ooooooooooxxxxxxxxxx
        // xxxxxxxxxxoooooooooo----------
        if ( ms == 0 ) {
            if ( param.verbose ) std::cout << "Type:Latch Left\n";
            if ( alpha::getLastAminoAcid(qkmers[qnkmer-1], param.kmer_size) == '*' ) {
                if ( me == rnkmer-1 && alpha::getLastAminoAcid(rkmers[rnkmer-1], param.kmer_size) == '*' ) type = LATCH;
                else {
                    if (param.verbose) std::cout << "Invalid latch - stop codon\n";
                }
            }
            else if ( sbjct_len < (int)latch_offset  ) {
                if ( param.verbose) std::cout << "BAD LATCH - SHORT:" << sbjct_len << "\n";
            }
            else { if (param.verbose) std::cout << "TRUE LATCH\n";type = LATCH;}
        }
        // 2. mismatch in query sequence
        // xxxxxxxxxxooooooooooxxxxxxxxxx
        // -----xxxxxoooooooooo----------
        else {
            if ( param.verbose ) std::cout << "Type:Spur Left\n";
//             std::cout << "Spur Left\n";
//             if ( (int)qs >= 2*spur_offset )
//                 std::cout << "BAD SPUR - start:" << qs << "\n";
//             else if ( (int)qs >= spur_offset ) {
//                 std::cout << "WEAK SPUR\n"; type = WEAK_SPUR;
//             }
//             else { std::cout << "TRUE SPUR\n"; type = SPUR; }
            //type = SPUR;
            if ( qs >= 2*spur_offset ) {
                if ( param.verbose ) std::cout << "BAD SPUR - start:" << qs << "\n";
            }
            else if ( qs >= spur_offset ) {
                if ( param.verbose ) std::cout << "WEAK SPUR\n";
                type = WEAK_SPUR;
            } else {
                if ( param.verbose ) std::cout << "TRUE SPUR\n";
                type = SPUR;
            }
            if ( type != NOTYPE && alpha::getLastAminoAcid(qkmers[qnkmer-1], param.kmer_size) == '*' ) {
                int offset = ms - qs; 
                if ( offset > 2*(int)spur_offset ) {
                    if ( param.verbose ) std::cout << "Stop codon with long offset:" << offset << "\n";
                    type = NOTYPE;
                }
            }

        }
    }

    // xxxxxxxxxxooooooooooxxxxxxxxxx
    // -------xxxooooooooooxxx-------
    else {
        if ( param.verbose ) std::cout << "Type:Frayed Rope\n";
        if ( qs >= 2*spur_offset || (qnkmer-qe) > 2*spur_offset ) { 
            if ( param.verbose ) std::cout << "BAD FR - start:" << qs << "\n";
        }
        else if ( qs >= spur_offset || (qnkmer-qe) > spur_offset ) {
            if ( param.verbose ) std::cout << "WEAK ROPE\n";
            type = WEAK_ROPE;
        } else {
            if ( param.verbose ) std::cout << "TRUE ROPE\n";
            type = ROPE;
        }
        if ( type != NOTYPE && alpha::getLastAminoAcid(qkmers[qnkmer-1], param.kmer_size) == '*' ) {
            int offset = rnkmer - me - (qnkmer-qe);
            if ( offset > 2*(int)spur_offset ) {
                if ( param.verbose ) std::cout << "Stop codon with long offset:" << offset << "\n";
                type = NOTYPE;
            }
        }
    }

    return type;
}

/* Bubble */
// xxxxxooooooooooxxxxx
// -----oooooooooo-----
bool assem::isBubble( size_t qs,
                      size_t qe,
                      size_t ms,
                      size_t me,
                      size_t qnkmer,
                      size_t rnkmer,
                      KmerId *qkmers,
                      KmerId *rkmers,
                      size_t spur_offset,
                      int &type,
                      Param &param )
{
    if ( qs != 0 || qe != qnkmer-1 ) return false;

    std::string sbjct = biostr::getSequenceString( &rkmers[ms], me-ms+1, param.kmer_size );
    std::string query = biostr::getSequenceString( &qkmers[qs], qe-qs+1, param.kmer_size );
    
    if ( sbjct == query ) { 
        if (param.verbose) std::cout << "Substring\n";
        type = SUBSTRING;
        return true;
    }

    type = BUBBLE;
    if ( alpha::getLastAminoAcid(qkmers[qnkmer-1], param.kmer_size) != '*' ) {
        if ( param.verbose ) std::cout << "Bubble\n";
        return true;
    }
    else {
        int offset = rnkmer - me - 1;
        if ( offset > 2*(int)spur_offset ) {
            if ( param.verbose ) 
                std::cout << "Stop condon with long offset:" << offset << "\tnot a bubble\n";
            type = NOTYPE;
            return false;
        }
        else {
            std::cout << "Stop condon with short offset:" << offset << "\tbubble\n";
        }
    }
    return true;
}

//
bool assem::isSpurRight( size_t qs,
                         size_t qe,
                         size_t ms,
                         size_t me,
                         size_t qnkmer,
                         size_t rnkmer,
                         KmerId *qkmers,
                         KmerId *rkmers,
                         size_t spur_offset,
                         int &type,
                         Param &param )
{
    if ( qs > 0 ) return false;

    int query_len = qnkmer + param.kmer_size - 1;
//     if ( (query_len-qe) > 2*spur_offset ) { 
//         if ( param.verbose ) std::cout << "Bad Spur Right:\tlength:" << query_len << "\tend:" << qe << "\n";;
//         return false;
//     }
    
//     if ( (query_len-qe) > spur_offset ) {
// //         if ( param.verbose ) std::cout << "Weak Spur\n";
// //         type = WEAK_SPUR;
//         return false;
//     } 
    if ( (qnkmer-1-qe) > 2*spur_offset ) { 
        if ( param.verbose ) std::cout << "Bad Spur Right:\tlength:" << query_len << "\tend:" << qe << "\n";;
        return false;
    }
    
    if ( (qnkmer-1-qe) > spur_offset ) {
        if ( param.verbose ) std::cout << "Weak Spur\n";
        type = WEAK_SPUR;
        //return false;
    } 

    else {
        if ( param.verbose ) std::cout << "True Spur Right\n";
        type = SPUR;
    }
    if ( type != NOTYPE && alpha::getLastAminoAcid(qkmers[qnkmer-1], param.kmer_size) == '*' ) {
        int offset = rnkmer - me - (qnkmer-qe);
        if ( offset > 2*(int)spur_offset ) {
            if ( param.verbose ) std::cout << "Stop codon with long offset:" << offset << "\n";
            type = NOTYPE;
            return false;
        }
    }
    return true;
}

bool assem::isSpurLeft(size_t qs,
                       size_t qe,
                       size_t ms,
                       size_t me,
                       size_t qnkmer,
                       size_t rnkmer,
                       KmerId *qkmers,
                       KmerId *rkmers,
                       size_t spur_offset,
                       int &type,
                       Param &param )
{
    if ( qe < qnkmer-1 ) return false;
    
    if ( qs >= 2*spur_offset ) {
        if ( param.verbose ) std::cout << "Bad Spur Left - start:" << qs << "\n";
        return false;
    }
  
    if ( qs >= spur_offset ) {
        if ( param.verbose ) std::cout << "Weak Spur\n";
        type = WEAK_SPUR;
        //return false;
    } else {
        if ( param.verbose ) std::cout << "True Spur Left\n";
        type = SPUR;
    }

    if ( type != NOTYPE && alpha::getLastAminoAcid(qkmers[qnkmer-1], param.kmer_size) == '*' ) {
        int offset = rnkmer -me - 1;
        if ( offset > 2*(int)spur_offset ) {
            if ( param.verbose ) std::cout << "Stop codon with long offset:" << offset << "\n";
            type = NOTYPE;
            return false;
        }
    }
    return true;
}

    
//       qs <= spur_offset
//       me >= rnkmer-1-spur_offset
//       query right end longer than sbjct
//    
// xxxxxxxxxxooooooooooxx
// -------xxxooooooooooxxxxxxxx-------
// xxxxxxxxxxoooooooooo----------
// ----------ooooooooooxxxxxxxxxx
bool assem::isLatchRight(size_t qs,
                         size_t qe,
                         size_t ms,
                         size_t me,
                         size_t qnkmer,
                         size_t rnkmer,
                         KmerId *qkmers,
                         KmerId *rkmers,
                         size_t spur_offset,
                         //int &type,
                         Param &param )
{

    //if ( qs > spur_offset || me < rnkmer-1-spur_offset ) return false;

    int qlen = qnkmer + param.kmer_size + 1;
    int rlen = rnkmer + param.kmer_size + 1;

    if ( qs > qlen*.5 || me < rlen*.5 ) {
        if ( param.verbose ) std::cout << "Invalid range - latch right\n";
        return false;
    }

    int soff = rnkmer - me + 1;
    int qoff = qnkmer - qe + 1;
    
    if ( qoff <= soff ) return false;

    if ( alpha::getLastAminoAcid(rkmers[rnkmer-1], param.kmer_size) == '*' ) {
        //if ( qe == qnkmer-1 && alpha::getLastAminoAcid(qkmers[qnkmer-1], param.kmer_size) == '*' ) {
        if ( alpha::getLastAminoAcid(qkmers[qnkmer-1], param.kmer_size) == '*' ) {
            if (param.verbose) std::cout << "Alternative stop site extension\n";
        }
        else {
            if (param.verbose) std::cout << "Invalid latch - stop codon\n";
            return false;
        }
    } 
    if ( param.verbose ) std::cout << "Latch Right\n";       
    return true;
}

// ----xxxooooooooooxxxxxxxxx
// xxxxxxxooooooooooxx-------
// ----------ooooooooooxxxxxxxxxx
// xxxxxxxxxxoooooooooo----------
bool assem::isLatchLeft(size_t qs,
                      size_t qe,
                      size_t ms,
                      size_t me,
                      size_t qnkmer,
                      size_t rnkmer,
                      KmerId *qkmers,
                      KmerId *rkmers,
                      size_t spur_offset,
                         //int &type,
                      Param &param )
{
    //if ( ms > spur_offset || qe < qnkmer-1-spur_offset ) return false;

    int qlen = qnkmer + param.kmer_size + 1;
    int rlen = rnkmer + param.kmer_size + 1;
    
    if ( ms > rlen*.5 || qe < qlen*.5 ) {
        if ( param.verbose ) std::cout << "Invalid range - latch left\n";
        return false;
    }


    int soff = rnkmer - me + 1;
    int qoff = qnkmer - qe + 1;
    
    if ( soff <= qoff ) return false;
    
    //if ( param.verbose ) std::cout << "Type:Latch Left\n";
    if ( alpha::getLastAminoAcid(qkmers[qnkmer-1], param.kmer_size) == '*' ) {
        //if ( me == rnkmer-1 && alpha::getLastAminoAcid(rkmers[rnkmer-1], param.kmer_size) == '*' ) {
        if ( alpha::getLastAminoAcid(rkmers[rnkmer-1], param.kmer_size) == '*' ) {
            if (param.verbose) std::cout << "Alternative stop site extension\n";
        }
        else {
            if (param.verbose) std::cout << "Invalid latch - stop codon\n";
            return false;
        }
    }
    if ( param.verbose ) std::cout << "Latch Left\n";
    return true;
}

// xxxxxxxxxxooooooooooxxxxxxxxxx
// -------xxxooooooooooxxx-------
bool assem::isFrayedRope(size_t qs,
                      size_t qe,
                      size_t ms,
                      size_t me,
                      size_t qnkmer,
                      size_t rnkmer,
                      KmerId *qkmers,
                      KmerId *rkmers,
                      size_t spur_offset,
                         int &type,
                      Param &param )
{
    if ( qs >= 2*spur_offset || (qnkmer-1-qe) > 2*spur_offset ) { 
        if ( param.verbose ) std::cout << "Bad frayed rope\n";
        return false;
    }
    
    else if ( qs >= spur_offset || (qnkmer-1-qe) > spur_offset ) {
        if ( param.verbose ) std::cout << "Weak Rope\n";
        type = WEAK_ROPE;
        //return false;
    } else {
        if ( param.verbose ) std::cout << "True Rope\n";
        type = ROPE;
    }
    if ( type != NOTYPE && alpha::getLastAminoAcid(qkmers[qnkmer-1], param.kmer_size) == '*' ) {
        int offset = rnkmer - me - (qnkmer-qe);
        if ( offset > 2*(int)spur_offset ) {
            if ( param.verbose ) std::cout << "Stop codon with long offset:" << offset << "\n";
            type = NOTYPE;
            return false;
        }
    }
    return true;
}

// int assem::determinePairEndLinkType( std::vector<int> &pos_vec,
//                                      IntPair &qrange,
//                                      IntPair &rrange,
//                                      KmerId *qkmers,
//                                      KmerId *rkmers,
//                                      size_t qnkmer,
//                                      size_t rnkmer,
//                                      size_t npairs,
//                                      int direction,
//                                      Param &param  )
// {
//     if ( direction == POSITIVE ) {
//         rrange = expandBound(blocks.back(), pos_vec, param);
//         if ( rrange.first == -1 || rrange.second == -1 ) {
//             if ( param.verbose ) std::cout << "Warning:Latch right fail\n";
//             return NOTYPE;
//         }
//         if ( param.verbose ) std::cout << "Latch right with pair-end\n";
//         qrange = IntPair( pos_vec[rrange.first], pos_vec[rrange.second] );
//         adjustRanges(qrange, rrange, qnkmer, rnkmer, param );
//         return LATCH;
//     }
//     else {
//         rrange = expandBound(blocks.front(), pos_vec, param);
//         if ( rrange.first == -1 || rrange.second == -1 ) {
//             if ( param.verbose ) std::cout << "Warning:Latch right fail\n";
//             return NOTYPE;
//         }
//         if ( param.verbose ) std::cout << "Latch right with pair-end\n";
//         qrange = IntPair( pos_vec[rrange.first], pos_vec[rrange.second] );
//         adjustRanges(qrange, rrange, qnkmer, rnkmer, param );
//         return LATCH;
//     }
// }

int assem::getPathAlignType( AlignSummary &summary )
{
//     //xxxxxxxxxxoooooooooo----------
//     //----------ooooooooooxxxxxxxxxx
//     if ( summary.lgap.first == 0 && summary.lgap.second > 0 && summary.egap.first > 0 && summary.egap.second == 0 ) 
//         return LATCH; // extend right

//     //----------ooooooooooxxxxxxxxxx
//     //xxxxxxxxxxoooooooooo----------
//     if ( summary.lgap.first > 0 && summary.lgap.second == 0 && summary.egap.first == 0 && summary.egap.second >0 ) 
//         return LATCH; // extend left
    
//     //xxxxxxxxxxooooooooooxxxxxxxxxx
//     //-----xxxxxooooooooooxxxxx-----
//     if ( summary.lgap.first == 0 && summary.lgap.second >= 0 && summary.egap.first == 0 && summary.egap.second >= 0 ) 
//         return MERGE;

    if ( summary.s2se.first >= summary.s1se.first && summary.s2se.second <= summary.s1se.second ) return MERGE;
    if ( summary.s2se.first  < summary.s1se.first && summary.s2se.second <  summary.s1se.second ) return LATCH_LEFT;
    if ( summary.s2se.first  > summary.s1se.first && summary.s2se.second >  summary.s1se.second ) return LATCH_RIGHT;
    
    return NOTYPE;
}

int assem::getMatchType( IntPair &qrange,
                         IntPair &rrange,
                         KmerId *qkmers,
                         KmerId *rkmers,
                         int    qnkmer,
                         int    rnkmer,
                         Param &param  )
{
    /* reference start and end postions */
    int ms = rrange.first;
    int me = rrange.second;
    if ( ms == -1 || me == -1 ) return NOTYPE;

    /* query start and end positions */
    int qs = qrange.first;
    int qe = qrange.second;

    int dqe = qnkmer-qe-1;
    int dme = rnkmer-me-1;

    int qlen = qnkmer + param.kmer_size + 1;
    int rlen = rnkmer + param.kmer_size + 1;
    
    // xxxxxooooxoooooxxxxx
    // -----ooooxooooo-----
    if ( qs == 0 && qe == qnkmer-1 ) {
        if ( param.verbose ) std::cout << "Bubble\n";
        return BUBBLE;
    }

    if ( qs == 0 ) {
        // xxxxxxxxxxooooooooooxxxxxxxxxx
        // ----------ooooooooooxxxxx-----
        if ( dqe <= dme ) {
            if ( param.verbose ) std::cout << "Spur right\n";
            return SPUR_RIGHT;
        } 

        // xxxxxxxxxxoooooooooo----------
        // ----------ooooooooooxxxxxxxxxx
        if ( param.verbose ) std::cout << "Simple extend right\n";
        return LATCH_RIGHT_EASY;
    }

    if ( qe == qnkmer-1 ) {
        // xxxxxxxxxxooooooooooxxxxxxxxxx
        // -----xxxxxoooooooooo----------
        if ( qs <= ms ) {
            if ( param.verbose ) std::cout << "Spur left\n";
            return SPUR_LEFT;
        }

        // ----------ooooooooooxxxxxxxxxx
        // xxxxxxxxxxoooooooooo----------
        if ( param.verbose ) std::cout << "Simple extend left\n";
        return LATCH_LEFT_EASY;
    }

    // xxxxxxxxxxooooooooooxxxxxxxxxx
    // -------xxxooooooooooxxx-------
    if ( ( qs > 0 && qe < qnkmer-1 ) &&
         ( qs <= ms && dqe <= dme ) ) {
        if ( param.verbose ) std::cout << "Frayed rope\n";
        return ROPE;
    }

    // xxxxxxxxxxxoxoxoooox----------
    // ----------xoxoxooooxxxxxxxxxxx
    if ( qs <= ms && dqe >= dme ) {
        if ( param.verbose ) std::cout << "Complex extend right\n";
        if ( qs > qlen*.3 || me < rlen*.7 ) {
            if ( param.verbose ) std::cout << "Invalid latch range\n";
            return NOTYPE;
        }
        return LATCH_RIGHT_DIFF;
    }

    // ----------xooxoxooxxxxxxxxxxxx
    // xxxxxxxxxxxooxoxooxx----------
    if ( qs >= ms && dqe <= dme ) {
        if ( param.verbose ) std::cout << "Complex extend Left\n";
        if ( ms > rlen*.3 || qe < qlen*.7 ) {
            if ( param.verbose ) std::cout << "Invalid latch range\n";
            return NOTYPE;
        }
        return LATCH_LEFT_DIFF;
    }
    return NOTYPE;
}

int assem::getPathPairType( std::vector<int> &pos_vec,
                            IntPair &qrange,
                            IntPair &rrange,
                            KmerId *qkmers,
                            KmerId *rkmers,
                            size_t qnkmer,
                            size_t rnkmer,
                            size_t latch_offset, 
                            size_t spur_offset,
                            Param &param  )
{
    if ( rrange.first == -1 || rrange.second == -1 ) return NOTYPE;

    /* reference start and end postions */
    size_t ms = rrange.first;
    size_t me = rrange.second;

    /* query start and end positions */
    size_t qs = qrange.first;
    size_t qe = qrange.second;

    int dqe = qnkmer-qe-1;
    int dme = rnkmer-me-1;

    
    // xxxxxooooxoooooxxxxx
    // -----ooooxooooo-----
    if ( qs == 0 && qe == qnkmer-1 ) {
        if ( param.verbose ) std::cout << "Bubble\n";
        return BUBBLE;
    }

    if ( qs == 0 ) {
        // xxxxxxxxxxooooooooooxxxxxxxxxx
        // ----------ooooooooooxxxxx-----
        if ( dqe <= dme ) {
            if ( param.verbose ) std::cout << "Spur right\n";
            return SPUR;
        } 

        // xxxxxxxxxxoooooooooo----------
        // ----------ooooooooooxxxxxxxxxx
        if ( param.verbose ) std::cout << "Simple extend right\n";
        return LATCH;
    }

    if ( qe == qnkmer-1 ) {
        // xxxxxxxxxxooooooooooxxxxxxxxxx
        // -----xxxxxoooooooooo----------
        if ( qs <= me ) {
            if ( param.verbose ) std::cout << "Spur left\n";
            return SPUR;
        }

        // ----------ooooooooooxxxxxxxxxx
        // xxxxxxxxxxoooooooooo----------
        if ( param.verbose ) std::cout << "Simple extend left\n";
        return LATCH;
    }

    // xxxxxxxxxxooooooooooxxxxxxxxxx
    // -------xxxooooooooooxxx-------
    if ( ( qs > 0 && qe < qnkmer-1 ) &&
         ( qs <= ms && dqe <= dme ) ) {
        if ( param.verbose ) std::cout << "Frayed rope\n";
        return ROPE;
    }

    // xxxxxxxxxxxoxoxoooox----------
    // ----------xoxoxooooxxxxxxxxxxx
    if ( qs <= ms && dqe >= dme ) {
        if ( param.verbose ) std::cout << "Complex extend right\n";
        return LATCH;
    }

    // ----------xooxoxooxxxxxxxxxxxx
    // xxxxxxxxxxxooxoxooxx----------
    if ( qs >= ms && dqe <= dme ) {
        if ( param.verbose ) std::cout << "Complex extend Left\n";
        return LATCH;
    }
    return NOTYPE;
}

// int assem::getPathPairType( std::vector<int> &pos_vec,
//                             IntPair &qrange,
//                             IntPair &rrange,
//                             KmerId *qkmers,
//                             KmerId *rkmers,
//                             size_t qnkmer,
//                             size_t rnkmer,
//                             size_t latch_offset, 
//                             size_t spur_offset,
//                             Param &param  )
// {
//     int type = NOTYPE;
    
//     /* reference start and end postions */
//     size_t ms = rrange.first;
//     size_t me = rrange.second;

//     /* query start and end positions */
//     size_t qs = qrange.first;
//     size_t qe = qrange.second;


//     if ( ms == -1 || me == -1 ) return NOTYPE;


//     //int query_len = qnkmer + param.kmer_size - 1;
//     //int sbjct_len = me - ms + param.kmer_size;

// //     std::string sbjct = biostr::getSequenceString( &rkmers[ms], me-ms+1, param.kmer_size );
// //     std::string query = biostr::getSequenceString( &qkmers[qs], qe-qs+1, param.kmer_size );
    
// //     if ( sbjct == query ) { 
// //         if (param.verbose) std::cout << "Substring\n";
// //         return SUBSTRING;
// //     }
//     if ( isBubble( qs, qe, ms, me, qnkmer, rnkmer, qkmers, rkmers, spur_offset, type, param ) ) return type;
//     if ( isSpurRight( qs, qe, ms, me, qnkmer, rnkmer, qkmers, rkmers, spur_offset, type, param ) ) return type;
//     if ( isSpurLeft( qs, qe, ms, me, qnkmer, rnkmer, qkmers, rkmers, spur_offset, type, param ) ) return type;
//     if ( isFrayedRope( qs, qe, ms, me, qnkmer, rnkmer, qkmers, rkmers, spur_offset, type, param  ) ) return type;

//     IntPair nrange;
//     nrange = getSimpleLatchRightRange(pos_vec);
//     if ( nrange.first != -1 && nrange.second != -1 && param.kmer_size + (nrange.second - nrange.first) >= latch_offset ) {
//         if ( param.verbose ) std::cout << "Simple latch right\n";
//         rrange = nrange; 
//         qrange.first = pos_vec[nrange.first];
//         qrange.second = pos_vec[nrange.second];
//         return LATCH;
//     }
//     nrange = getSimpleLatchLeftRange(pos_vec, qnkmer);
//     if ( nrange.first != -1 && nrange.second != -1 && param.kmer_size + (nrange.second - nrange.first) >= latch_offset ) {
//         if ( param.verbose ) std::cout << "Simple latch left\n";
//         rrange = nrange; 
//         qrange.first = pos_vec[nrange.first];
//         qrange.second = pos_vec[nrange.second];
//         return LATCH;
//     }

//     if ( isLatchRight( qs, qe, ms, me, qnkmer, rnkmer, qkmers, rkmers, spur_offset, param  ) ) return LATCH;
//     if ( isLatchLeft( qs, qe, ms, me, qnkmer, rnkmer, qkmers, rkmers, spur_offset, param ) ) return LATCH;


//     return NOTYPE;
// }



// Why doing this?
// A score is calculated from only match region
// But, entire region information is required for later merging.
AlignSummary assem::compareBySubstitution( std::string &query,
                                    std::string &sbjct,
                                    int qnkmer,
                                    int rnkmer,
                                    int qs,
                                    int qe,
                                    int rs,
                                    int re,
									Param &param  )
{
    if ( param.verbose ) std::cout << "Compare by Substitution\n";
    if ( query.size() != sbjct.size() ) {
        if ( param.verbose ) std::cout << "Not same length\n";
        return AlignSummary();
    }
    AlignSummary summary;

    summary.length = query.size();
    if ( param.verbose ) std::cout << "query len:" << query.size() << "\n";

    //-----------------
    // Alignment scores
    //-----------------
    summary.score = scoring::sumScore(query, sbjct, BLOSUM62);    
    int mcount = 0;
    for ( size_t i = 0; i < query.size(); i++ )
        if ( query[i] == sbjct[i] ) mcount++;
    summary.match = mcount;
    summary.mismatch = query.size() - mcount;
    summary.positive = scoring::countPositive(query, sbjct, BLOSUM62);
    summary.posrate = (double)summary.positive/summary.length;

    //------------
    // leading gap
    //------------
    if ( rs > 0 ) summary.lgap = IntPair( 0, rs );
    else if ( qs > 0 ) summary.lgap = IntPair( qs, 0 );

    //-------------
    // trailing gap
    //-------------
    if ( re < rnkmer-1 ) summary.egap = IntPair( 0, rnkmer-1 - re );
    else if ( qe < qnkmer-1 ) summary.egap = IntPair( qnkmer-1 - qe, 0 );

    //---------------------------------
    // Sequence begin and end positions
    //---------------------------------
    summary.s1se = IntPair(rs,re+param.kmer_size-1);
    summary.s2se = IntPair(qs,qe+param.kmer_size-1);

    //----------------------------
    // Alignment outer/inner bound
    //----------------------------
    int ne = re;
    if ( qs > 0 ) ne += qs;
    if ( qe < qnkmer-1 ) ne += (qnkmer-1-qe);
    summary.outer = IntPair(0, ne+param.kmer_size-1);

    if ( qs > 0 ) summary.range = summary.s2se;
    else if ( rs > 0 ) summary.range = summary.s1se;

    
    if ( param.verbose ) summary.print(std::cout);
    return summary;
}

AlignSummary assem::compareByAlignment( std::string query,
                                 std::string sbjct,
                                 int qnkmer,
                                 int rnkmer,
                                 int qs,
                                 int qe,
                                 int rs,
                                 int re,
								 Param &param )
{
    if ( param.verbose ) std::cout << "Compare by Alignment\n";
    GlobalAlignPair paln = GlobalAlignPair(sbjct,query);
    AlignSummary summary = paln.getSummary();

    if ( param.verbose ) {
        summary.print(std::cout);
//         std::cout << paln.getAlignment();
//         std::cout << "\t#ins:" << summary.ins.size() << "\t";
//         std::cout << "\t#del:" << summary.del.size() << "\t";
//         std::cout << "\t#mat:" << summary.match << "\t";
//         std::cout << "\t#pos:" << summary.positive << "\t";
//         std::cout << "\t#mms:" << summary.mismatch << "\n";
//         std::cout << "\tseq1:" << summary.s1se.first << "\t" << summary.s1se.second << "\n";
//         std::cout << "\tseq2:" << summary.s2se.first << "\t" << summary.s2se.second << "\n";
//         std::cout << "\tlgap:" << summary.lgap.first << "\t" << summary.lgap.second << "\n";
//         std::cout << "\tegap:" << summary.egap.first << "\t" << summary.egap.second << "\n";
//         std::cout << "\t#out:" << summary.outer.first << "\t" << summary.outer.second << "\n";
//         std::cout << "\t#in:" << summary.range.first << "\t" << summary.range.second << "\n";
//         std::cout << "\t#len:" << summary.range.second-summary.range.first+1 << "\n";
//         std::cout << "\t%pos:" << summary.posrate << "\n";
    }


    //--------------------------------
    // Query start position adjustment
    //--------------------------------
    //xxxxxxxxxxxxoooooooooo
    //------------oooooooooo
    if ( summary.lgap.second > 0 ) {
        summary.s2se.first  -= summary.lgap.second;
        summary.s2se.second -= summary.lgap.second;
    }


    
    //------------------------------------------------------
    // Leading gap adjustment
    // Alignment start of sbjct/query is not the first base.
    //------------------------------------------------------
    if ( rs > 0 ) summary.lgap.second += rs;
    else if ( qs > 0 ) summary.lgap.first += qs;


    //-----------------------------------------------------
    // Trailing gap adjustment
    // Alignment end of sbjct/query is not end of sequence.
    //-----------------------------------------------------
    if ( re < rnkmer-1 ) summary.egap.second += ( rnkmer-1 - re );
    else if ( qe < qnkmer-1 ) summary.egap.first += ( qnkmer-1 - qe );



    //-----------------------------------------------------
    // Sequence begin and end positions of sbjct and query.
    //-----------------------------------------------------
    if ( qs > 0 ) {
        summary.s2se.first  += qs;
        summary.s2se.second += qs;
    } else if ( rs > 0 ) {
        summary.s1se.first  += rs;
        summary.s1se.second += rs;
    }


 
//     int ne = re;
//     if ( qs > 0 ) ne += qs;
//     if ( qe < qnkmer-1 ) ne += (qnkmer-1-qe);
//     summary.outer = IntPair(0, ne);

    //---------------------------------
    // Alignment inner bound adjustment
    //---------------------------------
    if ( qs > 0 ) summary.range = summary.s2se;
    else if ( rs > 0 ) {
        summary.range.first  += rs;
        summary.range.second += rs;
    }

    //-----------------------------
    // Alignment outer bound adjust
    //-----------------------------
    if ( summary.range.second > summary.outer.second ) 
        summary.outer.second = summary.range.second;

    for ( AlignPosList::iterator it = summary.ilist.begin(); it != summary.ilist.end(); ++it ) {
        if ( rs > 0 ) it->ref_pos += rs;
        else if ( qs > 0 ) it->seq_pos += qs;
    }
    for ( AlignPosList::iterator it = summary.dlist.begin(); it != summary.dlist.end(); ++it ) {
        if ( rs > 0 ) it->ref_pos += rs;
        else if ( qs > 0 ) it->seq_pos += qs;
    }

    if ( param.verbose ) {
        std::cout << "Position adjustment\n";
        std::cout << "\tseq1:" << summary.s1se.first << "\t" << summary.s1se.second << "\n";
        std::cout << "\tseq2:" << summary.s2se.first << "\t" << summary.s2se.second << "\n";
        std::cout << "\tlgap:" << summary.lgap.first << "\t" << summary.lgap.second << "\n";
        std::cout << "\tegap:" << summary.egap.first << "\t" << summary.egap.second << "\n";
        std::cout << "\t#out:" << summary.outer.first << "\t" << summary.outer.second << "\n";
        std::cout << "\t#in:" << summary.range.first << "\t" << summary.range.second << "\n";
        std::cout << "\t#len:" << summary.range.second-summary.range.first+1 << "\n";
    }
        
    return summary;
}


/*
  Buggy???
  1  2  3 100
  9 10 11 12
 */
void assem::adjustRanges( IntPair &qrange,
                   IntPair &rrange,
                   int qnkmer,
                   int rnkmer,
				   Param &param )
{
    int qs = qrange.first;
    int qe = qrange.second;
    int ms = rrange.first;
    int me = rrange.second;
    if ( param.verbose ) {
        std::cout << "qs:" << qs << "\tqe:" << qe << "\n";
        std::cout << "ms:" << ms << "\tme:" << me << "\n";
    }
    if ( qs > 0 ) {
        if ( qs > ms ) {
            qs = qs-ms;
            ms = 0;
        }
        else {
            ms = ms-qs;                
            qs = 0;
        }
    }
    if ( qe < (int)qnkmer-1 ) {
        int qdiff = qnkmer-1 - qe;
        int rdiff = rnkmer-1 - me;
        if ( qdiff > rdiff ) {
            qe += rdiff;
            me = rnkmer-1;
        } else {
            qe = qnkmer-1;
            me += qdiff;
        }
    }
    // update range
    qrange = IntPair(qs, qe);
    rrange = IntPair(ms, me);
    
    if ( param.verbose ) {
        std::cout << "qs:" << qs << "\tqe:" << qe << "\n";
        std::cout << "ms:" << ms << "\tme:" << me << "\n";    
    }
}

void assem::adjustLatchRange( IntPair &qrange,
                              IntPair &rrange,
                              int qnkmer,
                              int rnkmer,
                              Param &param )
{
    int qs = qrange.first;
    int qe = qrange.second;
    int ms = rrange.first;
    int me = rrange.second;
    if ( param.verbose ) {
        std::cout << "qs:" << qs << "\tqe:" << qe << "\n";
        std::cout << "ms:" << ms << "\tme:" << me << "\n";
    }

    // Latch right
    if ( ms > qs ) {
        if ( qs > 0 ) {
            if ( param.verbose ) 
                std::cout << "Latch range adjust - start\n";
            ms -= qs;
            qs = 0;
        }
        if ( me < rnkmer-1 ) {
            if ( param.verbose ) 
                std::cout << "Latch range adjust - end\n";
            
            int rdiff = rnkmer-1-me;
            int qdiff = qnkmer-1-qe;            
            if ( qdiff > rdiff ) {
                qe += rdiff;
                me = rnkmer-1;
            } else {
                qe = qnkmer-1;
                me += qdiff;
            }
        }
    }

    // Latch left
    else {
        if ( ms > 0 ) {
            if ( param.verbose ) 
                std::cout << "Latch range adjust - start\n";
            qs -= ms;
            ms = 0;
        }
        if ( qe < qnkmer-1 ) {
            if ( param.verbose ) 
                std::cout << "Latch range adjust - end\n";

            me += (qnkmer-1-qe);
            qe = qnkmer-1;
        }
    }

    qrange = IntPair(qs, qe);
    rrange = IntPair(ms, me);
    
    if ( param.verbose ) {
        std::cout << "qs:" << qs << "\tqe:" << qe << "\n";
        std::cout << "ms:" << ms << "\tme:" << me << "\n";    
    }
}


void assem::Interchange(PathId &qid, 
                 PathId &mid, 
                 KmerId **qkmers, 
                 KmerId **rkmers, 
                 size_t &qnkmer, 
                 size_t &rnkmer, 
                 IntPair &qrange, 
                 IntPair &rrange, 
				 Param &param )
{
    if ( param.verbose ) {
        std::cout << "Before\n";
        std::cout << "qid:" << qid << "\tmid:" << mid << "\n";
        std::cout << "qnkmer:" << qnkmer << "\trnkmer:" << rnkmer << "\n";
        std::cout << "Query:" << biostr::getSequenceString(*qkmers, qnkmer, param.kmer_size) << "\n";
        std::cout << "Sbjct:" << biostr::getSequenceString(*rkmers, rnkmer, param.kmer_size) << "\n";
        std::cout << "Qrange:" << qrange.first << "\t" << qrange.second << "\n";
        std::cout << "Rrange:" << rrange.first << "\t" << rrange.second << "\n";
    }
    util::swap<size_t>(qnkmer, rnkmer);
    util::swapPtr<KmerId*>(qkmers, rkmers);
    util::swap<IntPair>(qrange, rrange);
    util::swap<PathId>(qid, mid);
    if ( param.verbose ) {
        std::cout << "After\n";
        std::cout << "qid:" << qid << "\tmid:" << mid << "\n";
        std::cout << "qnkmer:" << qnkmer << "\trnkmer:" << rnkmer << "\n";
        std::cout << "Query:" << biostr::getSequenceString(*qkmers, qnkmer, param.kmer_size) << "\n";
        std::cout << "Sbjct:" << biostr::getSequenceString(*rkmers, rnkmer, param.kmer_size) << "\n";
        std::cout << "Qrange:" << qrange.first << "\t" << qrange.second << "\n";
        std::cout << "Rrange:" << rrange.first << "\t" << rrange.second << "\n";
    }
}

// bool assem::pickMergePath( PosPathPairList &pos_paths,
//                     PathId &qid,
//                     PathId &mid,
//                     AlignSummary &summary,
//                     PathToAlnMap &path2aln_map,
//                     KmerToPathMap &pathid_map,
//                     ReadId *pairs, 
//                     size_t nreads,
// 					Param &param  )
// {
//     size_t qnkmer = path2aln_map[qid]->getKmerCount();
//     KmerId *qkmers = path2aln_map[qid]->getKmers();
//     size_t qnread = path2aln_map[qid]->getReadCount();
//     ReadId *qreads = path2aln_map[qid]->getReads();
//     int mk_size = filter::minSameKmerCount( qnkmer+param.kmer_size-1, param.kmer_size, 1-param.filter_score);
//     if ( mk_size < 1 ) mk_size = 1;
//     if ( param.verbose ) {
//         //std::cout << "filter score:" << param.filter_score << "\n";
//         std::cout << "Query length:" << qnkmer+param.kmer_size-1 << "\n";
//         std::cout << "Min same kmers:" << mk_size << "\n";
//     }

//     ReadIdArray qrids = ReadIdArray(qreads, qreads+qnread);
//     std::sort(qrids.begin(), qrids.end());
//     ReadIdArray qpair;
//     if ( param.pair_flag ) {
//         ReadIdList lpair;
//         for ( size_t i = 0; i < qnread; i++ ) 
//             if ( pairs[qreads[i]] != NOT_PAIR ) lpair.push_back( pairs[qreads[i]] );
//         qpair = ReadIdArray( lpair.begin(), lpair.end() );
//         std::sort(qpair.begin(), qpair.end());
//     }
//     for ( PosPathPairList::iterator it = pos_paths.begin(); it != pos_paths.end(); ++it ) {
//         mid = it->second;
//         size_t rnkmer = path2aln_map[mid]->getKmerCount();
//         KmerId *rkmers = path2aln_map[mid]->getKmers();
//         size_t rnread = path2aln_map[mid]->getReadCount();
//         ReadId *rreads = path2aln_map[mid]->getReads();
                
//         if ( param.verbose ) {
//             std::cout << "Sbjct:\t" << mid << "\t" << biostr::getSequenceString(rkmers, rnkmer, param.kmer_size) << "\n";
//         }

//         std::vector<int>::iterator jt;
//         //if ( param.verbose ) std::cout << "pos_vec size:" << (*it).first.size() << "\n";
//         IntPair rrange = getMatchRange(it->first, param);
//         if ( rrange.first == -1 || rrange.second == -1 ) {
//             if ( param.verbose ) std::cout << "No match range\n"; 
//             continue;
//         }
//         IntPair qrange; 
//         qrange.first = (*it).first[rrange.first];
//         qrange.second = (*it).first[rrange.second];
 
//         if ( qrange.second < qrange.first ) {
//             if ( param.verbose ) std::cout << "Query off ordered\n"; 
//             continue;
//         }

//         int mcount = countMatchKmers(it->first, rrange.first, rrange.second);
//         if ( param.verbose ) {
//             std::cout << "Reference length:" << rnkmer+param.kmer_size-1 << "\n";
//             std::cout << "sbjct kmer:" << mcount << "\n";
//         }
        
//         if ( param.pair_flag ) {
//             if ( param.verbose ) {
//                 std::cout << "Query pair read:" << qpair.size() << "\n";
//                 std::cout << "Query read nums:" << qnread << "\n";
//                 std::cout << "Sbjct read nums:" << rnread << "\n";
//             }

//             //std::sort( rreads, rreads+rnread );
//             ReadIdArray rrids = ReadIdArray( rreads, rreads+rnread );
//             std::sort( rrids.begin(), rrids.end() );

//             //ReadIdArray same = Set::Intersect<ReadId>(qrids, rrids, true);
//             //if ( param.verbose ) std::cout << "Query == Sbjct?:" << same.size() << "\n";

//             //double ts = mytime();
//             ReadIdArray comm = Set::Intersect<ReadId>(qpair, rrids, true);
//             if ( param.verbose ) {
//                 std::cout << "Paired reads:" << comm.size() << "\n";
//                 //std::cout << mytime() - ts << " sec\n";
//             }
//         }

//         bool indel_flag = indeled(it->first, rrange);        

//         bool swapped = false;
//         if ( qnkmer > rnkmer ) {
//             if ( param.verbose ) std::cout << "Now swapping\n";
//             Interchange(qid, mid, &qkmers, &rkmers, qnkmer, rnkmer, qrange, rrange, param);
//             if ( param.verbose ) {
//                 std::cout << "Query:" << biostr::getSequenceString(qkmers, qnkmer, param.kmer_size) << "\n";
//                 std::cout << "Sbjct:" << biostr::getSequenceString(rkmers, rnkmer, param.kmer_size) << "\n";
//             }
//             swapped = true;
//         }

//         //int type = getLinkType( qrange, rrange, qkmers, rkmers, qnkmer, rnkmer, param.latch_length, param.path_spur, param );
//         int type = getPathPairType( it->first, qrange, rrange, qkmers, rkmers, qnkmer, rnkmer, param.latch_length, param.path_spur, param );
//         //if ( param.verbose ) std::cout << "Type:" << type << "\n";


//         if ( type == NOTYPE ) { if ( param.verbose ) std::cout << "NOTYPE\n"; 
//             if ( swapped ) Interchange(qid, mid, &qkmers, &rkmers, qnkmer, rnkmer, qrange, rrange, param);
//             continue; 
//         }

//         //---------------------
//         // kmer based filtering
//         //---------------------
//         if ( type != LATCH && mcount < mk_size ) {
//             if ( param.verbose ) std::cout << "Filtered out\n"; 
//             if ( swapped ) Interchange(qid, mid, &qkmers, &rkmers, qnkmer, rnkmer, qrange, rrange, param);
//             continue;
//         }


// //         if ( type == SPUR || type == WEAK_SPUR || type == ROPE || type == WEAK_ROPE ) 
// //             adjustRanges(qrange, rrange, qnkmer, rnkmer, param);
// //         else if ( type == LATCH )
// //             adjustLatchRange(qrange, rrange, qnkmer, rnkmer, param );
//         adjustRanges(qrange, rrange, qnkmer, rnkmer, param);

//         std::string query = biostr::getSequenceString(&qkmers[qrange.first], qrange.second-qrange.first+1, param.kmer_size);
//         std::string sbjct = biostr::getSequenceString(&rkmers[rrange.first], rrange.second-rrange.first+1, param.kmer_size);
        

//         if ( (int)query.length() < param.latch_length ) {
//             if ( param.verbose ) std::cout << "Short latch overlap\n";
//             continue;
//         }

//         //AlignSummary summary = compareSequencePair(qkmers, rkmers, qrange, rrange, qnkmer, rnkmer, type, indel_flag);
//         //summary = compareSequencePair(qkmers, rkmers, qrange, rrange, qnkmer, rnkmer, type, indel_flag);
        
//         if ( !indel_flag )
//             summary = compareBySubstitution(query, sbjct, qnkmer, rnkmer, qrange.first, qrange.second, rrange.first, rrange.second, param);
//         else
//             summary = compareByAlignment(query, sbjct, qnkmer, rnkmer, qrange.first, qrange.second, rrange.first, rrange.second, param);

//         double score = summary.posrate;
//         if ( param.verbose ) std::cout << "Score:"<< score << "\n";
//         if ( type == LATCH ) {
//             if ( score < param.latch_score ) { if ( param.verbose ) std::cout << "Weak score\n"; 
//                 if ( swapped ) Interchange(qid, mid, &qkmers, &rkmers, qnkmer, rnkmer, qrange, rrange, param);
//                 continue; 
//             }
//         } else {
//             if ( score < param.merge_score ) { if ( param.verbose ) std::cout << "Weak score\n"; 
//                 if ( swapped ) Interchange(qid, mid, &qkmers, &rkmers, qnkmer, rnkmer, qrange, rrange,param);
//                 continue; 
//             }
//         }
//         //return std::pair<PathId, AlignSummary>(it->second, summary);
//         return 1;
//     }
//     //return std::pair<PathId, AlignSummary>( qid, AlignSummary() );
//     return 0;
// }


void assem::printGlobalAlignment( GlobalAlignPair &aln )
{
    aln.printAlignment(std::cout);
}

void assem::printLocalAlignment( LocalAlignPair &aln )
{
    aln.printAlignment(std::cout);
}

void assem::printPaths( PathToAlnMap &path2aln_map, Param &param  )
{
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {           std::cout << "Consensus:" << it->second->getConsensus() << "\n";
    }
}

void assem::dropMergedPaths( PathIdSet &merged_paths, 
                      PathToAlnMap &path2aln_map)
{
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ) {        
        if ( merged_paths.find(it->first) == merged_paths.end() ) {
            ++it;
        } else {
            delete it->second;
            path2aln_map.erase(it++);
        }
    }
}

// std::string assem::getSequence( SpaPath *spath, Param &param  )
// {
//     KmerId *kmers = spath->getKmers();
//     size_t nkmer = spath->getKmerCount();
//     std::string seq = biostr::getSequenceString(kmers, nkmer, param.kmer_size);
//     return seq;
// }

void assem::mergePairedPaths( PathToAlnMap &path2aln_map, 
                              InvertedIndex &iindex,
                              //char **seqs,
                              PathId *used_reads,
                              BitString *bstrs,
                              char *strands,
                              ReadId *pairs,
                              int nreads,
                              Param &param )
{
    if ( ! param.pair_flag || ! param.correction_flag ) return;

    double t0 = mytime();

    double new_merge = param.paired_score;
    double old_merge = param.merge_score;
    param.merge_score = new_merge;
    if ( param.verbose ) std::cout << "Temporary change of latch score for paired end path merging:" << old_merge << "\t" << new_merge << "\n";
    clusterPairedPaths(path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param);

    //param.merge_score = old_score;

    double new_latch = param.paired_score;
    double old_latch = param.latch_score;
    
    param.latch_score = new_latch;
    if ( param.verbose ) std::cout << "Temporary change of latch score for paired end path latching:" << old_latch << "\t" << new_latch << "\n";
    connectOverlappingPairedPaths( path2aln_map, bstrs, strands, pairs, used_reads, nreads, iindex, param );
    //connectLongOverlappingPaths(path2aln_map, bstrs, strands, pairs, used_reads, nreads, iindex, param );
    param.merge_score = old_merge;
    param.latch_score = old_latch;

    clusterPaths(path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param);

    if ( param.verbose ) std::cout << "\nMerging paths by paired paths:" << mytime() - t0 << " sec\n";
}
    
//     double t0 = mytime();
    
//     std::cout << "\nSTAGE 6:\n Correcting read placement ...\n";

//     int npair = 0;
//     int ngood = 0;
//     int npoor = 0;
//     int nweak = 0;
//     PathIdSet merged_paths; // Keep track of merged paths
//     Map<PathIdPair, bool>::Type scanned_pairs;
    
//     PathLengthMap plen_map = getPathLengths( path2aln_map );
//     for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
//         PathId sbjct_pid = it->second;
//         std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();

//         if ( merged_paths.find(sbjct_pid) != merged_paths.end() ) continue;
//         if ( param.verbose) 
//             std::cout << "\nSbjct:" << sbjct_pid << "\t" << sbjct << "\n";
        
//         PathIdArray pids = findPairedPaths( sbjct_pid, used_reads, pairs, path2aln_map[sbjct_pid]->getReads(), path2aln_map[sbjct_pid]->getReadCount(), param );
//         if ( pids.size() == 0 ) break;

//         for ( PathIdArray::iterator pt = pids.begin(); pt != pids.end(); ++pt ) {
//             if ( merged_paths.find( *pt ) != merged_paths.end() ) continue;
//             PathId query_pid = *pt;
//             std::string query = path2aln_map[query_pid]->getConsensusString();
//             if ( query.size() > sbjct_size() ) continue; // already scanned for sure???

//             std::cout << "\nQuery:" << query_pid << "\t" << query << "\n";
//             ReadIdArray qreads, sreads;
//             setPairedReads( qreads, sreads, query_pid, sbjct_pid, used_reads, pairs, path2aln_map );
//             if ( param.verbose ) std::cout << "Pair reads:" << qreads.size() << "\n";

//             ReadIdList pos_sreads, neg_sreads, pos_qreads, neg_qreads;
//             for ( size_t i = 0; i < sreads.size(); i++ ) {
//                 char sstr = strands[sreads[i]];
//                 char qstr = strands[pairs[sreads[i]]];
                
//                 if ( sstr == '+' && qstr == '-' ) {
//                     pos_sreads.push_back(sreads[i]);
//                     neg_preads.push_back(pairs[sreads[i]]);
//                 }
//                 else if ( sstr == '-' && qstr == '+' ) {
//                     neg_sreads.push_back(sreads[i]);
//                     pos_preads.push_back(pairs[sreads[i]]);
//                 }
//             }
            
//             if ( param.verbose ) std::cout << "sr+:" << pos_sreads.size() << "\t" << "sr-:" << neg_sreads.size() << "\t" << "qr+:" << pos_qreads.size() << "\t" << "qr-:" << neg_sreads.size() << "\n";

//             if ( ( pos_sreads.size() >= param.latch_support && neg_qreads.size() >= param.latch_support ) ||
//                  ( neg_sreads.size() >= param.latch_support && pos_qreads.size() >= param.latch_support ) ) {
                
//             }

//             if ( (int)sreads.size() < param.latch_support ) {
//                 if ( param.verbose ) std::cout << "Weak support\n";
//                 continue;
//             }

//             npair++;

//             int type;
//             AlignSummary summary;
//             bool success = alignMiddle( summary, sbjct, query, type, param );
            
//             if ( success ) ngood++;
//             else if ( summary.posrate > 0.5 ) nweak++;
//             else npoor++;
//         }
//     }    

//     std::cout << "#Pairs:" << npair << "\t#Good:" << ngood << "\tWeak:" << nweak << "\t#Poor:" << npoor << "\n";
//     std::cout << "Checking paired paths:" << mytime()-t0 << " sec\n";
// }

void assem::combinePaths( PathToAlnMap &path2aln_map, 
                   InvertedIndex &iindex,
                   //char **seqs,
                   PathId *used_reads,
                   BitString *bstrs,
                   char *strands,
                   ReadId *pairs,
                   int nreads,
				   Param &param )
{
    if ( ! param.merge_flag ) return;

    double t0 = mytime();
    
    //std::cout << "\nSTAGE 2:\nMerging and extending paths ...\n";
    std::cout << "\nClustering paths ...\n";

    PathIdSet paths;
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it)
        paths.insert(it->first);

    //--------------------------------------------------------------
    // Try to merge or latch similar pairs from all discovered paths
    //--------------------------------------------------------------
    PathIdSet success_paths;
    //linkPaths(paths, path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param, PATHMERGE);

    clusterPaths(path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param);
    //success_paths = connectPaths(path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param);
    //recruitPaths(success_paths, path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param);
    //mergePaths(paths, path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param, PATHMERGE, false, false);

    //mergePaths(paths, path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param, PATHMERGE, true, true);

//     paths.clear();
//     for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it)
//         paths.insert(it->first);
//     std::cout << "\nLatching paths ...\n";
//     success_paths = mergePaths(paths, path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param, PATHLATCH, false);

//     std::cout << "\nMerge/Latch paths ...\n";
//     mergePaths(success_paths, path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param, MERGELATCH, true);


    //success_paths = latchPaths( paths, path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param, PATHLATCH );

//     // Path extendsion
//     paths.clear();
//     for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it)
//         paths.insert(it->first);
//     success_paths = linkPaths(paths, path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param, PATHLATCH);

//     // Path merging/extension on success paths
//     success_paths = linkPaths(success_paths, path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param, MERGELATCH);
    
    
//     std::cout << "Path count:" << path2aln_map.size() << "\n";
    std::cout << "Path clustering:" << mytime()-t0 << " sec\n";

    //if ( param.verbose ) printPaths(path2aln_map, param);
}

std::multimap<int, PathId> assem::getSimilarPaths( PathId &query_pid, 
                                                   PathToAlnMap &path2aln_map, 
                                                   KmerToPathMap &pathid_map,
                                                   PathIdSet &merged_paths,
                                                   ReadId *pairs,
                                                   int nreads,
                                                   Param &param )
{
    KmerId *kmers = path2aln_map[query_pid]->getKmers();
    size_t  nkmer = path2aln_map[query_pid]->getKmerCount();    
    PathIdSet self; 
    self.insert(query_pid);

    if ( param.verbose ) 
        std::cout << "\nQuery:" << query_pid << "\t" << biostr::getSequenceString(kmers, nkmer, param.kmer_size) << "\n";
    
    //int mink = filter::minSameKmerCount( param.latch_length, param.kmer_size, 1-param.filter_score );
    //if ( mink < 1 ) mink = 1;
    //int mink = 1;

    //std::map<PathId, int> count_map = getKmerCountMap( self, kmers, nkmer, pathid_map, merged_paths );
    std::map<PathId, int> count_map = getKmerCountMap( self, kmers, nkmer, pathid_map, merged_paths );
    std::multimap<int, PathId> icount_map = util::invert_map<PathId, int>(count_map);

    return icount_map;
}

// PosPathPairList assem::getMergiblePaths( PathId &query_pid, 
//                                          PathToAlnMap &path2aln_map, 
//                                          KmerToPathMap &pathid_map,
//                                          PathIdSet &merged_paths,
//                                          ReadId *pairs,
//                                          int nreads,
//                                          Param &param )
// {
//     KmerId *kmers = path2aln_map[query_pid]->getKmers();
//     size_t  nkmer = path2aln_map[query_pid]->getKmerCount();    
//     PathIdSet self; 
//     self.insert(query_pid);
    
//     if ( param.verbose ) 
//         std::cout << "\nQuery:" << query_pid << "\t" << biostr::getSequenceString(kmers, nkmer, param.kmer_size) << "\n";
    
//     int mink = filter::minSameKmerCount( param.latch_length, param.kmer_size, 1-param.filter_score );
//     if ( mink < 1 ) mink = 1;

//     return searchSimilarPaths(self, kmers, nkmer, path2aln_map, pathid_map, merged_paths, mink, param);
// }


void assem::getMatchRanges( std::vector<IntPair> &rranges, 
                            std::vector<IntPair> &qranges, 
                            std::vector<int> &pos_vec,
                            PathId query_pid, 
                            PathId sbjct_pid, 
                            ReadId *pairs, 
                            int nreads, 
                            Param &param )
{
    IntPair mrange = getMatchRangeMax(pos_vec, param);
    IntPair lrange = getMatchRangeLeft(pos_vec, param);
    IntPair rrange = getMatchRangeRight(pos_vec, param);

    if ( mrange.first == -1 || mrange.second == -1 ) return;
    rranges.push_back(mrange); 
    qranges.push_back( IntPair( pos_vec[mrange.first], pos_vec[mrange.second] ) );

    if ( lrange != mrange ) {
        rranges.push_back(lrange);
        qranges.push_back( IntPair( pos_vec[lrange.first], pos_vec[lrange.second] ) );
    }

    if ( rrange != mrange ) {
        rranges.push_back(rrange);
        qranges.push_back( IntPair( pos_vec[rrange.first], pos_vec[rrange.second] ) );
    }
}

void assem::getMatchTypes(intVec &types, 
                   std::vector<IntPair> rranges, 
                   std::vector<IntPair> qranges, 
                   std::vector<int> &pos_vec,
                   PathId query_pid, 
                   PathId sbjct_pid, 
                   PathToAlnMap &path2aln_map,
                   Param &param )
{
    size_t qnkmer  = path2aln_map[query_pid]->getKmerCount();
    KmerId *qkmers = path2aln_map[query_pid]->getKmers();
    size_t rnkmer  = path2aln_map[sbjct_pid]->getKmerCount();
    KmerId *rkmers = path2aln_map[sbjct_pid]->getKmers();
    
    if ( param.verbose ) std::cout << "# blocks:" << rranges.size() << "\n";
    for ( size_t i = 0; i < rranges.size(); i++ ) {
        int type = getMatchType( qranges[i], rranges[i], qkmers, rkmers, qnkmer, rnkmer, param ); 
        types.push_back(type);
        if ( param.verbose ) {
            std::cout << "Block:" << i << "\ttype:" << type << "\tqrange:" << qranges[i].first << "," << qranges[i].second << "\trrange:" << rranges[i].first << "," << rranges[i].second << "\n";
        }
    }
}

void assem::determineMatchRangesAndTypes( intVec &types, 
                                   std::vector<IntPair> &rranges, 
                                   std::vector<IntPair> &qranges, 
                                   std::vector<int> &pos_vec,
                                   PathId query_pid, 
                                   PathId sbjct_pid, 
                                   ReadId *pairs, 
                                   int nreads, 
                                   PathToAlnMap &path2aln_map,
                                   Param &param )
{
    getMatchRanges(rranges, qranges, pos_vec, query_pid, sbjct_pid, pairs, nreads, param);
    getMatchTypes(types, rranges, qranges, pos_vec, query_pid, sbjct_pid, path2aln_map, param);
}
 
ReadIdArray assem::getCommonReads(PathId query_pid, 
                        PathId  sbjct_pid,
                        PathToAlnMap &path2aln_map, 
                        ReadId *pairs,
                           int nreads, 
                           Param &param )
{
    size_t qnread = path2aln_map[query_pid]->getReadCount();
    ReadId *qreads = path2aln_map[query_pid]->getReads();

    ReadIdList lpair;
    for ( size_t i = 0; i < qnread; i++ ) 
        if ( pairs[qreads[i]] != NOT_PAIR ) lpair.push_back( pairs[qreads[i]] );

    ReadIdArray qpair = ReadIdArray( lpair.begin(), lpair.end() );
    std::sort(qpair.begin(), qpair.end());

    size_t rnread = path2aln_map[sbjct_pid]->getReadCount();
    ReadId *rreads = path2aln_map[sbjct_pid]->getReads();
    ReadIdArray rrids = ReadIdArray( rreads, rreads+rnread );
    std::sort( rrids.begin(), rrids.end() );

    ReadIdArray comm = Set::Intersect<ReadId>(qpair, rrids, true);
    if ( param.verbose ) std::cout << "# pair reads:" << comm.size() << "\n";

    return comm;
}

// int getPairedReadCount( PathId query_pid, 
//                         PathId  sbjct_pid,
//                         PathToAlnMap &path2aln_map, 
//                         ReadId *pairs,
//                         int nreads )
// {
//     ReadIdArray reads = getCommonReads( query_pid, sbjct_pid, path2aln_map, pairs, nreads, param );
//     return reads.size();
// }

// bool __tryPairEndPathMerging( PathId query_pid, 
//                               PathId sbjct_pid, 
//                               intVec &types, 
//                               std::vector<IntPair> &ranges, 
//                               ReadIdArray &sbjpreads, 
//                               int direction, 
//                               PathToAlnMap &path2aln_map, 
//                               Param &param)
// {
//     // best score merging ???
//     // try latch first
//     // then, other merging
    
// }

    

// bool __tryOneEndPathMerging( PathId query_pid, 
//                              PathId sbjct_pid, 
//                              intVec &types, 
//                              std::vector<IntPair> &ranges, 
//                              PathToAlnMap &path2aln_map, 
//                              Param &param)
// {
//     // best score merging
// }

// bool assem::enoughKmerShare( int qnkmer, 
//                       IntPair rrange,
//                       std::vector<int> &pos_vec,
//                       int type, 
//                       Param &param )
// {
    
//     int min_nkmer = filter::minSameKmerCount( qnkmer+param.kmer_size-1, param.kmer_size, 1-param.filter_score);
//     if ( min_nkmer < 1 ) min_nkmer = 1;
//     if ( param.verbose ) {
//         std::cout << "Query length:" << qnkmer+param.kmer_size-1 << "\n";
//         std::cout << "Min same kmers (except latching):" << min_nkmer << "\n";
//     }
    
//     int mcount = countMatchKmers(pos_vec, rrange.first, rrange.second);
//     if ( type != LATCH && mcount < min_nkmer ) {
//         if ( param.verbose ) std::cout << "Not enough kmer sharing:" << mcount << "\n";
//         return false;
//     }
//     return true;
// }

// AlignSummary assem::alignPathPair( PathId query_pid, 
//                                    PathId sbjct_pid, 
//                                    int type,
//                                    IntPair &qrange, 
//                                    IntPair &rrange,
//                                    std::vector<int> &pos_vec,
//                                    PathToAlnMap &path2aln_map,
//                                    Param &param )
// {
//     AlignSummary summary;

//     KmerId *rkmers = path2aln_map[sbjct_pid]->getKmers();
//     size_t rnkmers = path2aln_map[sbjct_pid]->getKmerCount();
//     KmerId *qkmers = path2aln_map[query_pid]->getKmers();
//     size_t qnkmers = path2aln_map[query_pid]->getKmerCount();

//     if ( !enoughKmerShare( qnkmers, rrange, pos_vec, type, param ) ) return summary;

//     bool indel_flag;
//     if ( indeled(pos_vec, rrange) ) indel_flag = true;

//     adjustRanges(qrange, rrange, qnkmers, rnkmers, param);   

//     std::string query = biostr::getSequenceString(&qkmers[qrange.first], qrange.second-qrange.first+1, param.kmer_size);
//     std::string sbjct = biostr::getSequenceString(&rkmers[rrange.first], rrange.second-rrange.first+1, param.kmer_size);
//     if ( param.verbose ) {
//         std::cout << "Query:" << query << "\n";
//         std::cout << "Sbjct:" << sbjct << "\n";
//     }

//     if ( ! indel_flag ) 
//         summary = compareBySubstitution(query, sbjct, qnkmers, rnkmers, qrange.first, qrange.second, rrange.first, rrange.second, param);
//     else 
//         summary = compareByAlignment(query, sbjct, qnkmers, rnkmers, qrange.first, qrange.second, rrange.first, rrange.second, param);
    
//     return summary;
// }

// std::pair<int, AlignSummary> assem::getMaxAlign(PathId query_pid, 
//                                                 PathId sbjct_pid, 
//                                                 intVec &types, 
//                                                 std::vector<IntPair> &rranges, 
//                                                 std::vector<IntPair> &qranges, 
//                                                 std::vector<int> &pos_vec,
//                                                 PathToAlnMap &path2aln_map, 
//                                                 Param &param)
// {
//     double max_score = 0;
//     AlignSummary max_summary;
//     int which = -1;
    
//     for ( size_t i = 0; i < rranges.size(); i++ ) {
//         if ( param.verbose) std::cout << "Block:" << i << "\n";
//         if ( types[i] == NOTYPE ) {
//             if ( param.verbose ) std::cout << "Invalid type\n";
//             continue;
//         }
//         AlignSummary summary = alignPathPair(query_pid, sbjct_pid, types[i], qranges[i], rranges[i], pos_vec, path2aln_map, param );
//         if ( summary.posrate > max_score ) {
//             max_score = summary.posrate;
//             max_summary = summary;
//             which = i;
//         }
//     }
//     return std::pair<int, AlignSummary>(which, max_summary);
// }

// bool __goodScore( int type, 
//                   AlignSummary &summary,
//                   Param &param )
// {
//     if ( type != LATCH && summary.posrate < param.merge_score )  {
//         if ( param.verbose ) std::cout << "Weak score:" << summary.posrate << "\n";
//         return false;
//     }
//     if ( type == LATCH && summary.posrate < param.latch_score ) {
//         if ( param.verbose ) std::cout << "Weak score:" << summary.posrate << "\n";
//         return false;
//     }
//     return true;
// }

// bool __goodLatch( AlignSummary &summary, ReadIdArray &sbjpreads, int direction, Param &param )
// {
//     if ( ! param.pair_flag ) {
//         if ( summary.length < param.latch_length ) {
//             if ( param.verbose ) std::cout << "Short overlapping latch\n";
//             return false;
//         }
//     } else {
//         if ( summary.length < param.pairend_overlap ) {
//             if ( param.verbose ) std::cout << "Short overlapping latch\n";
//             return false;
//         }
//         if ( summary.length < param.latch_length && (int)sbjpreads.size() < param.overlap_support ) {
//             if ( param.verbose ) std::cout << "Short overalp with weak supported latch\n";
//             return false;
//         }
//         if ( direction == POSITIVE && summary.s2se.first > summary.s1se.first ) {
//             if ( param.verbose ) std::cout << "Direction conflicts\n";
//             return false;
//         }
//         if ( direction == NEGATIVE && summary.s1se.first > summary.s2se.first ) {
//             if ( param.verbose ) std::cout << "Direction conflicts\n";
//             return false;
//         }
//     }
//     return true;
// }

// bool assem::__mergeBestPath( PathId query_pid, 
//                              PathId sbjct_pid, 
//                              intVec &types, 
//                              std::vector<IntPair> &rranges, 
//                              std::vector<IntPair> &qranges, 
//                              std::vector<int> &pos_vec,
//                              ReadIdArray &sbjpreads, 
//                              int direction, 
//                              PathToAlnMap &path2aln_map, 
//                              KmerToPathMap &pathid_map,
//                              BitString *bstrs,
//                              char *strands, 
//                              ReadId *pairs,
//                              PathId *used_reads,
//                              InvertedIndex &iindex,
//                              PathIdSet &merged_paths,
//                              Param &param)
// {
//     std::pair<int, AlignSummary> result = getMaxAlign(query_pid, sbjct_pid, types, rranges, qranges, pos_vec, path2aln_map, param);
//     if ( result.first == -1 ) return false;
    
//     int type = types[result.first];
//     AlignSummary summary = result.second;
    
//     if ( ! __goodScore(type, summary, param) ) return false;
//     if ( type == LATCH && ! __goodLatch(summary, sbjpreads, direction, param) ) return false;

//     if ( param.verbose ) std::cout << "Strong score:" << summary.posrate << "\n";
//     mergePathPair( query_pid, sbjct_pid, summary, merged_paths, path2aln_map, pathid_map, bstrs, strands, pairs, used_reads, iindex, param );
//     return true;
// }

void assem::printPositions( std::vector<int> &pos_vec )
{
    std::cout << "Position vector:\n";
    int i = 0;
    while ( i < (int)pos_vec.size() ) 
        std::cout << pos_vec[i++] << " ";
    std::cout << "\n";
}

// intVec assem::remakePositionVector( PathId qpid,
//                                     PathId spid, 
//                                     PathToAlnMap &path2aln_map, 
//                                     Param &param )
    
// {
// //     int mink = filter::minSameKmerCount( param.latch_length, param.kmer_size, 1-param.latch_score );
// //     if ( mink < 1 ) mink = 1;

//     KmerId *qkmers = path2aln_map[qpid]->getKmers();
//     KmerId *skmers = path2aln_map[spid]->getKmers();
//     size_t  qnkmer = path2aln_map[qpid]->getKmerCount();
//     size_t  snkmer = path2aln_map[spid]->getKmerCount();
//     //return makePosFlags(qnkmer, snkmer, qkmers, skmers, mink, param);
// return makePosFlags(qnkmer, snkmer, qkmers, skmers, 1, param);
// }

// bool assem::__tryMergePath( PathId &query_pid, 
//                             PathId  sbjct_pid,
//                             intVec pos_vec,
//                             PathToAlnMap &path2aln_map, 
//                             KmerToPathMap &pathid_map,
//                             PathIdSet &merged_paths,
//                             ReadId *pairs,
//                             char *strands,
//                             int nreads,
//                             BitString *bstrs,
//                             PathId *used_reads,
//                             InvertedIndex &iindex,
//                             Param &param )
// {
//     bool swapped = false;

//     /* Longer path becomes subject */
//     if ( path2aln_map[query_pid]->getKmerCount() > path2aln_map[sbjct_pid]->getKmerCount()) {
//         util::swap<PathId>(query_pid, sbjct_pid);
//         pos_vec = remakePositionVector( query_pid, sbjct_pid, path2aln_map, param );
//         swapped = true;
//     }

//     if ( param.verbose ) printPositions(pos_vec);

//     /* kmer sharing regions in reference */
//     intVec types;
//     std::vector<IntPair> rranges, qranges;
//     determineMatchRangesAndTypes( types, rranges, qranges, pos_vec, query_pid, sbjct_pid, pairs, nreads, path2aln_map, param );    

//     if ( rranges.size() == 0 ) {
//         if ( swapped ) util::swap<PathId>(query_pid, sbjct_pid); // reset
//         return false;
//     }
    
//     /* pair-end information */
//     ReadIdArray sbjpreads; // reads in subjct that paired with reads in query
//     int direction;
//     if ( param.pair_flag ) {
//         sbjpreads = getCommonReads( query_pid, sbjct_pid, path2aln_map, pairs, nreads, param );
//         direction = getDirection(sbjpreads, pairs, strands, param );
//     }

//     /* Scan best alignable block and merge path anchoring the block */
//     bool success = __mergeBestPath( query_pid, sbjct_pid, types, rranges, qranges, pos_vec, sbjpreads, direction, path2aln_map, pathid_map, bstrs, strands, pairs, used_reads, iindex, merged_paths, param );

//     /* Set query to sbjct to find successive merging paths */
//     if ( success ) {
//         query_pid = sbjct_pid;
//         if ( param.verbose ) std::cout << "NewQuery:" << query_pid << "\n";
//         return true;
//     }
    
//     /* Reset swapped paths if failed to merge */
//     if ( !success && swapped ) util::swap<PathId>(query_pid, sbjct_pid);

//     return false;
// }

bool assem::alignMiddle( AlignSummary &summary,
                         std::string sbjct,
                         std::string query,
                         //int filter_nkmer, 
                         int &type,
                         Param &param )
{
    int filter_nkmer = countSameKmers( sbjct, query, param.filter_kmer );
    //int min_nkmer = filter::minSameKmerCount( query.length(), param.kmer_size, 1-param.filter_score );
    int min_filter = filter::minSameKmerCount( query.length(), param.filter_kmer, 1-param.filter_score );
    if ( min_filter < 1 ) min_filter = 1;
    int fk = param.filter_kmer;
    if ( param.verbose ) std::cout << "#Same " << fk << "-mer:" << filter_nkmer << "\t#Min " << fk << "-mer:" << min_filter << "\n";
    if ( filter_nkmer < min_filter ) {
        if ( param.verbose ) std::cout << "Weak kmer sharing\n";
        return false;
    }
    
    SemiGlobalAlign aln(sbjct, query, ANCHOR_CENTER, param.gap_ext, param.gap_open);
    summary = aln.getSummary();
    if ( param.verbose ) {
        aln.printAlignment(std::cout);
        summary.print(std::cout);
    }
    
    type = getPathAlignType(summary);
    if ( type == MERGE ) {
        if ( summary.posrate < param.merge_score ) {
            if ( param.verbose ) std::cout << "Weak Merge\n";
            return false;
        } else {
            if ( param.verbose ) std::cout << "Strong Merge\n";
            return true;
        }
    }
    else if ( type == LATCH_LEFT || type == LATCH_RIGHT ) {
        if ( summary.posrate < param.latch_score ) {
            if ( param.verbose ) std::cout << "Weak Latch\n";
            return false;
        } else {
            if ( param.verbose ) std::cout << "Strong Latch\n";
            return true;
        }
    }
    return false;
}

bool assem::checkSimpleLatch( AlignSummary &summary, 
                       KmerId *qkmers, 
                       KmerId *rkmers, 
                       int qnkmer,
                       int rnkmer, 
                       IntPair &qrange, 
                       IntPair &rrange, 
                       Param   &param )
{
    adjustRanges(qrange, rrange, qnkmer, rnkmer, param);

    std::string query = biostr::getSequenceString(&qkmers[qrange.first], qrange.second-qrange.first+1, param.kmer_size);
    std::string sbjct = biostr::getSequenceString(&rkmers[rrange.first], rrange.second-rrange.first+1, param.kmer_size);
    
    if ( query.size() != sbjct.size() ) {
        if ( param.verbose ) std::cout << "query size is not equal to sbjct size\n";
        return false;
    }
    
    if ( param.verbose ) {
        std::cout << "query:" << query << "\n";
        std::cout << "sbjct:" << sbjct << "\n";
    }
    summary = compareBySubstitution(query, sbjct, qnkmer, rnkmer, qrange.first, qrange.second, rrange.first, rrange.second, param);
    //if ( param.verbose ) summary.print(std::cout);
    
    if ( summary.posrate < param.latch_score ) {
        if ( param.verbose ) std::cout << "Weak latch by substitution\n";
        return false;
    }

    if ( param.pair_flag ) {
        if ( summary.length < param.pairend_overlap ) {
            if ( param.verbose ) std::cout << "Short latch overlap\n";
            return false;
        }
    }
    else {
        if ( summary.length < param.latch_length ) {
            if ( param.verbose ) std::cout << "Short latch overlap\n";
            return false;
        }
    }
    
    if ( param.verbose ) std::cout << "Strong latch by substitution\n";
    return true;
}

    
bool assem::checkLatchAlign( AlignSummary &summary, 
                             std::string &query,
                             std::string &sbjct,
                             int anchor,
                             Param   &param )
{
    if ( anchor == ANCHOR_LEFT ) {
        SemiGlobalAlign aln(sbjct, query, ANCHOR_LEFT, param.gap_ext, param.gap_open);
        summary = aln.getSummary();
        if ( param.verbose ) {
            aln.printAlignment(std::cout);
            summary.print(std::cout);
        }
    } else {
        SemiGlobalAlign aln(sbjct, query, ANCHOR_RIGHT, param.gap_ext, param.gap_open);
        summary = aln.getSummary();
        if ( param.verbose ) {
            aln.printAlignment(std::cout);
            summary.print(std::cout);
        }
    }
    
    if ( summary.posrate < param.latch_score ) {
        if ( param.verbose ) std::cout << "Weak Latch\n";
        return false;
    } else {
        if ( param.pair_flag ) {
            if ( summary.length < param.pairend_overlap ) {
                if ( param.verbose ) std::cout << "Short overlap\n";
                return false;
            }
        }
        else {
            if ( summary.length < param.latch_length ) {
                if ( param.verbose ) std::cout << "Short overlap\n";
                return false;
            }
        }
        if ( param.verbose ) std::cout << "Strong Latch\n";
        return true;
    }
}

bool assem::alignEnd( AlignSummary &summary,
                      std::string sbjct,
                      std::string query,
                      PathId query_pid,
                      PathId sbjct_pid,
                      std::set<IntPair> &bad_range,
                      PathToAlnMap &path2aln_map, 
                      Param &param, 
                      int anchor )
    
{
    KmerId *qkmers = path2aln_map[query_pid]->getKmers();
    KmerId *skmers = path2aln_map[sbjct_pid]->getKmers();
    size_t  qnkmer = path2aln_map[query_pid]->getKmerCount();
    size_t  snkmer = path2aln_map[sbjct_pid]->getKmerCount();
    

    std::vector<int> pvec = makePosFlags(qnkmer, snkmer, qkmers, skmers, 1, param);
    if ( param.verbose ) printPositions(pvec);

    IntPair srange;
    anchor == ANCHOR_LEFT ? 
        srange = getMatchRangeLeft(pvec, param) :
        srange = getMatchRangeRight(pvec, param) ;
    if ( srange.first == -1 || srange.second == -1 ) return false;
    IntPair qrange = IntPair( pvec[srange.first], pvec[srange.second] );
    
    if ( bad_range.find(qrange) != bad_range.end() ) {
        if ( param.verbose ) std::cout << "Already scanned failed range\n";
        return false;
    }

    if ( param.verbose ) {
        std::cout << "ss:" << srange.first << "\tse:" << srange.second << "\n";
        std::cout << "qs:" << qrange.first << "\tqe:" << qrange.second << "\n";
    }
    
    bool success = false;

    int type = getMatchType( qrange, srange, qkmers, skmers, qnkmer, snkmer, param );
    if ( param.verbose ) std::cout << "Anchor:" << anchor << "\tType:" << type << "\n";
    if ( anchor == ANCHOR_LEFT ) {
        //if ( type == LATCH_RIGHT_EASY || type == LATCH_RIGHT_DIFF || type == NOTYPE ) 
        if ( type != LATCH_LEFT_EASY && type != LATCH_LEFT_DIFF ) {
            if ( param.verbose ) std::cout << "Not a latch\n";
            return false;
        }
    }
    if ( anchor == ANCHOR_RIGHT ) {
        //if ( type == LATCH_LEFT_EASY || type == LATCH_LEFT_DIFF || type == NOTYPE ) 
        if ( type != LATCH_RIGHT_EASY && type != LATCH_RIGHT_DIFF )  {
            if ( param.verbose ) std::cout << "Not a latch\n";
            return false;
        }
    }

    std::string ssbj = biostr::getSequenceString(skmers+srange.first, srange.second-srange.first+1, param.kmer_size);
    std::string qsbj = biostr::getSequenceString(qkmers+qrange.first, qrange.second-qrange.first+1, param.kmer_size);
    if ( param.verbose ) {
        std::cout << "sbjct substr:" << ssbj << "\n";
        std::cout << "query substr:" << qsbj << "\n";
    }

    int filter_nkmer = countSameKmers( ssbj, qsbj, param.filter_kmer );
    int min_filter = filter::minSameKmerCount( param.latch_length, param.filter_kmer, 1-param.filter_score );
    if ( min_filter < 1 ) min_filter = 1;

    int fk = param.filter_kmer;
    if ( param.verbose ) std::cout << "#Same " << fk << "-mer:" << filter_nkmer << "\t#Min " << fk << "-mer:" << min_filter << "\n";
    if ( filter_nkmer < min_filter ) {
        if ( param.verbose ) std::cout << "Weak kmer sharing\n";
        return false;
    }

    if ( param.verbose ) std::cout << "Checking simple latch\n";
    success = checkSimpleLatch( summary, qkmers, skmers, qnkmer, snkmer, qrange, srange, param );
    if ( success ) return true;

    if ( param.verbose ) std::cout << "Checking latch based on alignment\n";
    success = checkLatchAlign( summary, query, sbjct, anchor, param );
    if ( success ) return true;  

    bad_range.insert(qrange);
    return false;
}


// std::pair<bool, int> assem::__tryMergePathPair( PathId &sbjct_pid,
//                                                 PathId query_pid, 
//                                                 //int filter_kmers,
//                                                 PathToAlnMap &path2aln_map, 
//                                                 //KmerToPathMap &pathid_map,
//                                                 PathIdSet &merged_paths,
//                                                 ReadId *pairs,
//                                                 char *strands,
//                                                 int nreads,
//                                                 BitString *bstrs,
//                                                 PathId *used_reads,
//                                                 InvertedIndex &iindex,
//                                                 Param &param,
//                                                 int join_type )
// {
//     bool swapped = false;

//     /* Longer path becomes subject */
//     if ( path2aln_map[query_pid]->getKmerCount() > path2aln_map[sbjct_pid]->getKmerCount()) {
//         util::swap<PathId>(query_pid, sbjct_pid);
//         swapped = true;
//     }

// //     IntVec pos_vec = remakePositionVector( query_pid, sbjct_pid, path2aln_map, param );
// //     if ( pos_vec.size() == 0 ) {
// //         if ( swapped ) util::swap<PathId>(query_pid, sbjct_pid);
// //         return false;
// //     }

// //     /* kmer sharing regions in reference */
// //     intVec types;
// //     std::vector<IntPair> rranges, qranges;
// //     determineMatchRangesAndTypes( types, rranges, qranges, pos_vec, query_pid, sbjct_pid, pairs, nreads, path2aln_map, param );    


//     std::string query = biostr::getSequenceString(path2aln_map[query_pid]->getKmers(), path2aln_map[query_pid]->getKmerCount(), param.kmer_size);
//     std::string sbjct = biostr::getSequenceString(path2aln_map[sbjct_pid]->getKmers(), path2aln_map[sbjct_pid]->getKmerCount(), param.kmer_size);
    
//     assert(query.size() > 0);
//     assert(sbjct.size() > 0);
//     if ( param.verbose ) {
//         std::cout << "\nsbjct:" << sbjct_pid << "\t" << sbjct << "\n";
//         std::cout << "query:" << query_pid << "\t" << query << "\n";
//     }


//     std::set<IntPair> bad_range;
//     AlignSummary summary;
//     int anchor = ANCHOR_CENTER;
//     int type = NOTYPE;
//     bool found = false;

//     if ( join_type == PATHMERGE || join_type == MERGELATCH ) {
//         if ( param.verbose ) std::cout << "Check middle\n";
//         found = alignMiddle(summary, sbjct, query, type, param );
//     }

// //     if ( found ) {
// //         if ( ( direction == POSITIVE && anchor == LATCH_LEFT  ) || 
// //              ( direction == NEGATIVE && anchor == LATCH_RIGHT ) ) { 
// //             if ( param.verbose ) std::cout << "\tDirectionality conflicts\n"; 
// //             found = false;
// //         }
// //     }

//     if ( join_type == PATHLATCH || join_type == MERGELATCH ) {
//         if ( !found ) {
//             if ( param.verbose ) std::cout << "Check left\n";
//             anchor = ANCHOR_LEFT;
//             found = alignEnd(summary, sbjct, query, query_pid, sbjct_pid, bad_range, path2aln_map, param, anchor );
//         }
//         if ( !found ) {
//             if ( param.verbose ) std::cout << "Check right\n";
//             anchor = ANCHOR_RIGHT;
//             found = alignEnd(summary, sbjct, query, query_pid, sbjct_pid, bad_range, path2aln_map, param, anchor );
//         }
//     }

//     if ( !found )  {
//         if ( swapped ) util::swap<PathId>(query_pid, sbjct_pid);
//         return std::pair<bool, int>(false, type);
//     }

//     /* pair-end information */
//     ReadIdArray sbjpreads; // reads in subjct that paired with reads in query
//     int direction = NONSENSE;
//     if ( param.pair_flag ) {
//         sbjpreads = getCommonReads( query_pid, sbjct_pid, path2aln_map, pairs, nreads, param );
//         direction = getDirection(sbjpreads, pairs, strands, param );
//     }
                
//     /* Check direction conflict */
//     if ( ( direction == POSITIVE && anchor == ANCHOR_LEFT  ) || 
//          ( direction == NEGATIVE && anchor == ANCHOR_RIGHT ) ) { 
//         if ( param.verbose ) std::cout << "\tDirectionality conflicts\n"; 
//         if ( swapped ) util::swap<PathId>(query_pid, sbjct_pid);
//         return std::pair<bool, int>(false, type);
//     }
    
//     /* Scan best alignable block and merge path anchoring the block */
//     //bool success = __mergeBestPath( query_pid, sbjct_pid, types, rranges, qranges, pos_vec, sbjpreads, direction, path2aln_map, pathid_map, bstrs, strands, pairs, used_reads, iindex, merged_paths, param );

//     //bool success = mergePathPair( query_pid, sbjct_pid, summary, merged_paths, path2aln_map, pathid_map, bstrs, strands, pairs, used_reads, iindex, param );
//     bool success = mergePathPair( query_pid, sbjct_pid, summary, merged_paths, path2aln_map, bstrs, strands, pairs, used_reads, iindex, param );
//     if ( !success ) return std::pair<bool, int>(false, type);

//     //query_pid = sbjct_pid;

//     return std::pair<bool, int>(true, type);

// //     /* Set query to sbjct to find successive merging paths */
// //     if ( success ) {
// //         query_pid = sbjct_pid;
// //         if ( param.verbose ) std::cout << "NewQuery:" << query_pid << "\n";
// //         return true;
// //     }
    
// //     /* Reset swapped paths if failed to merge */
// //     if ( !success && swapped ) util::swap<PathId>(query_pid, sbjct_pid);

// //     return false;
// }

bool assem::alreadyCompared( PathId sbjct_pid,
                             PathId query_pid, 
                             std::tr1::unordered_map<PathId, PathIdList> &compared_paths )
{
    if ( compared_paths.find(sbjct_pid) == compared_paths.end() ) return false;
    for ( PathIdList::iterator it = compared_paths[sbjct_pid].begin(); it != compared_paths[sbjct_pid].end(); ++it )
        if ( *it == query_pid ) return true;
    return false;
}

// bool assem::mergeSimilarPath( PathId &sbjct_pid, 
//                               //PathIdSet &skip_pids,
//                               PathToAlnMap &path2aln_map, 
//                               KmerToPathMap &pathid_map,
//                               PathIdSet &merged_paths,
//                               std::tr1::unordered_map<PathId, PathIdList> &compared_paths,
//                               std::tr1::unordered_map<PathId, bool> &updated_paths,
//                               ReadId *pairs,
//                               char *strands, 
//                               int nreads,
//                               BitString *bstrs,
//                               PathId *used_reads,
//                               InvertedIndex &iindex,
//                               Param &param )
// {
//     /* List is savend in order by the number of sharing kmers */
//     //PosPathPairList pos_paths = getMergiblePaths(sbjct_pid, path2aln_map, pathid_map, merged_paths, pairs, nreads, param);    
//     double t0 = mytime(); 

//     PathId curr_pid = sbjct_pid;
//     bool latch = false;

//     std::multimap<int, PathId> sim_paths = getSimilarPaths( sbjct_pid, path2aln_map, pathid_map, merged_paths, pairs, nreads, param );

//     if ( param.verbose ) std::cout << "# Similar paths:" << sim_paths.size() << "\t" << mytime()-t0 << "\n";
//     if ( sim_paths.size() ==  0 ) return false;

//     std::multimap<int, PathId>::reverse_iterator it;
//     for ( it = sim_paths.rbegin(); it != sim_paths.rend(); ++it ) {
//         if ( alreadyCompared(sbjct_pid, it->second, compared_paths) ) {
//             if ( param.verbose ) std::cout << "Already compared pair:" << sbjct_pid << "\t" << it->second << "\n";
//             continue;
//         }
        
//         std::pair<bool, int> result = __tryMergePathPair( sbjct_pid, it->second, it->first, path2aln_map, pathid_map, merged_paths, pairs, strands, nreads, bstrs, used_reads, iindex, param );
//         /* No further comparison in case of weak merge 
//            Weak latch will be compared again since extended path can merge it */
//         if ( !result.first && !latchableType(result.second) )
//             compared_paths[sbjct_pid].push_back(it->second);
//         if ( result.first && latchableType(result.second) ) latch = true;
//         /* Longer path found */
//         if ( result.first && sbjct_pid != curr_pid ) 
//             return true;
//     }
//     if ( latch ) return true; // search again because of path extension
//     return false;
// }


// bool assem::searchMostSimilarPath( PathId &query_pid, 
//                             PathId &sbjct_pid, 
//                             AlignSummary &summary,
//                             PathToAlnMap &path2aln_map, 
//                             KmerToPathMap &pathid_map,
//                             PathIdSet &merged_paths,
//                             ReadId *pairs,
//                             int nreads,
// 							Param &param )
// {
//     KmerId *kmers = path2aln_map[query_pid]->getKmers();
//     size_t  nkmer = path2aln_map[query_pid]->getKmerCount();    
//     PathIdSet self; 
//     self.insert(query_pid);
    
//     if ( param.verbose ) 
//         std::cout << "**Query:" << query_pid << "\t" << biostr::getSequenceString(kmers, nkmer, param.kmer_size) << "\n";
    
//     //int mink = getMinKmerCount( nkmer+param.kmer_size-1, param );
//     int mink = filter::minSameKmerCount( param.latch_length, param.kmer_size, 1-param.filter_score );
//     if ( mink < 1 ) mink = 1;
//     //if ( param.verbose ) std::cout << "Min kmer:" << mink << "\n";
//     PosPathPairList pos_paths = searchSimilarPaths(self, kmers, nkmer, path2aln_map, pathid_map, merged_paths, mink, param);
//     if ( pos_paths.size() == 0 ) {
//         return false;
//     }
//     if ( param.verbose ) std::cout << "#Similar paths:" << pos_paths.size() << "\n";
    
//     bool success = pickMergePath(pos_paths, query_pid, sbjct_pid, summary, path2aln_map, pathid_map, pairs, nreads, param);
//     if ( !success ) 
//         return false;

//     return true;
// }

void assem::makePathIdMap( KmerToPathMap &pathid_map,
                           PathToAlnMap &path2aln_map )
{
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
        PathId pid = it->first;
        SpaPath *spath = path2aln_map[pid];
        KmerId *kmers = spath->getKmers();
        size_t nkmer = spath->getKmerCount();
        addPathIds( pid, kmers, nkmer, pathid_map );
    }
}

void assem::makePathKmersMap( PathKmersMap &pkmers_map,
                              PathToAlnMap &path2aln_map,
                              int node_size,
                              int word_size )
{
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
        std::string sbjct = biostr::getSequenceString(path2aln_map[it->first]->getKmers(), path2aln_map[it->first]->getKmerCount(), node_size);
        KmerArray kmers = biostr::getKmers( sbjct, word_size );
        pkmers_map.insert( std::pair<PathId, KmerArray>( it->first, kmers ) );
    }
}

// PathIdSet assem::linkPaths(PathIdSet &path_ids,
//                            PathToAlnMap &path2aln_map, 
//                            InvertedIndex &iindex,
//                            PathId *used_reads,
//                            BitString *bstrs,
//                            char *strands,
//                            ReadId *pairs,
//                            int nreads,
//                            Param &param,
//                            int join_type )
// {

//     PathIdList success_paths;

//     size_t npath = path_ids.size();
//     if ( param.verbose ) std::cout << "#Paths to merge/extend:" << npath << "\n";
    
// //     /* Make k-mer to path-id membership map */
// //     KmerToPathMap pathid_map;
// //     makePathIdMap(pathid_map, path2aln_map );

//     /* Make a path length map */
//     //PathLengthMap plen_map = getPathLengths( path2aln_map, path_ids );
//     PathLengthMap plen_map = getPathLengths( path2aln_map );

// //     PathKmersMap pkmers_map;
// //     makePathKmersMap(pkmers_map, path2aln_map, param.kmer_size, param.filter_kmer);
    
//     /* auxilary variables for progress log */
//     double pratio = 1.0;
//     int    nproc  = 0;
//     double lt0    = mytime();

//     /* Keep track of merged paths */
//     PathIdSet merged_paths; 


// //     /* Keep track of compared path pairs */
//     std::tr1::unordered_map<PathId, PathIdList> compared_paths;
//     std::tr1::unordered_map<PathId, bool> updated_paths;
//     std::tr1::unordered_map<PathId, PathIdList> path_bins;

// //     std::tr1::unordered_map<PathId, bool> joined_paths;
// //     for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) 
// //         joined_paths.insert( std::pair<PathId, bool>( it->second, false ) );

//     //for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) 
//         //updated_paths.insert( std::pair<PathId, bool>( it->second, false) );

//     /* cluster centers */
//     PathIdSet clusters;
//     typedef std::tr1::unordered_map<KmerId, int> KmerCountMap;
//     typedef std::tr1::unordered_map<PathId, KmerCountMap> PathKmerCountMap;

//     PathKmerCountMap pkc_map;

//     KmerPathIndex path_index;
//     for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
//         KmerArray kmers = KmerArray( (*it).second->getKmers(),  (*it).second->getKmers() + (*it).second->getKmerCount() );
//         path_index.init(kmers);
//     }

//     //KmerPathTable path_index(param.filter_kmer);
//     if (param.verbose ) std::cout << "Size:" << path_index.getSize() << "\n";


//     /*---------------------------------------------------------------
//      * Starting from the longest path (path map is already in order).
//      * If a similar path found, merge/extend the pair. 
//      * Once the pair merged, repeat the process to further merging
//      *--------------------------------------------------------------*/
//     //size_t length_cutoff = 60;

//     size_t success = 0;
//     for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
//         //PathId sbjct_pid = it->second;
//         PathId query_pid = it->second;
//         //if ( it->first < length_cutoff ) break;
//         //if ( path_ids.find(sbjct_pid) == path_ids.end() ) continue;
//         //if ( path_ids.find(query_pid) == path_ids.end() ) continue;

//         KmerArray query_kmers = KmerArray( path2aln_map[query_pid]->getKmers(),  path2aln_map[query_pid]->getKmers() + path2aln_map[query_pid]->getKmerCount() );
//         //KmerArray query_kmers = pkmers_map[query_pid];
        
// //         /** empty clusters */
// //         if ( clusters.size() == 0 ) {         
// //         //if ( path_index.getSize() == 0 ) {
// //             clusters.insert( query_pid ); 
            
// // //             KmerCountMap count_map;
// // //             setKmerMap(count_map, query_kmers);
// // //             pkc_map[query_pid] = count_map;

// //             path_index.add( query_kmers, query_pid );

// //             continue;
// //         }

//         /* Print progress status at every percent of path proceeded */
//         nproc++;
//         if ( nproc/(double)npath*100 >= pratio ) {
//             fprintf( stdout, "%d/%d: %.2f sec\n", nproc, int(npath), mytime()-lt0);
//             pratio += 1;
//         }
        
//         /* Skip if merged previously */
//         //if ( merged_paths.find(sbjct_pid) != merged_paths.end() ) continue;
//         //if ( joined_paths[query_pid] ) continue;
        
//         if ( param.verbose ) std::cout << "\nQuery:" << query_pid << "\tLength:" << it->first << "\n";
        
//         double tm0 = mytime();

//         double tt0 = mytime();
//         //std::string query = biostr::getSequenceString(path2aln_map[query_pid]->getKmers(), path2aln_map[query_pid]->getKmerCount(), param.kmer_size);        
//         //if ( param.verbose ) std::cout << "string time:" << mytime() - tt0  << "\n";
//         tt0 = mytime();
//         //int min_filter = filter::minSameKmerCount( query.length(), param.filter_kmer, 1-param.filter_score );
//         size_t nkmers = path2aln_map[query_pid]->getKmerCount();
//         int min_filter = filter::minSameKmerCount( nkmers+param.kmer_size-1, param.filter_kmer, 1-param.filter_score );
//         if ( min_filter < 1 ) min_filter = 1;
        
//         if (param.verbose) std::cout << "Min k-mer:" << min_filter << "\n";//<< "\ttime:" << mytime()-tt0 << "\n";

//         tt0 = mytime();
// //         KmerArray query_kmers = pkmers_map[query_pid];
// //         std::multimap<size_t, PathId> share_kmers;
// //         double avgtime = 0;
// //         size_t i = 0;
// //         for ( PathKmerCountMap::iterator it = pkc_map.begin(); it != pkc_map.end(); ++it ) {
// //             double ttt0 = mytime();
// //             KmerCountMap count_map = it->second;
// //             size_t same_kmers = countSameKmers( count_map, query_kmers );
// //             avgtime += ( mytime()-ttt0 );
// //             if ( same_kmers >= min_filter ) share_kmers.insert( std::pair<size_t, PathId>( same_kmers, it->first ) );
// //             i++;
// //         }
// //         if ( param.verbose &&  i > 0 ) std::cout << "avg. time for shared kmer search:" << avgtime/i << "\n";
// //         if ( param.verbose ) std::cout << "# similar paths:" << share_kmers.size() << "/" << clusters.size() << "\ttime:" << mytime()-tt0 << "\n";
        
//         PathFreqMap path_freq;
// //         for ( PathIdSet::iterator it = clusters.begin(); it != clusters.end(); ++it )
// //             path_freq[*it] = 0;

//         path_index.findPaths( path_freq, query_kmers );
// //         for ( PathFreqMap::iterator pt = path_freq.begin(); pt != path_freq.end(); )
// //             if ( pt->second < min_filter ) path_freq.erase(pt++);
// //             else ++pt;
       
//         std::multimap<size_t, PathId> freq_map = util::sortByValue<PathId, size_t>(path_freq);

//         if ( param.verbose ) std::cout << "similar paths:" << freq_map.size() << "\tsearch time:" << mytime() -tt0 << "\n";
                
//         tt0 = mytime();

//         //size_t count = 0;
//         bool found = false;
//         PathId target = query_pid;
        
//         //if ( share_kmers.size() ) {
//         if ( freq_map.size() ) {
//             //for ( std::multimap<size_t, PathId>::reverse_iterator it = share_kmers.rbegin(); it != share_kmers.rend(); ++it ) {
//             //std::multimap<size_t, PathId>::reverse_iterator it = share_kmers.rbegin();
//             std::multimap<size_t, PathId>::reverse_iterator it = freq_map.rbegin();
            
//             if ( (int)it->first >= min_filter ) {
//                 //++count;
//                 //if ( min_filter <= it->first ) {
//                 path_bins[it->second].push_back( query_pid );
//                 target = it->second;
//                 found = true;
//                 //break;
//                 //}
//                 //}

//                 success++;
//                 if ( param.verbose ) std::cout << "success:" << target << "\tsame:" << it->first << "\tbin size:" << path_bins[target].size() << "\t# clusters:" << clusters.size() << "\t# kmers:" << path_index.getSize() << "\t time:" << mytime()-tm0 << " sec\n";
//                 continue;
//             } else found = false;
//         }


//         if ( !found ) {
//             clusters.insert( query_pid ); 
// //             KmerCountMap count_map;
// //             setKmerMap(count_map, query_kmers);
// //             pkc_map[query_pid] = count_map;

//             path_index.add( query_kmers, query_pid );

//             if ( param.verbose ) std::cout << "# update time:" << mytime()-tt0 << "\n";
//             if ( param.verbose ) std::cout << "fail:" <<  "\t# clusters:" << clusters.size() << "\t# kmers:" << path_index.getSize()  << "\t time:" << mytime()-tm0 << " sec\n";
//         }



// //         if ( found ) {
// //             if ( param.verbose ) std::cout << "success:" << target << "\tbin size:" << path_bins[target].size() << "\t# clutsters:" << clusters.size() << "\t time:" << mytime()-tm0 << " sec\n";
// //         } else {
// //             if ( param.verbose ) std::cout << "fail:" <<  "\t# clusters:" << clusters.size()  << "\t time:" << mytime()-tm0 << " sec\n";
// //         }


// //         std::string sbjct = biostr::getSequenceString(path2aln_map[sbjct_pid]->getKmers(), path2aln_map[sbjct_pid]->getKmerCount(), param.kmer_size);
// //         //KmerArray sbjct_kmers = biostr::getKmers( sbjct, param.filter_kmer );
// //         KmerArray sbjct_kmers = pkmers_map[sbjct_pid];
// //         std::tr1::unordered_map<KmerId, int> count_map;
// //         setKmerMap(count_map, sbjct_kmers);

// //         size_t count = 0;
// //         double tm0 = mytime();
// //         for ( PathLengthMap::reverse_iterator jt = it; jt != plen_map.rend(); ++jt ) {
// //             if ( jt == it ) continue;
// //             double tic = mytime();
// //             PathId query_pid = jt->second;
// //             //if ( merged_paths.find(query_pid) != merged_paths.end() ) continue;
// //             if ( joined_paths[query_pid] ) continue;
            
// //             count++;

// //             std::string query = biostr::getSequenceString(path2aln_map[query_pid]->getKmers(), path2aln_map[query_pid]->getKmerCount(), param.kmer_size);
// //             assert(query.size() >0);
// //             int min_filter = filter::minSameKmerCount( query.length(), param.filter_kmer, 1-param.filter_score );
// //             //size_t same_kmers = countSameKmers( pathid_map, path2aln_map, query_pid, sbjct_pid, param.kmer_size );
// //             //size_t same_kmers = countSameKmers( pathid_map, path2aln_map, query_pid, sbjct_pid, param.filter_kmer );
// //             //size_t same_kmers = countSameKmers( sbjct, query, param.filter_kmer );
// //             //KmerArray query_kmers = biostr::getKmers( query, param.filter_kmer );
// //             KmerArray query_kmers = pkmers_map[query_pid];
// //             size_t same_kmers = countSameKmers( count_map, query_kmers );

// //             if ( same_kmers >= min_filter ) {
// //                 path_bins[sbjct_pid].push_back( query_pid );
// //                 joined_paths[query_pid] = true;
// //             }
// // // //             if ( same_kmers == 0 ) {
// // // //                 if ( param.verbose ) std::cout << "No kmer sharing between " << sbjct_pid << " and " << query_pid << "\n";
// // // //                 continue;
// // // //             }
// // //             if ( same_kmers < min_filter ) {
// // //                 //if ( param.verbose ) std::cout << "Weak kmer sharing between " << sbjct_pid << " and " << query_pid << "\t" << same_kmers << "/" << min_filter << "\n";
// // //                 //if ( param.verbose ) std::cout << "Weak\t" << sbjct_pid << " vs " << query_pid << "\t" << same_kmers << "/" << min_filter << "\n";
// // //                 if ( param.verbose && count%100 == 0 ) std::cout << ".";
// // //                 continue;
// // //             } 
// // // //             else { 
// // // //                 if (param.verbose ) 
// // // //                     std::cout << "\n";
// // // //                 //<< "sbjct:" << sbjct_pid << "\tquery:" << query_pid << "\tsame:" << same_kmers << "\tmin:" << min_filter << "\n";
// // // //             }

// // //             //__tryMergePathPair( query_pid, sbjct_pid, filter_kmers, path2aln_map, pathid_map, merged_paths, pairs, strands, nreads, bstrs, used_reads, iindex, param );
// // //             std::pair<bool, int> result = __tryMergePathPair( query_pid, sbjct_pid, path2aln_map, pathid_map, merged_paths, pairs, strands, nreads, bstrs, used_reads, iindex, param, join_type );
// // //             if ( result.first ) {
// // //                 success_paths.push_back(sbjct_pid);
// // //                 // subjct might be changed after merging
// // //                 sbjct_kmers = biostr::getKmers( sbjct, param.filter_kmer );
// // //                 setKmerMap(count_map, sbjct_kmers);
// // //             }

// //             //double toc = mytime();
// //             //if ( param.verbose ) std::cout << "merge test time:" << toc-tic << " sec\n";
// //         }
// //         //if ( param.verbose ) std::cout << "\n# tested paths:" << count << "\t# success so far:" << success_paths.size() << "\tmerge time:" << mytime()-tm0 << " sec\n";
// //         if ( param.verbose ) std::cout << "# comparison:" << count << "\tbin size:" << path_bins[sbjct_pid].size() << "\t time:" << mytime()-tm0 << " sec\n";


// //         /* Recursively merge/extend paths */
// //         while(true) {
// // //             PathIdSet skip_pids;
// // //             for ( PathLengthMap::reverse_iterator jt = it; jt != plen_map.rend(); ++jt ) 
// // //                 skip_pids.insert(jt->second);

// //             //bool success = mergeSimilarPath(sbjct_pid, skip_pids, path2aln_map, pathid_map, merged_paths, pairs, strands, nreads, bstrs, used_reads, iindex, param);
// //             bool success = mergeSimilarPath(sbjct_pid, path2aln_map, pathid_map, merged_paths, compared_paths, updated_paths, pairs, strands, nreads, bstrs, used_reads, iindex, param);
// //             if ( !success ) break;
// //         }
//     }

//     if ( param.verbose ) std::cout << "Success:" << success << "\n";

//     /* Remove merged paths from the path list */
//     dropMergedPaths(merged_paths, path2aln_map);
    
//     if ( param.verbose ) printPaths(path2aln_map, param);

//     return PathIdSet( success_paths.begin(), success_paths.end() );
// }

bool assem::comparePathSize( const PathSize &one, const PathSize &two )
{
    //return one.size < two.size ? true : false;
    //if ( one.size < two.size) return true;
    //return false;
    return one.size < two.size;
}

void assem::mergeClusters( PathBins &path_bins,
                           PathToAlnMap &path2aln_map, 
                           BitString *bstrs,
                           PathId *used_reads,
                           char *strands,
                           ReadId *pairs,
                           Param &param )
{
    if ( param.verbose ) std::cout << "# Total bins:" << path_bins.size() << "\n";

    size_t no = 0;
    for ( PathBins::iterator it = path_bins.begin(); it != path_bins.end(); ++it ) {
        PathId sbjct_pid = it->first;
        PathAlignPairList path_aligns = it->second;
        assert( path2aln_map.find(sbjct_pid) != path2aln_map.end() );

        // copy old path
        //SpaPath *old_path = new SpaPath(*path2aln_map[sbjct_pid]);
        SpaPath *old_path = new SpaPath();
        *old_path = *path2aln_map[sbjct_pid];
        SpaPath *new_path = path2aln_map[sbjct_pid];

        ReadIdArray old_reads = ReadIdArray( old_path->getReads(), old_path->getReads() + old_path->getReadCount() );

        /* singleton */
        if ( path_aligns.size() == 0 ) continue;

        if ( param.verbose ) std::cout << "Cluster:" << ++no << "\t";
        if ( param.verbose ) std::cout << "center:" << sbjct_pid << "\t#members:" << path_aligns.size() << "\n";

        if ( param.verbose ) {
            std::cout << "old str:" << biostr::getSequenceString(new_path->getKmers(), new_path->getKmerCount(), param.kmer_size) << "\n";
            std::cout << "old nread:" << new_path->getReadCount() << "\n";
        }

        if ( param.verbose ) {
            std::cout << "Old alignment - Sbjct:" << sbjct_pid << "\n";
            new_path->printAlignment(std::cout, 100, bstrs);
            
            for ( PathAlignPairList::iterator jt = path_aligns.begin(); jt != path_aligns.end(); ++jt ) {
                PathId pid = jt->first;
                std::cout << "\nOld alignment - Query:" << pid << "\tGap:" << jt->second.lgap.first << "\t" << jt->second.egap.first << "\n";
                path2aln_map[pid]->printAlignment(std::cout, 100, bstrs);
            }
        }

        double tm0 = mytime();
        new_path->mergeCluster( path_aligns, path2aln_map, bstrs, param );
        if ( param.verbose ) std::cout << "A cluster merged:" << mytime()-tm0 << "\n";

        if ( param.verbose ) {
            std::cout << "new str:" << biostr::getSequenceString(new_path->getKmers(), new_path->getKmerCount(), param.kmer_size) << "\n";
            std::cout << "new nread:" << new_path->getReadCount() << "\n";
        }
        
        int error = 0;
        MSA msa(*new_path, sbjct_pid, used_reads, bstrs, strands, pairs, PILE|FLAT|TRIM|CODON, error, param);
        if ( error > 0 ) {
            if ( param.verbose ) std::cout << "MSA fail:\tcode:" << error << "\n";
            ReadIdArray new_reads = ReadIdArray( new_path->getReads(), new_path->getReads() + new_path->getReadCount() );            
            ReadIdArray invalids = Set::Difference<ReadId>( new_reads, old_reads, false );
            for ( size_t i = 0; i < invalids.size(); i++ ) used_reads[invalids[i]] = NOT_PATH;
            
            delete new_path;
            path2aln_map[sbjct_pid] = old_path;
            continue;
        }
        else delete old_path;
        
        if ( param.verbose ) std::cout << "Consensus:" << new_path->getConsensusString() << "\n";
        if ( param.verbose ) std::cout << "\nAlignment of cluster\n";
        if ( param.verbose ) new_path->printAlignment(std::cout, 100, bstrs);
    }
}


void assem::clusterPaths( PathToAlnMap &path2aln_map, 
                          InvertedIndex &iindex,
                          PathId *used_reads,
                          BitString *bstrs,
                          char *strands,
                          ReadId *pairs,
                          int nreads,
                          Param &param )
{
    double ts = mytime();

    setbuf(stdout, NULL); // no buffering
    if ( param.verbose) std::cout << "\nClustering paths\n";
    if ( param.verbose) std::cout << "#Total paths:" << path2aln_map.size() << "\n";

    
    /* Initialize Index */
    /* Avoid large hash miss */
    //size_t count = 1;
    KmerPathIndex path_index;
    double t0 = mytime();
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
        std::string query = (*it).second->getConsensusString();
        assert(query.size()>0);
        KmerArray kmers = biostr::getKmers( query, param.filter_kmer );
        path_index.init(kmers);
    }
    if ( param.verbose ) std::cout << "Path index init:" << mytime()-t0 << "\n";

    /* Keep track of merged paths */
    PathIdSet merged_paths; 

    /* Bins */
    PathBins path_bins;

    /* auxilary variables for progress log */
    double pratio = 1.0;
    int    nproc  = 0;
    double lt0    = mytime();
    size_t npath  = path2aln_map.size();

    /* Make a path length map */
    PathLengthMap plen_map = getPathLengths( path2aln_map );
    for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        PathId query_pid = it->second;

        std::string query = path2aln_map[query_pid]->getConsensusString();
        assert(query.size() > 0);
        KmerArray query_kmers = biostr::getKmers( query, param.filter_kmer );

        /* Print progress status at every percent of path proceeded */
        nproc++;
        if ( nproc/(double)npath*100 >= pratio ) {
            fprintf( stdout, "%d/%d: %.2f sec\n", nproc, int(npath), mytime()-lt0);
            pratio += 1;
        }
        
        if ( param.verbose ) {
            std::cout << "\nQuery:" << query_pid << "\tLength:" << it->first << "\n";
            std::cout << query << "\n";
        }
        
        size_t nkmers = path2aln_map[query_pid]->getKmerCount();
        //assert(nkmers>0);
        if ( nkmers == 0 ) {
            std::cout << "[Warning] zero kmers\n"; 
            std::cout << "Query:" << query_pid << "\tLength:" << it->first << "\n";
            std::cout << query << "\n";
            continue;
        }
        int min_filter = filter::minSameKmerCount( nkmers+param.kmer_size-1, param.filter_kmer, 1-param.filter_score );
        if ( min_filter < 1 ) min_filter = 1;
        if (param.verbose) std::cout << "Min k-mer:" << min_filter << "\n";

        double t0 = mytime();        
        PathFreqMap path_freq;
        path_index.findPaths( path_freq, query_kmers );
        std::multimap<size_t, PathId> freq_map = util::sortByValue<PathId, size_t>(path_freq);
        if ( param.verbose ) std::cout << "Scan paths time:" << mytime() -t0 << "\n";
                
        bool found = false;
        //PathId target = query_pid;
        size_t niter = 0;
        for ( std::multimap<size_t, PathId>::reverse_iterator it = freq_map.rbegin(); it != freq_map.rend(); ++it ) {
            if ( (int)it->first < min_filter ) break;

            PathId sbjct_pid = it->second;
            std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();
            assert(sbjct.size()>0);

            KmerArray sbjct_kmers = biostr::getKmers( sbjct, param.filter_kmer );
            //double tm0 = mytime();
            SemiGlobalAlign aln(sbjct, query, ANCHOR_CENTER, param.gap_ext, param.gap_open);
            AlignSummary summary = aln.getSummary();
            
            if ( param.verbose ) aln.printAlignment(std::cout);
            if ( param.verbose ) summary.print(std::cout);
            if ( param.verbose ) 
                std::cout << ++niter << "\tquery:" << query_pid << "\tsbjct:" << sbjct_pid 
                          << "\tlen1:" << query.size() << "\tlen2:" << sbjct.size() 
                          << "\t#kmer:" << it->first << "\tscore:" << summary.posrate << "\n";
            if ( summary.posrate >= param.merge_score ) {
                found = true; 
                merged_paths.insert( query_pid );
                path_bins[sbjct_pid].push_back( PathAlignPair(query_pid, summary) );
                break; // greedy
            }
        }
        if ( found ) {
            if ( param.verbose ) std::cout << "success:"  << mytime()-t0 << " sec\n";
        }

        if ( !found ) {
            path_bins.insert( std::pair<PathId, PathAlignPairList>( query_pid, PathAlignPairList() ) );
            path_index.add( query_kmers, query_pid );
            
            if ( param.verbose ) std::cout << "# update time:" << mytime()-t0 << "\n";
            if ( param.verbose ) std::cout << "fail:" <<  "\t# clusters:" << path_bins.size() << "\t# kmers:" << path_index.getSize()  << "\t time:" << mytime()-t0 << " sec\n";
        }
    }

    if ( param.verbose ) std::cout << "# Merged paths:" << merged_paths.size() << "\t" << mytime()-ts << " sec\n";

    ts = mytime();
    mergeClusters( path_bins, path2aln_map, bstrs, used_reads, strands, pairs, param );
    if ( param.verbose ) std::cout << "Paths merged:" << mytime()-ts << " sec\n";

    /* Remove merged paths from the path list */
    dropMergedPaths(merged_paths, path2aln_map);
    
    if ( param.verbose ) std::cout << "#Total paths:" << path2aln_map.size() << "\n";
}

void assem::clusterPairedPaths( PathToAlnMap &path2aln_map, 
                                InvertedIndex &iindex,
                                PathId *used_reads,
                                BitString *bstrs,
                                char *strands,
                                ReadId *pairs,
                                int nreads,
                                Param &param )
{
    double ts = mytime();
    
    setbuf(stdout, NULL); // no buffering
    if ( param.verbose ) std::cout << "\nClustering paths\n";
    if ( param.verbose ) std::cout << "#Total paths:" << path2aln_map.size() << "\n";

    
    /* Initialize Index */
    /* Avoid large hash miss */
    //size_t count = 1;
    KmerPathIndex path_index;
    double t0 = mytime();
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
        std::string query = (*it).second->getConsensusString();
        assert(query.size()>0);
        KmerArray kmers = biostr::getKmers( query, param.filter_kmer );
        path_index.init(kmers);
    }
    if ( param.verbose ) std::cout << "Path index init:" << mytime()-t0 << "\n";

    /* Keep track of merged paths */
    PathIdSet merged_paths; 

    /* Bins */
    PathBins path_bins;

    /* auxilary variables for progress log */
    double pratio = 1.0;
    int    nproc  = 0;
    double lt0    = mytime();
    size_t npath  = path2aln_map.size();

    /* Make a path length map */
    PathLengthMap plen_map = getPathLengths( path2aln_map );
    for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        PathId query_pid = it->second;

        std::string query = path2aln_map[query_pid]->getConsensusString();
        assert(query.size() > 0);
        KmerArray query_kmers = biostr::getKmers( query, param.filter_kmer );

        /* Print progress status at every percent of path proceeded */
        nproc++;
        if ( nproc/(double)npath*100 >= pratio ) {
            fprintf( stdout, "%d/%d: %.2f sec\n", nproc, int(npath), mytime()-lt0);
            pratio += 1;
        }
        
        if ( param.verbose ) {
            std::cout << "\nQuery:" << query_pid << "\tLength:" << it->first << "\n";
            std::cout << query << "\n";
        }
        
        size_t nkmers = path2aln_map[query_pid]->getKmerCount();
        //assert(nkmers>0);
        if ( nkmers == 0 ) {
            std::cout << "[Warning] zero kmers\n"; 
            std::cout << "\nQuery:" << query_pid << "\tLength:" << it->first << "\n";
            std::cout << query << "\n";
            continue;
        }
        int min_filter = filter::minSameKmerCount( nkmers+param.kmer_size-1, param.filter_kmer, 1-param.filter_score );
        if ( min_filter < 1 ) min_filter = 1;
        if (param.verbose) std::cout << "Min k-mer:" << min_filter << "\n";

        double t0 = mytime();        
        PathFreqMap path_freq;
        path_index.findPaths( path_freq, query_kmers );

        PathIdArray pids = findPairedPaths( query_pid, used_reads, pairs, path2aln_map[query_pid]->getReads(), path2aln_map[query_pid]->getReadCount(), param );
        std::tr1::unordered_map<PathId, bool> pair_paths;
        for ( size_t i = 0; i < pids.size(); i++ ) {
            pair_paths.insert( std::pair<PathId, bool>( pids[i], true ) );
        }

        for ( PathFreqMap::iterator pt = path_freq.begin(); pt != path_freq.end(); ) {
            if ( pair_paths.find( pt->first ) != pair_paths.end() ) ++pt;
            else path_freq.erase(pt++);
        }

        std::multimap<size_t, PathId> freq_map = util::sortByValue<PathId, size_t>(path_freq);
        if ( param.verbose ) std::cout << "Scan paths time:" << mytime() -t0 << "\n";

        bool found = false;

        //PathId target = query_pid;
        size_t niter = 0;
        for ( std::multimap<size_t, PathId>::reverse_iterator it = freq_map.rbegin(); it != freq_map.rend(); ++it ) {
            if ( (int)it->first < min_filter ) break;

            PathId sbjct_pid = it->second;
            std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();
            assert(sbjct.size()>0);

            KmerArray sbjct_kmers = biostr::getKmers( sbjct, param.filter_kmer );
            //double tm0 = mytime();
            SemiGlobalAlign aln(sbjct, query, ANCHOR_CENTER, param.gap_ext, param.gap_open);
            AlignSummary summary = aln.getSummary();
            
            if ( param.verbose ) 
                std::cout << ++niter << "\tquery:" << query_pid << "\tsbjct:" << sbjct_pid 
                          << "\tlen1:" << query.size() << "\tlen2:" << sbjct.size() 
                          << "\t#kmer:" << it->first << "\tscore:" << summary.posrate << "\n";
            if ( summary.posrate >= param.merge_score ) {
                found = true; 
                merged_paths.insert( query_pid );
                path_bins[sbjct_pid].push_back( PathAlignPair(query_pid, summary) );
                break; // greedy
            }
        }

        if ( found ) {
            if ( param.verbose ) std::cout << "success:"  << mytime()-t0 << " sec\n";
        }

        if ( !found ) {
            path_bins.insert( std::pair<PathId, PathAlignPairList>( query_pid, PathAlignPairList() ) );
            path_index.add( query_kmers, query_pid );
            
            if ( param.verbose ) std::cout << "# update time:" << mytime()-t0 << "\n";
            if ( param.verbose ) std::cout << "fail:" <<  "\t# clusters:" << path_bins.size() << "\t# kmers:" << path_index.getSize()  << "\t time:" << mytime()-t0 << " sec\n";
        }
    }

    if ( param.verbose ) std::cout << "# Merged paths:" << merged_paths.size() << "\t" << mytime()-ts << " sec\n";

    ts = mytime();
    mergeClusters( path_bins, path2aln_map, bstrs, used_reads, strands, pairs, param );
    if ( param.verbose ) std::cout << "Paths merged:" << mytime()-ts << " sec\n";

    /* Remove merged paths from the path list */
    dropMergedPaths(merged_paths, path2aln_map);
    
    if ( param.verbose ) std::cout << "#Total paths:" << path2aln_map.size() << "\n";
}

void assem::recruitPaths( PathIdSet &centroids,
                          PathToAlnMap &path2aln_map, 
                          InvertedIndex &iindex,
                          PathId *used_reads,
                          BitString *bstrs,
                          char *strands,
                          ReadId *pairs,
                          int nreads,
                          Param &param )
{
    double ts = mytime();

    setbuf(stdout, NULL); // no buffering
    if ( param.verbose ) std::cout << "\nPath recruitment\n";
    if ( param.verbose ) std::cout << "#Total paths:" << path2aln_map.size() << "\n";

    
    /* Initialize Index */
    /* Avoid large hash miss */
    //size_t count = 1;
    KmerPathIndex path_index;
    double t0 = mytime();
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
        std::string query = (*it).second->getConsensusString();
        assert(query.size()>0);
        KmerArray kmers = biostr::getKmers( query, param.filter_kmer );
        path_index.init(kmers);
    }
    if ( param.verbose ) std::cout << "Path index init:" << mytime()-t0 << "\n";

    size_t max_center = 0;
    for ( PathIdSet::iterator it = centroids.begin(); it != centroids.end(); ++it ) {
        std::string query = path2aln_map[*it]->getConsensusString();
        KmerArray kmers = biostr::getKmers( query, param.filter_kmer );
        PathId pid = *it;
        path_index.add(kmers, pid);
        if ( query.size() > max_center ) max_center = query.size();
    }
    if ( param.verbose ) std::cout << "Path index construction:" << mytime()-t0 << "\n";


    /* Keep track of merged paths */
    PathIdSet merged_paths; 

    /* Bins */
    PathBins path_bins;

    /* auxilary variables for progress log */
    double pratio = 1.0;
    int    nproc  = 0;
    double lt0    = mytime();
    size_t npath  = path2aln_map.size();

    /* Make a path length map */
    PathLengthMap plen_map = getPathLengths( path2aln_map );
    for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        PathId query_pid = it->second;

        std::string query = path2aln_map[query_pid]->getConsensusString();
        assert(query.size() > 0);
        if ( query.size() >= max_center ) continue;
        
        KmerArray query_kmers = biostr::getKmers( query, param.filter_kmer );

        /* Print progress status at every percent of path proceeded */
        nproc++;
        if ( nproc/(double)npath*100 >= pratio ) {
            fprintf( stdout, "%d/%d: %.2f sec\n", nproc, int(npath), mytime()-lt0);
            pratio += 1;
        }
        
        if ( param.verbose ) {
            std::cout << "\nQuery:" << query_pid << "\tLength:" << it->first << "\n";
            std::cout << query << "\n";
        }
        
        size_t nkmers = path2aln_map[query_pid]->getKmerCount();
        assert(nkmers>0);
        int min_filter = filter::minSameKmerCount( nkmers+param.kmer_size-1, param.filter_kmer, 1-param.filter_score );
        if ( min_filter < 1 ) min_filter = 1;
        if (param.verbose) std::cout << "Min k-mer:" << min_filter << "\n";

        double t0 = mytime();        
        PathFreqMap path_freq;
        path_index.findPaths( path_freq, query_kmers );
        std::multimap<size_t, PathId> freq_map = util::sortByValue<PathId, size_t>(path_freq);
        if ( param.verbose ) std::cout << "Scan paths time:" << mytime() -t0 << "\n";
                
        bool found = false;
        //PathId target = query_pid;
        size_t niter = 0;
        for ( std::multimap<size_t, PathId>::reverse_iterator it = freq_map.rbegin(); it != freq_map.rend(); ++it ) {
            if ( (int)it->first < min_filter ) break;

            PathId sbjct_pid = it->second;
            if ( query_pid == sbjct_pid ) continue;
            
            std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();
            assert(sbjct.size()>0);

            KmerArray sbjct_kmers = biostr::getKmers( sbjct, param.filter_kmer );
            //double tm0 = mytime();
            SemiGlobalAlign aln(sbjct, query, ANCHOR_CENTER, param.gap_ext, param.gap_open);
            AlignSummary summary = aln.getSummary();
            
            if ( param.verbose ) 
                std::cout << ++niter << "\tquery:" << query_pid << "\tsbjct:" << sbjct_pid 
                          << "\tlen1:" << query.size() << "\tlen2:" << sbjct.size() 
                          << "\t#kmer:" << it->first << "\tscore:" << summary.posrate << "\n";
            if ( summary.posrate >= param.merge_score ) {
                found = true; 
                merged_paths.insert( query_pid );
                path_bins[sbjct_pid].push_back( PathAlignPair(query_pid, summary) );
                break; // greedy
            }
        }
        if ( found ) {
            if ( param.verbose ) std::cout << "success:"  << mytime()-t0 << " sec\n";
        }

        if ( !found ) {
//             path_bins.insert( std::pair<PathId, PathAlignPairList>( query_pid, PathAlignPairList() ) );
//             path_index.add( query_kmers, query_pid );
            
//             if ( param.verbose ) std::cout << "# update time:" << mytime()-t0 << "\n";
            if ( param.verbose ) std::cout << "fail:" <<  "\t# clusters:" << path_bins.size() << "\t# kmers:" << path_index.getSize()  << "\t time:" << mytime()-t0 << " sec\n";
        }
    }

    if ( param.verbose ) std::cout << "\n# Paths to merge:" << merged_paths.size() << "\t" << mytime()-ts << " sec\n";

    ts = mytime();
    mergeClusters( path_bins, path2aln_map, bstrs, used_reads, strands, pairs, param );
    if ( param.verbose ) std::cout << "Paths merged:" << mytime()-ts << " sec\n";

    /* Remove merged paths from the path list */
    dropMergedPaths(merged_paths, path2aln_map);
    
    if ( param.verbose ) std::cout << "#Total paths:" << path2aln_map.size() << "\n";
}

void assem::findLatchablePath( size_t min_filter,
                               PathId sbjct_pid,
                               KmerArray sbjct_kmers,
                               KmerPathIndex &path_index,
                               PathIdSet &skip_pids,
                               PathToAlnMap &path2aln_map, 
                               PathAlignPairList &palign_list,
                               ReadId *pairs,
                               char *strands,
                               int nreads, 
                               int direction,
                               Param &param, 
                               size_t filter_kmer )
{
    PathFreqMap path_freq;
    double t0 = mytime();
    path_index.findPaths( path_freq, sbjct_kmers );
    std::multimap<size_t, PathId> freq_map = util::sortByValue<PathId, size_t>(path_freq);
    if ( param.verbose ) std::cout << "\nSimilar path:" << path_freq.size() << "\tscan time:" << mytime() -t0 << "\n";
    
    double aln_time = 0;

    std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();

    int anchor = (direction == LEFT) ? ANCHOR_LEFT : ANCHOR_RIGHT;

    PathId max_path = NOT_PATH;
    int max_pair = -1;
    int max_same =-1;
    AlignSummary max_summ;
    for ( std::multimap<size_t, PathId>::reverse_iterator it = freq_map.rbegin(); it != freq_map.rend(); ++it ) {
        if ( it->first < min_filter ) break;
        PathId query_pid = it->second;

        if ( skip_pids.find( query_pid ) != skip_pids.end() ) continue;
        
        std::string query = path2aln_map[query_pid]->getConsensusString();

        AlignSummary summary;
        std::set<IntPair> bad_range;
        double tas = mytime();
        bool found = alignEnd(summary, sbjct, query, query_pid, sbjct_pid, bad_range, path2aln_map, param, anchor ); 
        aln_time += mytime()-tas;

        if ( param.verbose ) std::cout << "Found? " << found << "\n";
        if ( !found ) continue;
        /* pair-end information */

        ReadIdArray sbjpreads; // reads in subjct that paired with reads in query
        int orientation = NONSENSE;
        if ( param.pair_flag ) {
            sbjpreads = getCommonReads( query_pid, sbjct_pid, path2aln_map, pairs, nreads, param );
            orientation = getDirection(sbjpreads, pairs, strands, param );

            if ( param.verbose ) std::cout << "common:" << sbjpreads.size() << "\tdirection:" << orientation << "\n";
            /* Check orientation conflict */
            if ( ( orientation == POSITIVE && anchor == ANCHOR_LEFT  ) || 
                 ( orientation == NEGATIVE && anchor == ANCHOR_RIGHT ) ) { 
                if ( param.verbose ) std::cout << "\tDirectionality conflicts\n"; 
                continue;
            }
        }

        IntPair overlap = findOverlap( query_pid, sbjct_pid, path2aln_map, direction, param );
        if ( param.verbose && overlap.second ) std::cout << "overlap:" << overlap.second << "\n";
        if ( param.pair_flag ) {
            if ( (int)sbjpreads.size() > max_pair ) {
                max_path = query_pid;
                max_pair = sbjpreads.size();
                max_summ = summary;
            }
        }
        else {
            if ( overlap.second > max_same ) {
                max_path = query_pid;
                max_same = overlap.second;
                max_summ = summary;
            }
        }
        if ( found ) break;
    }

    if ( param.verbose ) std::cout << "Alignment time:" << aln_time << "\n";

    if ( max_path == NOT_PATH ) {
        if ( param.verbose ) std::cout << "No latchable path found\n";
        return;
    }

    palign_list.push_back( PathAlignPair( max_path, max_summ ) );
    skip_pids.insert( max_path );
    std::string nsbj = path2aln_map[max_path]->getConsensusString();
    sbjct_kmers = biostr::getKmers( nsbj, filter_kmer );
    sbjct_kmers = (direction == LEFT) ? KmerArray(sbjct_kmers.begin(), sbjct_kmers.begin()+sbjct_kmers.size()/2) : KmerArray(sbjct_kmers.begin()+sbjct_kmers.size()/2, sbjct_kmers.end());
    findLatchablePath( min_filter, max_path, sbjct_kmers, path_index, skip_pids, path2aln_map, palign_list, pairs, strands, nreads, direction, param, filter_kmer );
}

void assem::connectOverlappingPairedPaths( PathToAlnMap &path2aln_map, 
                                           BitString *bstrs,
                                           char *strands,
                                           ReadId *pairs, 
                                           PathId *used_reads,
                                           int nreads, 
                                           InvertedIndex &iindex,
                                           Param &param )
{
    double t0 = mytime();

    setbuf(stdout, NULL); // no buffering
    //if ( param.verbose ) 
    std::cout << "\nConnecting overlapping paired-end paths ...\n";

    PathIdSet merged_paths; // Keep track of merged paths
    PathIdSet successPaths; // Keep track of new extended pathso

    Map<PathIdPair, bool>::Type scanned_pairs;

    //-------------------------
    // k-mer to path id mapping
    //-------------------------
    KmerToPathMap pathid_map;
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) 
        addPathIds( it->first, path2aln_map[it->first]->getKmers(), path2aln_map[it->first]->getKmerCount(), pathid_map );

    /* auxilary variables for progress log */
    double pratio = 1.0;
    int    nproc  = 0;
    double lt0    = mytime();
    size_t npath  = path2aln_map.size();

    /* Make a path length map */
    PathLengthMap plen_map = getPathLengths( path2aln_map );
    for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        /* Print progress status at every percent of path proceeded */
        nproc++;
        if ( nproc/(double)npath*100 >= pratio ) {
            fprintf( stdout, "%d/%d: %.2f sec\n", nproc, int(npath), mytime()-lt0);
            pratio += 1;
        }

        PathId sbjct_pid = it->second;
        std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();
        //if ( it->first < minlen ) break;
        
        if ( merged_paths.find( sbjct_pid ) != merged_paths.end() ) continue;

        if ( param.verbose ) {
            std::cout << "\nCenter:" << sbjct_pid << "\tLength:" << it->first << "\n";
            std::cout << sbjct << "\n";
        }

        /* Recursive latching */
        while(true) {
            
            PathIdArray pids = findPairedPaths( sbjct_pid, used_reads, pairs, path2aln_map[sbjct_pid]->getReads(), path2aln_map[sbjct_pid]->getReadCount(), param );
            if ( pids.size() == 0 ) break;

            PathIdList reduced;
            for ( PathIdArray::iterator pt = pids.begin(); pt != pids.end(); ++pt ) {
                PathIdPair ppair = sbjct_pid > *pt ? PathIdPair(*pt, sbjct_pid) : PathIdPair(sbjct_pid, *pt);
                if ( scanned_pairs.find( ppair ) == scanned_pairs.end() && 
                     merged_paths.find( *pt ) == merged_paths.end() )
                    reduced.push_back( *pt );
            }
            
            pids = PathIdArray( reduced.begin(), reduced.end() );
            if ( pids.size() == 0 ) break;


            std::list<PairedPath> ppaths = makePairedPaths( pids, sbjct_pid, pairs, used_reads, strands, path2aln_map, param );
            if ( ppaths.size() == 0 ) break;

            std::multimap<int, PairedPath> orders;
            orderPathsByOverlap( orders, ppaths, param );
            
            bool success = false;
            for ( std::multimap<int, PairedPath>::reverse_iterator ot = orders.rbegin(); ot != orders.rend(); ++ot ) {
                PathIdPair ppair = sbjct_pid > ot->second.pid ? PathIdPair(ot->second.pid, sbjct_pid) : PathIdPair(sbjct_pid, ot->second.pid);
//                 if ( scanned_pairs.find( ppair ) != scanned_pairs.end() ) {
//                     if ( param.verbose ) std::cout << "Scanned early - skip\n";
//                     continue;
//                 }

                PairedPath best_pair = ot->second;
                if (param.verbose) {
                    std::cout << "Max Query:" << best_pair.pid << "\t" << path2aln_map[best_pair.pid]->getConsensusString() << "\n";
                    std::cout << "Max Overlap:" << best_pair.overlap << "\n";
                    std::cout << "Max Support:" << best_pair.match_reads.size() << "\n";
                    std::cout << "Max Links:" << best_pair.latch_reads.size() << "\n";
                    std::cout << "Max Direction:" << best_pair.direction << "\n";
                }
                
                if ( best_pair.pid == NOT_PATH ) {
                    if ( param.verbose ) std::cout << "Not a path\n"; 
                    scanned_pairs.insert( std::pair<PathIdPair, bool>( ppair, true) );
                    continue;
                }
                if ( best_pair.overlap < param.pairend_overlap ) {
                    if ( param.verbose ) std::cout << "Weak overlap\n"; 
                    scanned_pairs.insert( std::pair<PathIdPair, bool>( ppair, true) );
                    continue;
                }
                
                if ( (int)best_pair.match_reads.size() < param.overlap_support ) {
                    if ( param.verbose ) std::cout << "Weak support\n"; 
                    scanned_pairs.insert( std::pair<PathIdPair, bool>( ppair, true) );
                    continue;
                }

                success = linkPairedPath( sbjct_pid, best_pair, pathid_map, path2aln_map, bstrs, strands, pairs, used_reads, nreads, merged_paths, successPaths, iindex, param.kmer_size, param );
                if ( success ) {
                    scanned_pairs.insert( std::pair<PathIdPair, bool>( ppair, true) );
                    if ( param.verbose ) std::cout << "Search recursive latch for path:" << sbjct_pid << "\n";
                    break;
                }
            }
            if ( !success ) break;
        }
    }

    dropMergedPaths(merged_paths, path2aln_map);

    std::cout << "# Merged paths:" << merged_paths.size() << "\n";
    std::cout << "Path count:" << path2aln_map.size() << "\n";
    std::cout << "Latching overlapping paths:" << mytime()-t0 << " sec\n";
}

void assem::connectLongOverlappingPaths(PathToAlnMap &path2aln_map, 
                                        BitString *bstrs,
                                        char *strands,
                                        ReadId *pairs, 
                                        PathId *used_reads,
                                        int nreads, 
                                        InvertedIndex &iindex,
                                        Param &param )
{
    double t0 = mytime();

    setbuf(stdout, NULL); // no buffering
    //if ( param.verbose ) 
    std::cout << "\nConnecting long overlapping paths ...\n";

    if ( param.verbose ) std::cout << "#Total paths:" << path2aln_map.size() << "\n";

    PathIdSet merged_paths; // Keep track of merged paths
    PathIdSet success_paths; // Keep track of new extended pathso

    Map<PathIdPair, bool>::Type scanned_pairs;

    //-------------------------
    // k-mer to path id mapping
    //-------------------------
    KmerToPathMap pathid_map;
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) 
        addPathIds( it->first, path2aln_map[it->first]->getKmers(), path2aln_map[it->first]->getKmerCount(), pathid_map );

    double pratio = 1.0;
    int    nproc  = 0;
    double lt0    = mytime();
    size_t npath  = path2aln_map.size();

    /* make a path length map */
    PathLengthMap plen_map = getPathLengths( path2aln_map );
    for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        nproc++;
        if ( nproc/(double)npath*100 >= pratio ) {
            fprintf( stdout, "%d/%d: %.2f sec\n", nproc, int(npath), mytime()-lt0);
            pratio += 1;
        }

        PathId sbjct_pid = it->second;
        std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();
        if ( merged_paths.find( sbjct_pid ) != merged_paths.end() ) continue;

        if ( param.verbose ) {
            std::cout << "\ncenter:" << sbjct_pid << "\tlength:" << it->first << "\n";
            std::cout << sbjct << "\n";
        }

        size_t srch_nkmer = 2 * param.latch_length - param.kmer_size + 1;
        if ( path2aln_map[sbjct_pid]->getKmerCount() < srch_nkmer ) 
            srch_nkmer = path2aln_map[sbjct_pid]->getKmerCount();

        //if ( path2aln_map[sbjct_pid]->getKmerCount() < srch_nkmer ) continue;
        
        /* recursive latching left */
        while(true) {
//             PathIdSet self;
//             self.insert(sbjct_pid);
            KmerId *kmers = path2aln_map[sbjct_pid]->getKmers();
            // size_t  nkmer = path2aln_map[sbjct_pid]->getKmerCount();
            
            //PosPathPairList pos_paths = searchSimilarPaths(self, kmers, nkmer/2, path2aln_map, pathid_map, merged_paths, 1, param);
            PathIdSet candidates = searchOverlapCandidates(sbjct_pid, kmers, srch_nkmer, path2aln_map, pathid_map, merged_paths, 1, param, LEFT);
            //PathIdSet candidates = searchOverlapCandidates(sbjct_pid, kmers, nkmer, srch_nkmer, path2aln_map, pathid_map, merged_paths, 1, param, LEFT);
            if ( param.verbose ) std::cout << "# similar paths (left):" << candidates.size() << "\n";

            PathIdList reduced;
            for ( PathIdSet::iterator pt = candidates.begin(); pt != candidates.end(); ++pt ) {
                PathIdPair ppair = sbjct_pid > *pt ? PathIdPair(*pt, sbjct_pid) : PathIdPair(sbjct_pid, *pt);
                if ( //scanned_pairs.find( ppair ) == scanned_pairs.end() && 
                     merged_paths.find( *pt ) == merged_paths.end() )
                    reduced.push_back( *pt );
            }
            
            candidates = PathIdSet( reduced.begin(), reduced.end() );
            if ( candidates.size() == 0 ) break;

            bool success = false;
            for ( PathIdSet::iterator jt = candidates.begin(); jt != candidates.end(); ++jt ) {
                PathId query_pid = *jt;
                if ( param.verbose ) std::cout << "query:" << *jt << "\t" << path2aln_map[*jt]->getConsensusString() << "\n";
                IntPair found = checkLongOverlap(query_pid, sbjct_pid, path2aln_map, param.kmer_size, LEFT, param );
                if ( found.first != LEFT || found.second < param.latch_length ) {
                    PathIdPair ppair = sbjct_pid > *jt ? PathIdPair(*jt, sbjct_pid) : PathIdPair(sbjct_pid, *jt);
                    scanned_pairs.insert( std::pair<PathIdPair, bool>( ppair, true) );
                    continue;
                }
                ReadIdArray latch_reads;
                std::vector<int> read_inits;
                success = connectPathPair(bstrs, strands, pairs, path2aln_map, pathid_map, merged_paths, success_paths, sbjct_pid, query_pid, used_reads, latch_reads, read_inits, found.second, NEGATIVE, iindex, param );
                if ( success ) break;
            }
            if ( !success ) break;
        }
        
        /* recursive latching right */
        while(true) {
//             PathIdSet self;
//             self.insert(sbjct_pid);
            KmerId *kmers = path2aln_map[sbjct_pid]->getKmers();
            size_t  nkmer = path2aln_map[sbjct_pid]->getKmerCount();

            //PosPathPairList candidates = searchSimilarPaths(self, &kmers[(int)ceil(nkmer/2)], nkmer/2, path2aln_map, pathid_map, merged_paths, 1, param);
            PathIdSet candidates = searchOverlapCandidates(sbjct_pid, &kmers[nkmer-srch_nkmer], srch_nkmer, path2aln_map, pathid_map, merged_paths, 1, param, RIGHT);
            //PathIdSet candidates = searchOverlapCandidates(sbjct_pid, kmers, nkmer, srch_nkmer, path2aln_map, pathid_map, merged_paths, 1, param, RIGHT);
            if ( param.verbose ) std::cout << "# similar paths (right):" << candidates.size() << "\n";
            
            PathIdList reduced;
            for ( PathIdSet::iterator pt = candidates.begin(); pt != candidates.end(); ++pt ) {
                PathIdPair ppair = sbjct_pid > *pt ? PathIdPair(*pt, sbjct_pid) : PathIdPair(sbjct_pid, *pt);
                if ( //scanned_pairs.find( ppair ) == scanned_pairs.end() && 
                     merged_paths.find( *pt ) == merged_paths.end() )
                    reduced.push_back( *pt );
            }
            
            candidates = PathIdSet( reduced.begin(), reduced.end() );
            if ( candidates.size() == 0 ) break;

            bool success = false;
            for ( PathIdSet::iterator jt = candidates.begin(); jt != candidates.end(); ++jt ) {
                PathId query_pid = *jt;
                //if ( merged_paths.find( query_pid ) != merged_paths.end() ) continue;
                if ( param.verbose ) std::cout << "query:" << *jt << "\t" << path2aln_map[*jt]->getConsensusString() << "\n";

                IntPair found = checkLongOverlap(query_pid, sbjct_pid, path2aln_map, param.kmer_size, RIGHT, param );
                if ( found.first != RIGHT || found.second < param.latch_length ) {
                    PathIdPair ppair = sbjct_pid > *jt ? PathIdPair(*jt, sbjct_pid) : PathIdPair(sbjct_pid, *jt);
                    scanned_pairs.insert( std::pair<PathIdPair, bool>( ppair, true) );
                    continue;
                }
                ReadIdArray latch_reads;
                std::vector<int> read_inits;
                success = connectPathPair(bstrs, strands, pairs, path2aln_map, pathid_map, merged_paths, success_paths, sbjct_pid, query_pid, used_reads, latch_reads, read_inits, found.second, POSITIVE, iindex, param );
                if ( success ) break;
            }
            if ( !success ) break;
        }
        dropPathIds( sbjct_pid, path2aln_map[sbjct_pid]->getKmers(), path2aln_map[sbjct_pid]->getKmerCount(), pathid_map );
    }

    dropMergedPaths(merged_paths, path2aln_map);

    std::cout << "# Merged paths:" << merged_paths.size() << "\n";
    std::cout << "Path count:" << path2aln_map.size() << "\n";
    std::cout << "Latching long overlapping paths:" << mytime()-t0 << " sec\n";
}

void assem::connectPaths(PathToAlnMap &path2aln_map, 
                         BitString *bstrs,
                         char *strands,
                         ReadId *pairs, 
                         PathId *used_reads,
                         int nreads, 
                         InvertedIndex &iindex,
                         Param &param )
{
    if ( !param.extend_flag ) return;
    double t0 = mytime();

    //std::cout << "\nStage 4:\nConnecting long/short ovelapping paths ....\n";
    std::cout << "\nExtending paths ....\n";
    if ( param.verbose ) std::cout << "#Total paths:" << path2aln_map.size() << "\n";

    if ( param.pair_flag ) {
        connectPairedPathsByBridgingReads( used_reads, path2aln_map, bstrs, strands, pairs, iindex, param );
        if ( param.extend_read_flag ) {
            connectPairedReadsToPath( used_reads, strands, pairs, bstrs, path2aln_map, param );
        }
        connectOverlappingPairedPaths( path2aln_map, bstrs, strands, pairs, used_reads, nreads, iindex, param );
    }

    connectLongOverlappingPaths(path2aln_map, bstrs, strands, pairs, used_reads, nreads, iindex, param );

    //if ( successpaths.size() ) {
    //double ts = mytime();
    //if ( param.verbose ) 
    //std::cout << "\nFinal path clustering\n";
    std::cout << "\nClustering paths\n";
    clusterPaths(path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param);
    //std::cout << "Final path clustering:" << mytime()-ts << " sec\n";
    // }
    std::cout << "Connecting paths:" << mytime()-t0 << " sec\n";
}

// void assem::tunePaths(PathToAlnMap &path2aln_map, 
//                          BitString *bstrs,
//                          char *strands,
//                          ReadId *pairs, 
//                          PathId *used_reads,
//                          int nreads, 
//                          InvertedIndex &iindex,
//                          Param &param )
// {
//     if ( !param.tune_flag ) return;
//     double t0 = mytime();
//     connectOverlappingPairedPaths( path2aln_map, bstrs, strands, pairs, used_reads, nreads, iindex, param );
//     clusterPaths(path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param);
//     std::cout << "\nTuning paths:" << mytime()-t0 << " sec\n";
// }

void assem::connectPairedReadsToPath( PathId *used_reads, 
                                      char *strands, 
                                      ReadId *pairs, 
                                      BitString *bstrs, 
                                      PathToAlnMap &path2aln_map, 
                                      Param &param )
{
    double t0 = mytime();

    setbuf(stdout, NULL); // no buffering
    //if ( param.verbose ) 
        std::cout << "\nAttaching reads to a path end\n";
  
    size_t count = 0;

    KmerToPathMap pathid_map;
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) 
        addPathIds( it->first, path2aln_map[it->first]->getKmers(), path2aln_map[it->first]->getKmerCount(), pathid_map );

    /* auxilary variables for progress log */
    double pratio = 1.0;
    int    nproc  = 0;
    double lt0    = mytime();
    size_t npath  = path2aln_map.size();

    PathLengthMap plen_map = getPathLengths( path2aln_map );
    for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        /* Print progress status at every percent of path proceeded */
        nproc++;
        if ( nproc/(double)npath*100 >= pratio ) {
            fprintf( stdout, "%d/%d: %.2f sec\n", nproc, int(npath), mytime()-lt0);
            pratio += 1;
        }

        PathId sbjct_pid = it->second;
        std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();
        if ( param.verbose) 
            std::cout << "\nSbjct:" << sbjct_pid << "\t" << sbjct << "\n";

        ReadIdArray pos_reads = findPairedReads( sbjct_pid, used_reads, pairs, strands, path2aln_map[sbjct_pid]->getReads(), path2aln_map[sbjct_pid]->getReadCount(), param, '+' );
        ReadIdArray neg_reads = findPairedReads( sbjct_pid, used_reads, pairs, strands, path2aln_map[sbjct_pid]->getReads(), path2aln_map[sbjct_pid]->getReadCount(), param, '-' );

        if ( param.verbose )
            std::cout << "+reads:" << pos_reads.size() << "\t-reads:" << neg_reads.size() << "\n";

        bool lflag = false;
        if ( (int)pos_reads.size() >= param.latch_support ) {
            lflag = attachReadsToPathEnd( pos_reads, sbjct_pid, path2aln_map, pathid_map, bstrs, used_reads, RIGHT, param ); // right path, prepend reads to left end
        }
        bool rflag = false;
        if ( (int)neg_reads.size() >= param.latch_support ) {
            if ( sbjct[sbjct.size()-1] == '*' ) continue;
            rflag = attachReadsToPathEnd( neg_reads, sbjct_pid, path2aln_map, pathid_map, bstrs, used_reads, LEFT, param); // left path, append reads to right end
        }

        if ( param.verbose ) std::cout << "flags:" << lflag << "," << rflag << "\n";
        if ( !lflag && !rflag ) continue;

        count++;

        SpaPath *old_path = new SpaPath(*path2aln_map[sbjct_pid]);
        SpaPath *new_path = path2aln_map[sbjct_pid];
        int error = 0;
		MSA msa(*new_path, sbjct_pid, used_reads, bstrs, strands, pairs, PILE|FLAT|TRIM|CODON, error, param);
        
        if ( error ) {
            std::cout << "Error occurred\n";
            delete new_path;
            path2aln_map[sbjct_pid] = old_path;
            for ( size_t i = 0; i < pos_reads.size(); i++ ) 
                if ( used_reads[pos_reads[i]] == sbjct_pid ) used_reads[pos_reads[i]] = NOT_PATH;
            for ( size_t i = 0; i < neg_reads.size(); i++ ) 
                if ( used_reads[neg_reads[i]] == sbjct_pid ) used_reads[neg_reads[i]] = NOT_PATH;
            continue;
        } else delete old_path;

        if ( param.verbose ) {
            std::cout << "\nAlignment\n";
            path2aln_map[sbjct_pid]->printAlignment(std::cout, 100, bstrs);
        }

    }
    //if ( param.verbose ) 
    std::cout << "Extendend paths:" << count << "\n";
    //if ( param.verbose ) 
    std::cout << "Attaching reads to path end:" << mytime()-t0 << " sec\n";
}

bool assem::attachReadsToPathEnd( ReadIdArray &reads,
                                  PathId pid,
                                  PathToAlnMap &path2aln_map,
                                  KmerToPathMap &pathid_map,
                                  BitString *bstrs,
                                  PathId *used_reads, 
                                  int direction, // path direction
                                  Param &param )
{
    ReadIdArray nreads = getLatchableReads( reads, bstrs, pid, path2aln_map, pathid_map, param );

    PathReadMap latch_map;
    makePathReadMap( nreads, pid, bstrs, latch_map, path2aln_map, param, direction );

    std::vector<int> starts;
    ReadIdArray good_reads;
    for ( PathReadMap::iterator it = latch_map.begin(); it != latch_map.end(); ++it ) {
        if ( it->first != pid ) continue;
        
        for ( ReadPosList::iterator jt = it->second.begin(); jt != it->second.end(); ++jt ) {
            int start = (direction==RIGHT)? -1 * jt->spos : jt->rpos;
            starts.push_back( start );
            good_reads.push_back( jt->rid );

            if ( param.verbose ) std::cout << "ReadId:" << jt->rid << "\tStart:" << start << "\tDirection:" << direction << "\n";
        }
    }

    if ( (int)good_reads.size() < param.latch_support ) {
        if ( param.verbose ) std::cout << "Weak one side read latch:" << good_reads.size() << "\n";
        return false;
    }

    path2aln_map[pid]->addReads( &good_reads[0], &starts[0], good_reads.size() );
    
    return true;
}

ReadIdArray assem::getLatchableReads( ReadIdArray &reads, 
                                      BitString *bstrs,
                                      PathId pid,
                                      PathToAlnMap &path2aln_map,
                                      KmerToPathMap &pathid_map,
                                      //int direction,
                                      Param &param )
{
    PathIdSet merged_paths;
    ReadIdList good_reads;
    //KmerId *rkmers = path2aln_map[pid]->getKmers();
    //size_t rnkmers = path2aln_map[pid]->getKmerCount();

    for ( size_t i = 0; i < reads.size(); i++ ) {
        std::string str = bstrs[reads[i]].toString();

        KmerArray kmers = biostr::getKmers( str, param.kmer_size );
        PosPathPairList pos_paths = searchSimilarPaths(PathIdSet(), &kmers[0], kmers.size(), path2aln_map, pathid_map, merged_paths, 1, true, param);

        for ( PosPathPairList::iterator jt = pos_paths.begin(); jt != pos_paths.end(); ++jt ) {
            if ( jt->second == pid ) {
                good_reads.push_back(reads[i]);
                break;
            }
        }
    }
    return ReadIdArray( good_reads.begin(), good_reads.end() );
}                               

ReadIdArray assem::findPairedReads( PathId &pid,
                                    PathId *used_reads, 
                                    ReadId *pairs, 
                                    char *strands, 
                                    ReadId *reads,
                                    size_t nreads,
                                    Param &param,
                                    char strand )
{
    ReadIdList preads;
    for ( size_t i = 0; i < nreads; i++ ) {
        ReadId pread = pairs[reads[i]];
        if ( pread == NOT_PAIR ) continue;
        if ( used_reads[pread] != NOT_PATH ) continue;
        if ( strands[pread] != strand ) continue;
        if ( strands[reads[i]] == strand ) continue; // do not allow same strand

        preads.push_back(pread);
    }
    return ReadIdArray( preads.begin(), preads.end() );
}
                             
// pathidset assem::connectpaths( pathtoalnmap &path2aln_map, 
//                                invertedindex &iindex,
//                                pathid *used_reads,
//                                bitstring *bstrs,
//                                char *strands,
//                                readid *pairs,
//                                int nreads,
//                                param &param )
// {
//     pathidlist success_paths;

//     double ts = mytime();
//     if ( param.verbose ) std::cout << "\nlatching paths\n";
//     if ( param.verbose ) std::cout << "#total paths:" << path2aln_map.size() << "\n";
    

//     double scan_time = 0;
//     double alig_time = 0;
//     double merg_time = 0;

//     size_t filter_kmer = param.kmer_size > param.filter_kmer ? param.kmer_size : param.filter_kmer ;

//     /* initialize index */
//     /* avoid large hash miss */
//     size_t count = 1;
//     size_t minlen = 50;
    
//     kmerpathindex path_index;
//     double t0 = mytime();
//     for ( pathtoalnmap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
//         std::string query = (*it).second->getconsensusstring();
//         assert(query.size()>0);
//         kmerarray kmers = biostr::getkmers( query, filter_kmer );
//         path_index.init(kmers);
//     }
//     if ( param.verbose ) std::cout << "path index init:" << mytime()-t0 << "\n";

//     for ( pathtoalnmap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
//         std::string query = (*it).second->getconsensusstring();
//         kmerarray kmers = biostr::getkmers( query, filter_kmer );
//         pathid pid = it->first;
//         path_index.add(kmers, pid);
//     }
//     if ( param.verbose ) std::cout << "path index construction:" << mytime()-t0 << "\n";
    
//     /* keep track of merged paths & success paths */
//     pathidset merged_paths; 

//     /* auxilary variables for progress log */
//     double pratio = 1.0;
//     int    nproc  = 0;
//     double lt0    = mytime();
//     size_t npath  = path2aln_map.size();

//     /* make a path length map */
//     pathlengthmap plen_map = getpathlengths( path2aln_map );
//     for ( pathlengthmap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
//         pathid sbjct_pid = it->second;
//         //if ( it->first < minlen ) break;

//         if ( merged_paths.find( sbjct_pid ) != merged_paths.end() ) continue;

//         std::string sbjct = path2aln_map[sbjct_pid]->getconsensusstring();
//         assert(sbjct.size() > 0);
//         kmerarray sbjct_kmers = biostr::getkmers( sbjct, filter_kmer );

//         /* remove itself from index */
//         path_index.remove( sbjct_kmers, sbjct_pid );

//         /* print progress status at every percent of path proceeded */
//         nproc++;
//         if ( nproc/(double)npath*100 >= pratio ) {
//             fprintf( stdout, "%d/%d: %.2f sec\n", nproc, int(npath), mytime()-lt0);
//             pratio += 1;
//         }
        
//         if ( param.verbose ) {
//             std::cout << "\ncenter:" << sbjct_pid << "\tlength:" << it->first << "\n";
//             std::cout << sbjct << "\n";
//         }
        
//         pathalignpairlist left_pairs, right_pairs;        
//         pathidset skip_pids;
//         if ( param.verbose ) std::cout << "scanning left\n";
//         findlatchablepath( 1, sbjct_pid, kmerarray(sbjct_kmers.begin(), sbjct_kmers.begin()+sbjct_kmers.size()/2), path_index, skip_pids, path2aln_map, left_pairs, pairs, strands, nreads, left, param, filter_kmer );
//         if ( param.verbose ) std::cout << "scanning right\n";
//         findlatchablepath( 1, sbjct_pid, kmerarray(sbjct_kmers.begin()+sbjct_kmers.size()/2, sbjct_kmers.end()), path_index, skip_pids, path2aln_map, right_pairs, pairs, strands, nreads, right, param, filter_kmer );
//         if ( param.verbose ) std::cout << "scanning done\n";
//         if ( left_pairs.size() == 0 && right_pairs.size() == 0 ) {
//             if ( param.verbose ) std::cout << "zero latchable paths\n";
//             continue;
//         }
        
//         if ( param.verbose ) std::cout << "left:" << left_pairs.size() << "\tright:" << right_pairs.size() << "\n";
        
//         bool found = false;
//         if ( left_pairs.size() || right_pairs.size() ) found = true;
        
//         PathAlignPairList::iterator jt; 
//         while (1) {
//             if ( left_pairs.size() == 0 ) break;
//             //*jt = left_pairs.front();
//             //PathId query_pid = jt->first;
//             //AlignSummary summary = jt->second;
//             PathId query_pid = left_pairs.front().first;
//             AlignSummary summary = left_pairs.front().second;
//             if ( param.verbose ) std::cout << "Query:" << query_pid << "\n";
//             double tms = mytime();
//             path2aln_map[sbjct_pid]->join( path2aln_map[query_pid], summary.range.first, summary.lgap.first, summary.egap.first, summary.ilist, summary.dlist, param.kmer_size, iindex, bstrs, param );
//             merg_time += mytime()-tms;

//             std::string query = path2aln_map[query_pid]->getConsensusString();
//             KmerArray query_kmers = biostr::getKmers( query, param.filter_kmer );
//             path_index.remove( query_kmers, query_pid );

//             left_pairs.pop_front();
//             merged_paths.insert(query_pid);
//         }

//         if ( param.verbose ) std::cout << "Go to right\n";
//         while (1) {
//             if (right_pairs.size() == 0 ) break;
//             //*jt = right_pairs.front();
//             //PathId query_pid = jt->first;
//             //AlignSummary summary = jt->second;
//             PathId query_pid = right_pairs.front().first;
//             AlignSummary summary = right_pairs.front().second;

//             if ( param.verbose ) std::cout << "Query:" << query_pid << "\n";
//             double tms = mytime();
//             path2aln_map[sbjct_pid]->join( path2aln_map[query_pid], summary.range.first, summary.lgap.first, summary.egap.first, summary.ilist, summary.dlist, param.kmer_size, iindex, bstrs, param );
//             merg_time += mytime()-tms;

//             std::string query = path2aln_map[query_pid]->getConsensusString();
//             KmerArray query_kmers = biostr::getKmers( query, filter_kmer );
//             path_index.remove( query_kmers, query_pid );

//             right_pairs.pop_front();
//             merged_paths.insert(query_pid);
//         }
        
//         if ( found ) {
//             int error = 0;
//             MSA msa(*path2aln_map[sbjct_pid], sbjct_pid, used_reads, bstrs, strands, pairs, PILE|FLAT|TRIM|CODON, error, param);
//             msa.printAlignment(std::cout, *path2aln_map[sbjct_pid], 100);
            
//             success_paths.push_back(sbjct_pid);
//         }
//     }

//     if ( param.verbose ) std::cout << "Paths latched:" << mytime()-ts << " sec\n";
//     if ( param.verbose ) std::cout << "Join time:" << merg_time << "\n";
//     /* Remove merged paths from the path list */
//     dropMergedPaths(merged_paths, path2aln_map);
    
//     if ( param.verbose ) std::cout << "#Total paths:" << path2aln_map.size() << "\n";

//     return PathIdSet( success_paths.begin(), success_paths.end() );
// }

PathIdSet assem::mergePaths(PathIdSet &path_ids,
                            PathToAlnMap &path2aln_map, 
                            InvertedIndex &iindex,
                            PathId *used_reads,
                            BitString *bstrs,
                            char *strands,
                            ReadId *pairs,
                            int nreads,
                            Param &param,
                            int join_type,
                            bool fixed,
                            bool topdown)
{
    double ts = mytime();

    double time_merge = 0;

    PathIdList success_paths;

    size_t npath = path_ids.size();

    setbuf(stdout, NULL); // no buffering
    if ( param.verbose ) std::cout << "#Paths to merge/extend:" << npath << "\n";
    if ( param.verbose ) std::cout << "#Total paths:" << path2aln_map.size() << "\n";

    /* Initialize Index */
    /* Avoid large hash miss */
    //size_t count = 1;
    KmerPathIndex path_index;
    double t0 = mytime();
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
        size_t nk = (*it).second->getKmerCount();
        assert(nk>0);
        //std::cout << "#kmers:" << nk << "\t";
        KmerId* ks = (*it).second->getKmers();
        std::string query = biostr::getSequenceString(ks, nk, param.kmer_size);
        assert(query.size()>0);
//         if ( param.verbose ) {
//             std::cout << count++ << "\t" << it->first;
//             std::cout << str << "\n";
//         }
        
//         std::string query = (*it).second->getConsensusString();
//         assert( query.size() > 0 );
        KmerArray kmers = biostr::getKmers( query, param.filter_kmer );
        path_index.init(kmers);
    }
    if ( param.verbose ) std::cout << "Path index init:" << mytime()-t0 << "\n";

    /* Use this centroid only */
    if ( fixed ) {
        t0 = mytime();
        for ( PathIdSet::iterator it = path_ids.begin(); it != path_ids.end(); ++it ) {
            size_t nk = path2aln_map[*it]->getKmerCount();
            assert(nk>0);
            //std::cout << "#kmers:" << nk << "\t";
            KmerId* ks = path2aln_map[*it]->getKmers();
            std::string query = biostr::getSequenceString(ks, nk, param.kmer_size);
            assert(query.size()>0);
            //std::string query = path2aln_map[*it]->getConsensusString();
            KmerArray kmers = biostr::getKmers( query, param.filter_kmer );
            PathId pid = *it;
            path_index.add(kmers, pid);
        }
        if ( param.verbose ) std::cout << "Path index construction:" << mytime()-t0 << "\n";
    }

    /* auxilary variables for progress log */
    double pratio = 1.0;
    int    nproc  = 0;
    double lt0    = mytime();

    /* Keep track of merged paths */
    PathIdSet merged_paths; 

    /* Bins */
    PathBins path_bins;

    size_t min_length = 60;

    size_t success = 0;
    /* Make a path length map */
    PathLengthMap plen_map = getPathLengths( path2aln_map );
    for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        PathId query_pid = it->second;

        if (topdown && merged_paths.find(query_pid) != merged_paths.end() ) continue;

        std::string query = path2aln_map[query_pid]->getConsensusString();
        assert(query.size() > 0);
        KmerArray query_kmers = biostr::getKmers( query, param.filter_kmer );

        if ( topdown && query.size() < min_length ) break;

        /* Print progress status at every percent of path proceeded */
        nproc++;
        if ( nproc/(double)npath*100 >= pratio ) {
            fprintf( stdout, "%d/%d: %.2f sec\n", nproc, int(npath), mytime()-lt0);
            pratio += 1;
        }
        
        //std::string query = biostr::getSequenceString(path2aln_map[query_pid]->getKmers(), path2aln_map[query_pid]->getKmerCount(), param.kmer_size);                
        if ( param.verbose ) {
            std::cout << "\nQuery:" << query_pid << "\tLength:" << it->first << "\n";
            std::cout << query << "\n";
        }
        
        size_t nkmers = path2aln_map[query_pid]->getKmerCount();
        assert(nkmers>0);
        int min_filter;
        if ( topdown ) min_filter = 1;
        else {
            if ( join_type == PATHMERGE ) {
                min_filter = filter::minSameKmerCount( nkmers+param.kmer_size-1, param.filter_kmer, 1-param.filter_score );
                if ( min_filter < 1 ) min_filter = 1;
            }
            else min_filter = 1;
        }

        if (param.verbose) std::cout << "Min k-mer:" << min_filter << "\n";
        
        PathFreqMap path_freq;
//         for ( PathBins::iterator it = path_bins.begin(); it != path_bins.end(); ++it )
//             path_freq[it->first] = 0;
        
        double t0 = mytime();
        path_index.findPaths( path_freq, query_kmers );
        std::multimap<size_t, PathId> freq_map = util::sortByValue<PathId, size_t>(path_freq);
        if ( param.verbose ) std::cout << "Scan paths time:" << mytime() -t0 << "\n";
                
        bool found = false;
        //PathId target = query_pid;

        t0 = mytime();
        size_t npaths = 0;
        for ( std::multimap<size_t, PathId>::reverse_iterator it = freq_map.rbegin(); it != freq_map.rend(); ++it ) {
            if ( (int)it->first < min_filter ) break;
            if ( it->second == query_pid ) continue;
            npaths++;
        }


        if ( param.verbose ) std::cout << "# good paths:" << npaths << "\n";
        //if ( npaths ) found = true;


        size_t niter = 0;
        for ( std::multimap<size_t, PathId>::reverse_iterator it = freq_map.rbegin(); it != freq_map.rend(); ++it ) {
            if ( (int)it->first < min_filter ) break;
            if ( it->second == query_pid ) continue;

            PathId sbjct_pid = it->second;
            std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();
            assert(sbjct.size()>0);

            if ( topdown ) {
                int mink = filter::minSameKmerCount( sbjct.size(), param.filter_kmer, 1-param.filter_score );
                if ( mink < 1 ) mink = 1;
                if ( (int)it->first < mink ) continue;
            }

            KmerArray sbjct_kmers = biostr::getKmers( sbjct, param.filter_kmer );

            //double tm0 = mytime();
            AlignSummary summary;
            if ( topdown ) {
                SemiGlobalAlign aln(query, sbjct, ANCHOR_CENTER, param.gap_ext, param.gap_open);
                summary = aln.getSummary();
            } else {
                SemiGlobalAlign aln(sbjct, query, ANCHOR_CENTER, param.gap_ext, param.gap_open);
                summary = aln.getSummary();
            }
            if ( param.verbose ) std::cout << ++niter << "\tpid1:" << query_pid << "\tpid2:" << sbjct_pid << "\tlen1:" << query.size() << "\tlen2:" << sbjct.size() << "\t#kmer:" << it->first << "\tscore:" << summary.posrate << "\n";
            if ( summary.posrate >= param.merge_score ) {
                found = true; 
                success++;
                if ( topdown ) {
                    merged_paths.insert( sbjct_pid );
                    path_index.remove( sbjct_kmers, sbjct_pid );
                } else {
                    merged_paths.insert( query_pid );
                    break;
                }
            }
            //std::pair<bool, int> result = __tryMergePathPair( query_pid, sbjct_pid, path2aln_map, merged_paths, pairs, strands, nreads, bstrs, used_reads, iindex, param, join_type );
        }

        /*
        for ( std::multimap<size_t, PathId>::reverse_iterator it = freq_map.rbegin(); it != freq_map.rend(); ++it ) {
            if ( it->first < min_filter ) break;
            PathId sbjct_pid = it->second;
            std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();
            assert(sbjct.size()>0);
            KmerArray sbjct_kmers = biostr::getKmers( sbjct, param.filter_kmer );

            double tm0 = mytime();
            std::pair<bool, int> result = __tryMergePathPair( query_pid, sbjct_pid, path2aln_map, merged_paths, pairs, strands, nreads, bstrs, used_reads, iindex, param, join_type );
            double tm1 = mytime();
            time_merge += (tm1-tm0);
            
            if ( result.first ) {
                success_paths.push_back(sbjct_pid);
                path_bins[sbjct_pid].push_back( query_pid );
                target = sbjct_pid;
                found = true;
                success++;

                std::string nstr = path2aln_map[sbjct_pid]->getConsensusString();
                assert(nstr.size()>0);
                if ( sbjct != nstr ) {
                    double tu = mytime();
                    path_index.remove( sbjct_kmers, sbjct_pid );

                    sbjct_kmers = biostr::getKmers( nstr, param.filter_kmer );
                    path_index.add( sbjct_kmers, sbjct_pid );
                    if ( param.verbose ) std::cout << "udpate centroid:" << mytime() - tu << "\n";
                }
                break;
            }
        }


        if ( found ) {
            if ( param.verbose ) std::cout << "success:" << target << "\tbin size:" << path_bins[target].size() << "\t# clusters:" << path_bins.size() << "\t# kmers:" << path_index.getSize() << "\ttime:" << mytime()-t0 << " sec\n";
            //             //continue;
        } 
        */
        
        if ( found ) {
            if ( param.verbose ) std::cout << "success:"  << mytime()-t0 << " sec\n";
        }

        if ( !found ) {
            if ( !fixed ) {
                path_bins.insert( std::pair<PathId, PathAlignPairList>( query_pid, PathAlignPairList() ) );
                path_index.add( query_kmers, query_pid );
            }
            if ( param.verbose ) std::cout << "# update time:" << mytime()-t0 << "\n";
            if ( param.verbose ) std::cout << "fail:" <<  "\t# clusters:" << path_bins.size() << "\t# kmers:" << path_index.getSize()  << "\t time:" << mytime()-t0 << " sec\n";
        }


        if ( topdown ) 
            path_index.remove( query_kmers, query_pid );

    }

    if ( param.verbose ) std::cout << "\n# Success:" << success << "\n";

    /* Remove merged paths from the path list */
    if ( param.verbose ) std::cout << "# Merged:" << merged_paths.size() << "\n";
    dropMergedPaths(merged_paths, path2aln_map);
    
    if ( param.verbose ) std::cout << "#Total paths:" << path2aln_map.size() << "\n";

    double time_total = mytime()-ts;
   
    if ( param.verbose ) std::cout << "Elapsed:" << mytime()-ts << " sec\tMerge time:" << time_total - time_merge << "\n";

    //if ( param.verbose ) printPaths(path2aln_map, param);

    return PathIdSet( success_paths.begin(), success_paths.end() );
}

PathIdSet assem::latchPaths(PathIdSet &path_ids,
                            PathToAlnMap &path2aln_map, 
                            InvertedIndex &iindex,
                            PathId *used_reads,
                            BitString *bstrs,
                            char *strands,
                            ReadId *pairs,
                            int nreads,
                            Param &param,
                            int join_type )
{

    PathIdList success_paths;

    size_t npath = path_ids.size();
    if ( param.verbose ) std::cout << "#Paths to merge/extend:" << npath << "\n";
    
    PathLengthMap plen_map = getPathLengths( path2aln_map );
    
    /* auxilary variables for progress log */
    //double pratio = 1.0;
    //int    nproc  = 0;
    //double lt0    = mytime();

    /* Keep track of merged paths */
    PathIdSet merged_paths; 


//     /* Keep track of compared path pairs */
    std::tr1::unordered_map<PathId, PathIdList> compared_paths;
    std::tr1::unordered_map<PathId, bool> updated_paths;
    std::tr1::unordered_map<PathId, PathIdList> path_bins;

    //std::list<PathSize> path_sizes;
    std::multimap<size_t, PathId> path_sizes;
    KmerPathIndex path_index;
    size_t count = 0;
    double t0 = mytime();
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
        PathId pid = it->first;
        KmerArray kmers = KmerArray( (*it).second->getKmers(), (*it).second->getKmers() + (*it).second->getKmerCount() );
        path_index.add( kmers, pid );
        path_sizes.insert( std::pair<size_t, PathId> ( kmers.size()+param.kmer_size-1, pid ) );
        if ( ++count%10000 == 0 && param.verbose ) std::cout << ".";
    }
    if ( param.verbose ) std::cout << "\n";
    if ( param.verbose ) std::cout << "Built index:" << mytime() - t0 << "\n";

    count = 0;
    std::multimap<size_t, PathId>::reverse_iterator pt;
    while ( path_sizes.size() ) {
        pt = path_sizes.rbegin();
        PathId query_pid = pt->second;
        //size_t size = pt->first;
        //if ( param.verbose ) std::cout << "Curr:" << query_pid << "\tsize:" << size << "\n";

        KmerId *kmers = path2aln_map[query_pid]->getKmers();
        size_t nkmers = path2aln_map[query_pid]->getKmerCount();

        //double tt0 = mytime();
        KmerArray query_kmers = KmerArray( kmers, kmers+nkmers );
        PathFreqMap path_freq;       
        path_index.findPaths( path_freq, query_kmers );
//         for ( PathFreqMap::iterator pt = path_freq.begin(); pt != path_freq.end(); )
//             if ( pt->second < min_filter ) path_freq.erase(pt++);
//             else ++pt;

       std::multimap<size_t, PathId> freq_map = util::sortByValue<PathId, size_t>(path_freq);
        //if ( param.verbose ) std::cout << "similar paths:" << freq_map.size() << "\tsearch time:" << mytime() -tt0 << "\n";

        if ( ++count%10000 == 0 && param.verbose ) std::cout << ".";

        for ( size_t i = 0; i < nkmers; i++ ) {
            path_index.remove( kmers[i], query_pid );

        }
        path_sizes.erase((++pt).base());
        //path_sizes.erase(pt++);
    }

    if ( param.verbose ) std::cout << "\n";

    return PathIdSet( success_paths.begin(), success_paths.end() );
}

int assem::getUnusedReadCount( PathId *used_reads,
                        int nreads )
{
    int count = 0;
    for ( int i = 0; i < nreads; i++ )
        if ( used_reads[i] == NOT_PATH ) 
            count++;
    return count;
}

bool assem::latchableType(int type)
{
    switch(type) {
    case LATCH: 
        return true;
    case WEAK_LATCH:
        return true;
    case LATCH_LEFT:
        return true;
    case LATCH_LEFT_EASY:
        return true;
    case LATCH_LEFT_DIFF:
        return true;
    case LATCH_RIGHT_EASY:
        return true;
    case LATCH_RIGHT:
        return true;
    case LATCH_RIGHT_DIFF:
        return true;
    default:
        return false;
    }
}

bool assem::mergibleType(int type)
{
    switch(type) {
    case SUBSTRING:
        return true;
    case BUBBLE: 
        return true;
    case SPUR_RIGHT:
        return true;
    case SPUR_LEFT:
        return true;
    case ROPE:
        return true;
    default:
        return false;
    }
}

IntPair assem::__extractSbjectRange(std::vector<int> &pos_vec,
                                   Param &param,
                                   int anchor )
{
    IntPair rrange = IntPair(-1, -1);
    switch (anchor) {
    case ANCHOR_CENTER:
        rrange = getMatchRangeMax(pos_vec, param); 
        break;
    case ANCHOR_LEFT:
        rrange = getMatchRangeLeft(pos_vec, param);
        break;
    default:
        rrange = getMatchRangeRight(pos_vec, param);
    }
    return rrange;
}

bool assem::__validRange(IntPair &rrange,
                         std::vector<int> &pos_vec,
                         Param &param)
{
    if ( rrange.first == -1 ||  rrange.second == -1 ) return false;
    
    if ( !ordered(pos_vec, rrange) ) { 
        if ( param.verbose ) std::cout << "Not in order\n"; 
        return false; 
    }
    if ( indeled(pos_vec, rrange) ) { 
        if ( param.verbose ) std::cout << "Indeled\n"; 
        return false; 
    }
    return true;
}

bool assem::__alignReadToPath(int &type,
                              AlignSummary &summary,
                              IntPair &qrange,
                              IntPair &rrange,
                              KmerId *qkmers,
                              size_t qnkmers,
                              KmerId *skmers,
                              size_t snkmers,
                              Param &param )

{    
    adjustRanges(qrange, rrange, qnkmers, snkmers, param);
    
    std::string query = biostr::getSequenceString(&qkmers[qrange.first], qrange.second-qrange.first+1, param.kmer_size);
    std::string sbjct = biostr::getSequenceString(&skmers[rrange.first], rrange.second-rrange.first+1, param.kmer_size);
    
    if ( param.verbose ) {
        std::cout << "Query:" << query << "\n";
        std::cout << "Sbjct:" << sbjct << "\n";
    }
    
    //AlignSummary 
    summary = compareBySubstitution(query, sbjct, qnkmers, snkmers, qrange.first, qrange.second, rrange.first, rrange.second, param);
    if ( latchableType(type) ) {
        if ( summary.posrate < param.latch_score ) {
            if ( param.verbose ) std::cout << "Weak score\n";
            return false;
        }
    } else {
        if ( summary.posrate < param.merge_score ) {
            if ( param.verbose ) std::cout << "Weak score\n";
            return false;
        }
    }
    return true;
}

bool assem::__goodReadToJoin( std::vector<int> &pos_vec,
                              IntPair &qrange,
                              IntPair &rrange,
                              int &type,
                              AlignSummary &summary,
                              KmerId *qkmers,
                              size_t qnkmers,
                              KmerId *skmers,
                              size_t snkmers,
                              Param &param,
                              int anchor )
{
    rrange = __extractSbjectRange(pos_vec, param, anchor);
    if ( !__validRange(rrange, pos_vec, param) ) return false;

    qrange = IntPair( pos_vec[rrange.first], pos_vec[rrange.second] );
    type = getMatchType( qrange, rrange, qkmers, skmers, qnkmers, snkmers, param );
    
    if ( param.verbose ) 
        std::cout << "rrange:" << rrange.first << ":" << rrange.second 
                  << "\tqrange:" << qrange.first << ":" << qrange.second << "\n"
                  << "Type:" << type << "\n";
    
    if ( type == NOTYPE ) return false;

    return __alignReadToPath(type, summary, qrange, rrange, qkmers, qnkmers, skmers, snkmers, param);
}

void assem::__mergeToPath( ReadId curr_rid,
                           PathId *used_reads,
                           PathToAlnMap &path2aln_map,
                           PathId sbjct_pid,
                           IntPair &qrange,
                           IntPair &rrange,
                           AlignSummary &summary,
                           int &merged,
                           //PathIdSet &dirty_paths,
                           ReadPathPosMap &added_reads,
                           Param &param)
{
    ++merged;
    used_reads[curr_rid] = sbjct_pid;
    
//     int qs = qrange.first;
//     int rs = rrange.first;
    
//     int spos = rs;
//     if ( qs > 0 ) spos = (rs-qs);
    
//     if ( param.verbose ) std::cout << "Merged to :" << sbjct_pid << "\tSpos:" << spos << "\n";
    
//     SpaPath *spath = path2aln_map[sbjct_pid];
//     spath->addRead( curr_rid, spos );
    
//     //dirty_paths.insert(sbjct_pid);

    int spos = 0;
    if ( summary.lgap.first > 0 ) spos = -1 * summary.lgap.first;
    else spos = summary.range.first;
    
    if ( param.verbose ) std::cout << "Merged to :" << sbjct_pid << "\tSpos:" << spos << "\n";

    added_reads[sbjct_pid].push_back( std::pair<ReadId, int>( curr_rid, spos ) );
}

void assem::recruit( BitString *bstrs,
                     PathId *used_reads,
                     ReadId urid,
                     PathToAlnMap &path2aln_map,
                     KmerToPathMap &pathid_map,
                     std::list<LatchRead> &latch_reads,
                     int &merged,
                     ReadPathPosMap &added_reads,
                     //PathIdSet &dirty_paths,
                     Param &param
                     )
{
    double tic = mytime();

    PosPathPairList pp_list;
    std::list<IntPair> qranges,  rranges;
        
    std::string str = bstrs[urid].toString();
    if ( (int)str.size() < param.kmer_size || !seq::validSequence(str) ) {
        if (param.verbose) std::cout << "Invalid Read Sequence\n"; 
        return;
    }

    if ( param.verbose ) std::cout << "\nReadId:" << urid << "\tSequence:" << str << "\n";
    KmerArray kmers = biostr::getKmers( str, param.kmer_size );

    int mink = 1;
    if ( param.pair_flag ) {
        mink = filter::minSameKmerCount( str.length(), param.kmer_size, 1-READ_FILTER );
        if ( mink < 1 ) mink = 1;
    }
    PathIdSet merged_paths;
    PosPathPairList pos_paths = searchSimilarPaths(PathIdSet(), &kmers[0], kmers.size(), path2aln_map, pathid_map, merged_paths, mink, true, param);
    if ( param.verbose ) std::cout << "Min kmer:" << mink << "\t# Similar paths:" << pos_paths.size() << "\n";

    bool merge = false;
    PosPathPairList::iterator jt;
    int i = 0;
    for ( jt = pos_paths.begin(); jt != pos_paths.end(); ++jt ) {
        if (param.verbose) {
            std::cout << "\n" << ++i << "\tPathId:" << jt->second << "\n";
            std::cout << "Sequence:" << path2aln_map[jt->second]->getConsensus() << "\n";
        }

        KmerId *rkmers = path2aln_map[jt->second]->getKmers();
        size_t rnkmers = path2aln_map[jt->second]->getKmerCount();

        if ( rnkmers < kmers.size() ) {
            if (param.verbose) std::cout << "Too short reference\n";
            continue;
        }

        IntPair qrange, rrange;
        int type;
        AlignSummary summary;

        /* check max block */
        if ( __goodReadToJoin( jt->first, qrange, rrange, type, summary, &kmers[0], kmers.size(), rkmers, rnkmers, param, ANCHOR_CENTER ) ) {
            if ( mergibleType(type) ) {
                __mergeToPath(urid, used_reads, path2aln_map, jt->second, qrange, rrange, summary, merged, added_reads, param);
                merge = true;
                if ( param.verbose ) std::cout << "time:" << mytime() - tic << " sec\n";//tpaths:" << pos_paths.size() << "\n";
                break;
            }
            
            else if ( type == LATCH_LEFT_EASY || type == LATCH_RIGHT_EASY ) {
                pp_list.push_back(*jt); qranges.push_back(qrange); rranges.push_back(rrange);
                continue;
            }
        }
        
//         /* check left block */
//         if ( __goodReadToJoin( jt->first, qrange, rrange, type, &kmers[0], kmers.size(), rkmers, rnkmers, param, ANCHOR_LEFT ) ) {
//             if ( type == LATCH_LEFT_EASY ) {
//                 pp_list.push_back(*jt); qranges.push_back(qrange); rranges.push_back(rrange);
//                 continue;
//             }
//         }
        
//         /* check right block */
//         if ( __goodReadToJoin( jt->first, qrange, rrange, type, &kmers[0], kmers.size(), rkmers, rnkmers, param, ANCHOR_RIGHT ) ) {
//             if ( type == LATCH_RIGHT_EASY ) {
//                 pp_list.push_back(*jt); qranges.push_back(qrange); rranges.push_back(rrange);
//                 continue;
//             }
//         }
    }

    if ( !merge ) {
        PosPathPairList::iterator pt = pp_list.begin();
        std::list<IntPair>::iterator qt = qranges.begin();
        std::list<IntPair>::iterator rt = rranges.begin();
        for ( ; ; ) {
            if ( pt == pp_list.end() ) break;
            latch_reads.push_back( LatchRead(urid, pt->second, qt->first, rt->first ) );
            ++pt; ++qt; ++rt;
        }
    }
    if ( param.verbose ) {
        if ( !merge ) {
            std::cout << "# Latchable paths:" << pp_list.size() << "\n";
            std::cout << "time:" << mytime() - tic << " sec\n";//tpaths:" << pos_paths.size() << "\n";
            //std::cout << "\n";
        }
    }
}

//         //IntPair rrange = getMatchRange(jt->first, param);
//         /* Inspect only maximum block */
//         IntPair rrange = getMatchRangeMax(jt->first, param);
//         KmerId *rkmers = path2aln_map[jt->second]->getKmers();
//         size_t rnkmers = path2aln_map[jt->second]->getKmerCount();

//         if ( rnkmers < kmers.size() ) {
//             if (param.verbose) std::cout << "Too short reference\n";
//             continue;
//         }
//         std::string region = biostr::getSequenceString(&rkmers[rrange.first], rrange.second-rrange.first+1, param.kmer_size);
//         if ( param.verbose ) {
//             std::cout << "Share:" << region << "\n";
//         }

//         IntPair qrange;
//         qrange.first  = (*jt).first[rrange.first];
//         qrange.second = (*jt).first[rrange.second];

//         //int type = getLinkType( qrange, rrange, &kmers[0], rkmers, kmers.size(), rnkmers, param.latch_length, param.read_spur, param );
//         //int type = getPairType( qrange, rrange, &kmers[0], rkmers, kmers.size(), rnkmers, param.latch_length, param.read_spur, param );
//         //int type = getPathPairType( jt->first, qrange, rrange, &kmers[0], rkmers, kmers.size(), rnkmers, param.latch_length, param.read_spur, param );
//         int type = getMatchType( qrange, rrange, &kmers[0], rkmers, kmers.size(), rnkmers, param );
//         if ( param.verbose ) {
//             std::cout << "rrange:" << rrange.first << ":" << rrange.second 
//                       << "\tqrange:" << qrange.first << ":" << qrange.second << "\n";
//            std::cout << "Type:" << type << "\n";
//         }
        
//         if ( type == NOTYPE ) continue;

//         if ( !ordered(jt->first, rrange) ) { if ( param.verbose ) std::cout << "Not in order\n"; continue; }
//         if ( indeled(jt->first, rrange) ) { if ( param.verbose ) std::cout << "Indeled\n"; continue; }
            

//         int qnkmers = kmers.size();
// //         if ( type == SPUR || type == WEAK_SPUR || type == ROPE || type == WEAK_ROPE ) 
// //             adjustRanges(qrange, rrange, qnkmers, rnkmers, param);
// //         else if ( type == LATCH ) adjustLatchRange(qrange, rrange, qnkmers, rnkmers, param );
// //        adjustRanges(qrange, rrange, qnkmers, rnkmers, param);

//         std::string query = biostr::getSequenceString( &kmers[qrange.first], qrange.second-qrange.first+1, param.kmer_size);
//         std::string sbjct = biostr::getSequenceString(&rkmers[rrange.first], rrange.second-rrange.first+1, param.kmer_size);

//         if ( param.verbose ) {
//             std::cout << "Query:" << query << "\n";
//             std::cout << "Sbjct:" << sbjct << "\n";
//         }

//         AlignSummary summary = compareBySubstitution(query, sbjct, qnkmers, rnkmers, qrange.first, qrange.second, rrange.first, rrange.second, param);
//         //if ( type == LATCH ) {
//         if ( latchableType(type) ) {
//             if ( summary.posrate < param.latch_score ) { if ( param.verbose ) std::cout << "Weak score\n"; continue; }
//         } else {
//             if ( summary.posrate < param.merge_score ) { if ( param.verbose ) std::cout << "Weak score\n"; continue; }
//         }

//         success = true;
//         pp_list.push_back(*jt);
//         qranges.push_back(qrange);
//         rranges.push_back(rrange);

//         //if ( type == LATCH ) latch = true;
//         if ( latchableType(type) ) latch = true;
//         else { merge = true; mid = jt->second; break; }
//     }

//     if ( merge ) {
//         ++merged;
//         used_reads[urid] = mid;

//         PosPathPair pp = pp_list.back();
//         int qs = qranges.back().first;
//         int rs = rranges.back().first;

//         int spos = rs;
//         if ( qs > 0 ) spos = (rs-qs);
            
//         if ( param.verbose ) std::cout << "PathId:" << pp.second << "\tSpos:" << spos << "\n";
            
//         SpaPath *spath = path2aln_map[pp.second];
//         spath->addRead( (ReadId)urid, spos );


//         dirty_paths.insert(pp.second);

//         if (param.verbose) std::cout << "time:" << mytime() - tic << " sec\tpaths:" << pos_paths.size() << "\tmerged\n";

//         return;
//     } 
//     if ( latch ) {
//         PosPathPairList::iterator pt = pp_list.begin();
//         std::list<IntPair>::iterator qt = qranges.begin();
//         std::list<IntPair>::iterator rt = rranges.begin();
//         for ( ; ; ) {
//             if ( pt == pp_list.end() ) break;
//             latch_reads.push_back( LatchRead(urid, pt->second, qt->first, rt->first ) );
//             ++pt; ++qt; ++rt;
//         }
//     }
//     if ( param.verbose ) {
//         std::cout << "time:" << mytime() - tic << " sec\tpaths:" << pos_paths.size() << "\n";
//         std::cout << "\n";
//     }
// }

void assem::addReadsToPath( PathId pid,
                            ReadPathPosList rposs,
                            PathToAlnMap &path2aln_map,
                            BitString *bstrs,
                            PathId *used_reads,
                            Param &param )
{
    double t0 = mytime();

    std::vector<ReadId> nreads;
    std::vector<int> sposs;
    int nrid = 0;

    //Profile profile;
    ProfileVector *vprof = path2aln_map[pid]->getProfileVector();
    //if ( param.verbose ) std::cout << "Profile vector size:" << vprof->getSize() << "\n";

    //vprof->convert(profile);
    Profile profile = vprof->convert();

    if ( profile.ncol == 0 ) {
        std::cout << "Zero length profile\n";
        std::cout << "Path:" << pid << "\n";
        std::cout << "vector size:" << vprof->getSize() << "\n";
        std::cout << "# Reads:" << path2aln_map[pid]->getReadCount() << "\n";
        std::cout << "Consensus:" << path2aln_map[pid]->getConsensus() << "\n";
        std::cout << "New Reads:" << rposs.size() << "\n";
    }

    std::string consensus = path2aln_map[pid]->getConsensusString();
    for ( ReadPathPosList::iterator jt = rposs.begin(); jt != rposs.end(); ++jt ) {
        ReadId rid = jt->first;
        int pos = jt->second;
        std::string read = bstrs[rid].toString();
        if ( param.verbose ) std::cout << "Adding new read:" << rid << "\tat:" << pos << "\t" << read << "\n";

        if ( consensus[consensus.size()-1] == '*' && read[read.size()-1] == '*' ) {
            int qend = pos + read.size()-1;
            int send = consensus.size()-1;
            if ( qend < send ) {
                if ( param.verbose ) std::cout << "Another stop codon in another pos\n";
                used_reads[rid] = NOT_PATH;
                continue;
            }
        }
        
        for ( size_t i = 0; i < read.size(); i++ ) {
            char aa = read[i];
            int na = alpha::AsciiToAA[(int)aa];
            if ( na < 1 || na > 27 ) {
                std::cerr << "Invalid aa:" << aa << "\tnum:" << na << "\tcol:" << i << "\t" << read << "\n";
                exit(1);
            }
            int mpos = pos + i;
            if ( mpos < 0 ) continue;
            if ( mpos > (int)profile.ncol ) {
                break; // temporary solution
            }
            profile.matrix[mpos][na-1]++;
        }

        used_reads[rid] = pid;
        nreads.push_back(rid);
        sposs.push_back(pos);
        nrid++;
    }
    if ( nrid > 0 ) {
        path2aln_map[pid]->addReads( &nreads[0], &sposs[0], nrid );
        path2aln_map[pid]->resetProfile();
        path2aln_map[pid]->setProfile(profile);
        std::string ncon = profile.makeConsensus();
        if ( consensus != ncon ) {
            if ( param.verbose ) {
                std::cout << "Old consensus:" << consensus << "\n";
                std::cout << "New consensus:" << ncon << "\n";
            }
            KmerArray kmers = biostr::getKmers( ncon, param.kmer_size );
            path2aln_map[pid]->updateConsensusSequence( ncon.c_str() );
            path2aln_map[pid]->updateKmers( &kmers[0], kmers.size() );
        }
        else 
            if ( param.verbose ) std::cout << "Same consensus\n";

        if ( param.verbose ) {
            std::cout << "\nAlignment\n";
            path2aln_map[pid]->printAlignment(std::cout, 100, bstrs);
        }
    }
    if ( param.verbose ) std::cout << "\nRead joined:" << mytime()-t0 << "\n";
}

void assem::simplePathClean( ReadPathPosMap &added_reads,                            
                             PathIdSet &merged_paths, 
                             PathToAlnMap &path2aln_map,
                             PathId *used_reads,
                             BitString *bstrs,
                             char *strands, 
                             ReadId *pairs, 
                             Param &param)
{
    //double t0 = mytime();
    if ( param.verbose ) std::cout << "\n# Dirty paths:" << added_reads.size() << "\n";
    
    for ( ReadPathPosMap::iterator it = added_reads.begin(); it != added_reads.end(); ++it ) {
        if ( param.verbose ) std::cout << "\nDirty path:" << it->first << "\t#new reads:" << it->second.size() << "\n";

        if ( merged_paths.find(it->first) != merged_paths.end() ) {
            if ( param.verbose ) std::cout << "Merged -> skip\n"; 
            continue;
        }
        if ( path2aln_map.find(it->first) == path2aln_map.end() ) {
            if ( param.verbose) std::cout << "Path does not exist\n";
            //exit(1);
            continue;
        }

        addReadsToPath( it->first, it->second, path2aln_map, bstrs, used_reads, param );
    }
    //if ( param.verbose ) 
    // std::cout << "Path cleaned:" << mytime()-t0 << "\n";
}

void assem::cleanDirtyPaths( PathIdSet &dirty_paths,
                             PathIdSet &merged_paths, 
                             PathToAlnMap &path2aln_map,
                             PathId *used_reads,
                             BitString *bstrs,
                             char *strands, 
                             ReadId *pairs, 
                             Param &param)
{
    double t0 = mytime();
    if ( param.verbose ) std::cout << "\n# Dirty paths:" << dirty_paths.size() << "\n";
    
    for ( PathIdSet::iterator it = dirty_paths.begin(); it != dirty_paths.end(); ++it ) {
        if ( param.verbose ) std::cout << "Dirty path:" << *it << "\n";
        if ( merged_paths.find(*it) != merged_paths.end() ) {
            std::cout << "Merged -> skip\n"; 
            continue;
        }
        if ( path2aln_map.find(*it) == path2aln_map.end() ) {
            if ( param.verbose) std::cout << "Path does not exist\n";
            //exit(1);
            continue;
        }
        SpaPath *spath = path2aln_map[*it];
        int error = 0;
		MSA msa(*spath, *it, used_reads, bstrs, strands, pairs, PILE|FLAT|TRIM|CODON, error, param);
        
        if ( error ) {
            if ( param.verbose ) std::cout << "MSA fail:\tcode:" << error << "\n";
            
            path2aln_map[*it]->resetReads(used_reads, NOT_PATH);
            path2aln_map.erase(*it);
            delete spath;
            continue;
        }
    }
    dirty_paths.clear();
    if (param.verbose) std::cout << "Path cleaned:" << mytime()-t0 << " sec\n";
}

void assem::recruitReads( PathToAlnMap &path2aln_map,
                   InvertedIndex &iindex,
                   BitString *bstrs,
                   char *strands,
                   ReadId *pairs,
                   PathId *used_reads,
                   size_t nreads,
				   Param &param )
{
    if ( !param.recruit_flag ) return;

    double t0 = mytime();
    
    setbuf(stdout, NULL); // no buffering
    //std::cout << "\nSTAGE 3:\nRecruiting reads ...\n";

    int unused = getUnusedReadCount( used_reads, nreads );
    std::cout << "Unused read:" << unused << "\n";
    if ( unused == 0 ) return ;
        

    //-------------------------
    // k-mer to path id mapping
    //-------------------------
    KmerToPathMap pathid_map;
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) 
        addPathIds( it->first, path2aln_map[it->first]->getKmers(), path2aln_map[it->first]->getKmerCount(), pathid_map );

    
    /*=============================================================================
     * If a read can be merged to an existing path, join it to the path.
     * If a read can be latched into existing paths, record the read and the paths.
     *===========================================================================*/
    double pratio = 1.0;
    int nproc = 0;

    ReadPathPosMap added_reads;
    //PathIdSet dirty_paths;
    std::list<LatchRead> latch_reads; // record latchable reads to be proceed later
    int merged = 0; // keep track of no. of reads merged
    double lt0 = mytime();
    for ( size_t urid = 0; urid < nreads; urid++ ) {
        if ( used_reads[urid] != NOT_PATH ) continue;

        nproc++;
        if ( nproc/(double)unused*100 >= pratio ) {
            fprintf( stdout, "%d/%d: %.2f sec\n", nproc, unused, mytime()-lt0);
            pratio += 1;
        }
        recruit( bstrs, used_reads, urid, path2aln_map, pathid_map, latch_reads, merged, added_reads, param);
    }

    if ( param.verbose ) std::cout << "\nScanned unused reads:" << mytime()-t0 << " sec\n";

    /* Merge reads on paths */
    PathIdSet merged_paths;
    //cleanDirtyPaths(dirty_paths, path2aln_map, used_reads, bstrs, strands, pairs, param);
    simplePathClean(added_reads, merged_paths, path2aln_map, used_reads, bstrs, strands, pairs, param);

    //-------------------------------
    // Latch paths with read supports
    //------------------------------- 
    PathIdSet successPaths;
    PathIdSet dirty_paths;
    if ( !param.pair_flag ) {
        extendPathsByBridgingReads( latch_reads, used_reads, path2aln_map, bstrs, strands, pairs, pathid_map, merged_paths, successPaths, dirty_paths, iindex, param );
//     else
//         extendPairedPathsByBridgingReads( used_reads, path2aln_map, bstrs, strands, pairs, pathid_map, merged_paths, successPaths, dirty_paths, iindex, param );

        /* Merge reads remained after path joining */
        cleanDirtyPaths(dirty_paths, merged_paths, path2aln_map, used_reads, bstrs, strands, pairs, param);
        
        dropMergedPaths(merged_paths, path2aln_map);
        pathid_map.clear();
    }


//     //--------------------------------------------------------------------
//     // Try to merge or latch similar pairs from read supported latch paths
//     //--------------------------------------------------------------------
//     double lt = mytime();
//     linkPaths( successPaths, path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param, MERGELATCH );
//     std::cout << "Paths merged/extended with newly extended paths:" << mytime()-lt << " sec\n";

	/* Recuit summary */
	int unused2 = getUnusedReadCount( used_reads, nreads );
	std::cout << "Path count:" << path2aln_map.size() << "\n";
	std::cout << "Recruited:" << unused-unused2 << "\n";
	std::cout << "Unrecruited:" << unused2 << "\n";
	std::cout << "Read recruiment:" << mytime()-t0 << " sec\n";        
}



ReadIdList assem::getUnusedReads( PathId *used_reads,
                           int nreads )
{
    ReadIdList unused_reads; 
    
    for ( int i = 0; i < nreads; i++ )
        if ( used_reads[i] == NOT_PATH ) 
            unused_reads.push_back(i);
    return unused_reads;
}

ReadIdArray assem::getReadIds( std::list<ReadPos> &rpos_list )
{
    ReadIdList rids;
    std::list<ReadPos>::iterator it = rpos_list.begin();
    for ( ; it != rpos_list.end(); ++it )
        rids.push_back(it->rid);
    return ReadIdArray( rids.begin(), rids.end() );
}

std::tr1::unordered_map<ReadId, int> assem::getReadStartPos( ReadPosList rp, int which )
{
    std::tr1::unordered_map<ReadId, int> pmap;
    for ( ReadPosList::iterator it = rp.begin(); it != rp.end(); ++it )
        // sbjct: xxxxxxxxxxooooo
        //  read:           oooooxxxxx
        if ( which == LEFT ) 
            pmap.insert( std::pair<ReadId, int>( it->rid, it->rpos ) );
        else 
        // sbjct:      oooooxxxxxxxxxx
        //  read: xxxxxooooo
            pmap.insert( std::pair<ReadId, int>( it->rid, it->spos ) );
    return pmap;
}

void assem::makePathPairs( PathPairMap &pair_map,
                           PathIdList &left_paths, 
                           PathIdList &right_paths )
{
    for ( PathIdList::iterator lt = left_paths.begin(); lt != left_paths.end(); ++lt ) 
        for ( PathIdList::iterator rt = right_paths.begin(); rt != right_paths.end(); ++rt ) 
            pair_map[*lt].insert(*rt);
}

void assem::pairUpPaths( PathPairMap &pair_map,
                         PathReadMap &latch_left,
                         PathReadMap &latch_right,
                         std::list<LatchRead> &latch_reads,
                         Param &param )
{
    if ( latch_reads.size() == 0 ) return;

    double t0 = mytime();
    if ( param.verbose ) std::cout << "\nParing paths...\n";
    if ( param.verbose ) std::cout << "Latching reads:" << latch_reads.size() << "\n";


    PathIdList left_paths, right_paths;
    std::list<LatchRead>::iterator it = latch_reads.begin();
    ReadId curr = it->rid;
    for ( ; it != latch_reads.end(); ++it ) {
        /* Update if working read id is not identical to the new read id */ 
        if ( it->rid != curr ) {
            if ( param.verbose ) std::cout << "Read:" << curr << "\t";
            makePathPairs(pair_map, left_paths, right_paths);
            if ( param.verbose ) std::cout << "path paired\n";
            curr = it->rid;
            /* reset */
            left_paths.clear(); right_paths.clear();
        }

        /* latch right  wrt sbjct */
        // sbjct: xxxxxxxxxxoooooooooo
        // read :           ooooooooooxxxxxxxxxx
        if ( it->rpos > 0 ) {
            left_paths.push_back(it->pid);
            latch_right[it->pid].push_back( ReadPos(it->rid, it->spos, it->rpos) );
        }
        else {
            right_paths.push_back(it->pid);
            latch_left[it->pid].push_back( ReadPos(it->rid, it->spos, it->rpos) );
        }
    }
    if ( param.verbose ) std::cout << "Read:" << curr << "\t";
    makePathPairs(pair_map, left_paths, right_paths);
    if ( param.verbose ) std::cout << "path paired\n";
                
    if ( param.verbose ) {
        std::cout << "Path Pairs:" << pair_map.size() << "\n";
        std::cout << "Latch Right:" << latch_right.size() << "\n";
        std::cout << "Latch Left:" << latch_left.size() << "\n";
//         PathPairMap::iterator pt = pair_map.begin();
//         for ( ; pt != pair_map.end(); ++pt ) {
//             std::cout << pt->first << "\t";
//             for ( PathIdSet::iterator st = (*pt).second.begin(); st != (*pt).second.end(); ++st )
//                 std::cout << *st << " ";
//             std::cout << "\n";
//         }
        std::cout << "Path paired:" << mytime()-t0 << "\n";
    }
}

void assem::commonReadSet( PathId left,
                    PathPairMap &pair_map,
                    PathReadMap &latch_left,
                    PathReadMap &latch_right,
                    std::map<int, PathId> &commpid_map, 
                    std::tr1::unordered_map<PathId, ReadIdArray> &commset_map,
					Param &param )
{
    ReadIdArray lrids = getReadIds(latch_right[left]);
    std::sort(lrids.begin(), lrids.end());
    PathIdSet pset = pair_map[left];
    for ( PathIdSet::iterator pt = pset.begin(); pt != pset.end(); ++pt ) {
        ReadIdArray rrids = getReadIds(latch_left[*pt]);
        std::sort(rrids.begin(), rrids.end());
        ReadIdArray comm = Set::Intersect(rrids, lrids, true);
        if ( param.verbose ) std::cout << "common:" << comm.size() << "\n";
        
        commpid_map.insert( std::pair<int, PathId>( comm.size(), *pt ) );
        commset_map.insert( std::pair<PathId, ReadIdArray>( *pt, comm ) );
    }
}

void assem::validLatchReads( PathId *used_reads,
                             BitString *bstrs,
                             std::string &lstr,
                             std::string &rstr,
                             PathId left,
                             PathId right,
                             PathReadMap &latch_left,
                             PathReadMap &latch_right,
                             ReadIdArray &comm_reads,
                             ReadIdArray &good_reads,
                             std::vector<int> &read_inits,
                             int &overlap,
                             Param &param
                             )
{
    //bool setflag = false;
    //std::multi_map<int, ReadId> overlap_rids;
    std::map<int, std::list<ReadId> > overlap_rids;
    std::map<int, int> overlap_nums;

    // sbjct: xxxxxxxxxxxxxxxxxxxx        
    // lrids:           ------------
    //                     ----------
    // sbjct:           xxxxxxxxxxxxxxxxxxxx        
    // rrids:  ------------
    //             ----------
    std::tr1::unordered_map<ReadId, int> left_rids  = getReadStartPos( latch_right[left], LEFT );   // start pos in left path
    std::tr1::unordered_map<ReadId, int> right_rids  = getReadStartPos( latch_left[right], RIGHT ); // start pos in read
    for ( ReadIdArray::iterator cit = comm_reads.begin(); cit != comm_reads.end(); ++cit ) {
        if ( used_reads[*cit] != NOT_PATH ) continue;

        std::string read = bstrs[*cit].toString();
        assert(read.size());

        if ( left_rids.find(*cit) == left_rids.end() || right_rids.find(*cit) == right_rids.end() ) continue;
        
        if ( param.verbose ) std::cout << "Read:" << *cit << "\t" << read << "\tleft:" << left_rids[*cit] << "\tright:" << right_rids[*cit] << "\n";


        // ----------*-----
        //           | lsub
        int len = lstr.length()-left_rids[*cit];
        if ( len < 0 ) continue;
        //if (param.verbose) std::cout << "len:" << len << "\tsize:" << lstr.size() << "\n";
        if (param.verbose) std::cout << "llen:" << len << "\trlen:" << right_rids[*cit]+1 << "\tmlen:" << (int)read.size()-len-(right_rids[*cit]+1) << "\n"; //<< "\tsize:" << lstr.size() << "\n";
        if ( left_rids[*cit] + len > (int)lstr.size() ) {
            if ( param.verbose ) std::cout << "Off read:" << read << "\n"; continue;
        }
        
        assert(left_rids[*cit]>=0);
        assert(right_rids[*cit]>=0);
        assert( (int)lstr.length() > left_rids[*cit]);
        assert( (int)read.length() > right_rids[*cit]);
        std::string lsub = lstr.substr(left_rids[*cit], lstr.length()-left_rids[*cit]);
        std::string rsub = read.substr(right_rids[*cit], read.length()-right_rids[*cit]);
        if (param.verbose) std::cout << "lsub:" << lsub << "\n";
        if (param.verbose) std::cout << "rsub:" << rsub << "\n";

        // case 1. positive overlap
        // xxxxxxxxxxxxxxxooooooooooo                       lpath
        //                |  lsub   |
        //                ooooooooooooooo                   read
        //                     |olap|
        //                     ooooooooooxxxxxxxxxxxxxxx    rpath
        // case 2. negative overlap
        // xxxxxxxxxxxxxxxoooooooooo
        //                    |lsub|
        //                    ooooooxxxxxxxoooooo
        //                          |-olap|
        //                                 ooooooxxxxxxxxxxxxx
        
//         if ( !setflag ) {
//             overlap = lsub.length() - right_rids[*cit];
//             setflag = true;
//         }
        overlap = lsub.length() - right_rids[*cit];
        overlap_rids[overlap].push_back(*cit);
        if (overlap_nums.find(overlap) == overlap_nums.end())
            overlap_nums.insert(std::pair<int, int>(overlap, 0));
        overlap_nums[overlap]++;

        if (param.verbose) std::cout << "Overlap:" << overlap << "\n";
    }

    if (overlap_nums.size() == 0 ) return;

    /* no conflict */
    if ( overlap_rids.size() == 1 ) {
        if ( param.verbose ) std::cout << "No conflicting overlap\n";

        std::map<int, std::list<ReadId> >::iterator it = overlap_rids.begin();
        for ( std::list<ReadId>::iterator lt = it->second.begin(); lt != it->second.end(); ++lt ) {
            good_reads.push_back( *lt );
            read_inits.push_back( left_rids[*lt] );
            overlap = it->first;
        }
        return;
    }

    /* positive overlap */
    std::map<int, std::list<ReadId> >::reverse_iterator it;// = overlap_rids.begin();
    for ( it = overlap_rids.rbegin(); it != overlap_rids.rend(); ++it ) {
        if ( it->first > 0 && (int)it->second.size() >= param.latch_support ) {
            if ( param.verbose ) std::cout << "Positive overlap:" << it->first << "\tsize:" << it->second.size() << "\n";
            for ( std::list<ReadId>::iterator lt = it->second.begin(); lt != it->second.end(); ++lt ) {
                good_reads.push_back( *lt );
                read_inits.push_back( left_rids[*lt] );
                overlap = it->first;
            }
            return;
        }
    }


   
    /***
        do not allow conflicint extension
    */

    if ( param.verbose ) {
        std::cout << "Conflicting overlaps\n";
        for ( it = overlap_rids.rbegin(); it != overlap_rids.rend(); ++it ) {
            std::cout << it->first << "\t" << it->second.size() << "\n";
        }
    }
    

//     /* overlap  <= 0 */
//     std::multimap<int, int> count_map = util::sortByValue<int, int>(overlap_nums);
//     std::map<int,int>::reverse_iterator it = count_map.rbegin();
//     int major_olap = it->second;
    
//     for ( ; it != count_map.rend(); ++it )
//         if ( param.verbose ) std::cout << "overlap:" << it->second << "\tcount:" << it->first << "\n";

//     if ( param.verbose ) std::cout << "lstr:" << lstr << "\n";
//     if ( param.verbose ) std::cout << "rstr:" << rstr << "\n";
    
    
//     std::list<ReadId> major_rids = overlap_rids[major_olap];
//     for ( std::list<ReadId>::iterator lt = major_rids.begin(); lt != major_rids.end(); ++lt ) {
        

//         std::string lsub = lstr.substr(left_rids[*lt], lstr.length()-left_rids[*lt]);
//         std::string read = bstrs[*lt].toString();

//         if ( param.verbose ) std::cout << "rid:" << *lt << "\tpos:" << left_rids[*lt] << "\t:read:" << read << "\tlsub:" << lsub << "\n";

//         bool good = true;
//         std::string rsub;

        
//         if ( major_olap >= 0 ) {
//             int leftover = read.length()-lsub.length();
//             if ( param.verbose ) std::cout << "leftover:" << leftover << "\n";
//             if ( leftover < 0 ) { good = false; continue; }
//             if ( major_olap + leftover <= (int)rstr.length() ) 
//                 rsub = rstr.substr(major_olap,  read.length()-lsub.length());
//             else good = false;
//         }
//         else {
//             int middle = -1*major_olap;
//            if ( param.verbose )  std::cout << "middle:" << middle << "\n";
//             if ( lsub.length()+middle > read.length() ) good = false;
//             int leftover = read.length()-lsub.length()+major_olap;
//             if ( param.verbose ) std::cout << "leftover:" << leftover << "\n";
//             if ( leftover > (int)rstr.length() ) good = false;
//             //rsub = read.substr(lsub.length(), -1*major_olap) + rstr.substr(0, read.length()-lsub.length()+major_olap);
//             rsub = read.substr(lsub.length(), -1*major_olap) + rstr.substr(0, leftover);
//         }
        
//         if ( !good ) continue;
//         if ( param.verbose ) std::cout << "lsub:" << lsub << "\trsub:" << rsub << "\n";
//         std::string nstr = lsub + rsub;
//         if ( param.verbose ) {
//             std::cout << "Read:" << read << "\n";
//             std::cout << "nStr:" << nstr << "\n";
//         }

//         int npos = scoring::countPositive(read, nstr, BLOSUM62);
//         double score = (double)npos/read.length();
//         if ( param.verbose ) std::cout << "Score:" << score << "\n";
//         if ( score < param.latch_score ) {
//             if ( param.verbose ) std::cout << "Weak score\n"; 
//             continue;
//         }
        
//         good_reads.push_back( *lt );
//         read_inits.push_back( left_rids[*lt] );
//     }

// //         bool good = true;
// //         std::string rsub;
// //         if ( overlap >= 0 ) {
// //             int len = read.length()-lsub.length();
// //             if ( len < 0 ) { good = false; continue; }
// //             if ( overlap + len <= (int)rstr.length() ) 
// //                 rsub = rstr.substr(overlap,  read.length()-lsub.length());
// //             else good = false;
// //         }
// //         else {
// //             int len = -1*overlap;
// //             if ( lsub.length()+len > read.length() ) good = false;
// //             len = read.length()-lsub.length()+overlap;
// //             if ( len > (int)rstr.length() ) good = false;
// //             rsub = read.substr(lsub.length(), -1*overlap) + rstr.substr(0, read.length()-lsub.length()+overlap);
// //         }
        
// //         if ( !good ) continue;
// //         std::string nstr = lsub + rsub;
// //         if ( param.verbose ) {
// //             std::cout << "Read:" << read << "\n";
// //             std::cout << "nStr:" << nstr << "\n";
// //         }

// //         int npos = scoring::countPositive(read, nstr, BLOSUM62);
// //         double score = (double)npos/read.length();
// //         if ( param.verbose ) std::cout << "Score:" << score << "\n";
// //         if ( score < param.latch_score ) {
// //             if ( param.verbose ) std::cout << "Weak score\n"; 
// //             continue;
// //         }
        
// //         good_reads.push_back( *cit );
// //         read_inits.push_back( left_rids[*cit] );
// //     }
}


bool assem::connectPathPair( BitString *bstrs,
                             char *strands,
                             ReadId *pairs,
                             PathToAlnMap &path2aln_map,
                             KmerToPathMap &pathid_map,
                             PathIdSet &merged_paths,
                             PathIdSet &success_paths,
                             PathId sbjct_pid,
                             PathId query_pid,
                             PathId *used_reads,
                             ReadIdArray &good_reads,
                             std::vector<int> &read_inits,
                             int overlap,
                             int direction,
                             InvertedIndex &iindex,
                             Param &param
                             )
{
    //if ( overlap <= 0 ) return false;

    SpaPath *sbj_path = path2aln_map[sbjct_pid];
    KmerId *sbj_kmers = sbj_path->getKmers();
    size_t sbj_nkmer = sbj_path->getKmerCount();
    std::string sbj_str = biostr::getSequenceString(sbj_kmers, sbj_nkmer, param.kmer_size);
    int error = 0;
        
    SpaPath *qry_path = path2aln_map[query_pid];
    KmerId *qry_kmers = qry_path->getKmers();
    size_t qry_nkmer = qry_path->getKmerCount();
    std::string qry_str = biostr::getSequenceString(qry_kmers, qry_nkmer, param.kmer_size);

    if ( param.verbose ) {
        std::cout << "sbjct:" << sbj_str << "\n";
        std::cout << "query:" << qry_str << "\n";
    }

    dropPathIds( sbjct_pid, sbj_kmers, sbj_nkmer, pathid_map );
    dropPathIds( query_pid, qry_kmers, qry_nkmer, pathid_map );
    
    bool success = false;
    AlignPosList il, dl;

    
    if ( il.size() > 0 || dl.size() > 0 ) {
        if ( param.verbose ) std::cout << "Gappy joining\n";
    }

    if ( overlap > 0 ) {
        // pivot - sbjct (longer sequence)
        // xxxxxxxxxxxxxxxxxxxxooooo           sbjct (+)
        //                     oooooxxxxxxxxxx query (-)
        if ( direction == POSITIVE ) {
            if ( param.verbose ) std::cout << "Positive direction latching\n";
            int init = sbj_str.length()-1*overlap;
            int tgap = qry_str.length()-1*overlap;
            sbj_path->join(qry_path, init, 0, tgap, il, dl, param.kmer_size, iindex, bstrs, param);
        }
        
        // xxxxxxxxxxooooo                     query (+)
        //           oooooxxxxxxxxxxxxxxxxxxxx sbjct (-)
        else {
            if ( param.verbose ) std::cout << "Negative direction latching\n";
            int lgap = qry_str.length()-1*overlap;
            sbj_path->join(qry_path, 0, lgap, 0, il, dl, param.kmer_size, iindex, bstrs, param);
        }
    }
    else {
        if ( overlap == NOT_OVERLAP ) return false;
        if ( param.verbose ) std::cout << "MINUS OFFSET:" << overlap << "\n";
        std::string mid = std::string(-1*overlap, 'X');

        // pivot - sbjct (longer sequence)
        // xxxxxxxxxxxxxxxxxxxxxxxx           sbjct (+)
        //                            xxxxxxxxxxxxxxx query (-)
        if ( direction == POSITIVE ) {
            if ( param.verbose ) std::cout << "Positive direction latching\n";
            int init = sbj_str.length()-1*overlap;
            int tgap = qry_str.length();
            sbj_path->append(mid, param.kmer_size);
            sbj_path->join(qry_path, init, 0, tgap, il, dl, param.kmer_size, iindex, bstrs, param);
        }
        
        // xxxxxxxxxxxxxxx                     query (+)
        //                   xxxxxxxxxxxxxxxxxxxxxxxxx sbjct (-)
        else {
            if ( param.verbose ) std::cout << "Negative direction latching\n";
            int lgap = qry_str.length()-1*overlap;
            qry_path->append(mid, param.kmer_size);
            sbj_path->join(qry_path, 0, lgap, 0, il, dl, param.kmer_size, iindex, bstrs, param);
        }
    }
    
    if ( param.verbose ) std::cout << "New Seq:";
    KmerId *rqry_kmers = sbj_path->getKmers();
    size_t rqry_nkmer = sbj_path->getKmerCount();
    std::string rqry_str = biostr::getSequenceString(rqry_kmers, rqry_nkmer, param.kmer_size);
    if ( param.verbose ) std::cout << rqry_str << "\n";
    success = true;


    if (!success) {
        // reset to original mapping 
        addPathIds( sbjct_pid,  sbj_kmers, sbj_nkmer, pathid_map );
        addPathIds( query_pid, qry_kmers, qry_nkmer, pathid_map );
        return false;
    }

    addPathIds( sbjct_pid, sbj_path->getKmers(), sbj_path->getKmerCount() , pathid_map );
    
    merged_paths.insert(query_pid);
    success_paths.insert(sbjct_pid);
    
    sbj_path->addReads( &good_reads[0], &read_inits[0], good_reads.size() );
    
    error = 0;

    SpaPath *old_path = new SpaPath( *sbj_path );
    SpaPath *new_path = sbj_path;
	MSA msa(*new_path, sbjct_pid, used_reads, bstrs, strands, pairs, PILE|FLAT|TRIM|CODON, error, param);
    if ( error ) {
        if ( param.verbose ) std::cout << "MSA fail:\tcode:" << error << "\n";
        //dropPathIds( sbjct_pid, new_path->getKmers(), new_path->getKmerCount() , pathid_map );
        path2aln_map[query_pid]->resetReads(used_reads, NOT_PATH);
        //path2aln_map.erase(sbjct_pid);
        delete new_path;
        path2aln_map[sbjct_pid] = old_path;
        return false;
    }
    else delete old_path;

    if ( param.verbose ) {
        std::cout << "\nAlignment\n";
        msa.printAlignment(std::cout, *new_path, 100);
        std::cout << "\nProfile\n";
        msa.printProfile(std::cout);
    }

    return true;
}

bool assem::latch( BitString *bstrs,
            char *strands,
            ReadId *pairs,
            PathToAlnMap &path2aln_map,
            KmerToPathMap &pathid_map,
            PathIdSet &merged_paths,
            PathIdSet &success_paths,
            PathId left,
            PathId right,
            PathId *used_reads,
            ReadIdArray &good_reads,
            std::vector<int> &read_inits,
            int overlap,
            InvertedIndex &iindex,
			Param &param
            )
{
    SpaPath *lspath = path2aln_map[left];
    KmerId *lkmers = lspath->getKmers();
    size_t lnkmer = lspath->getKmerCount();
    std::string lstr = biostr::getSequenceString(lkmers, lnkmer, param.kmer_size);
    int error = 0;
        
    SpaPath *rspath = path2aln_map[right];
    KmerId *rkmers = rspath->getKmers();
    size_t rnkmer = rspath->getKmerCount();
    std::string rstr = biostr::getSequenceString(rkmers, rnkmer, param.kmer_size);

//     ReadId *rrids = rspath->getReads();
//     size_t rnum = rspath->getReadCount();
//     if ( param.verbose ) {
//         std::cout << "Right Reads\n";
//         for ( size_t i = 0; i < rnum; i++ )
//             std::cout << rrids[i] << " ";
//         std::cout << "\n";
//     }


    dropPathIds( left,  lkmers, lnkmer, pathid_map );
    dropPathIds( right, rkmers, rnkmer, pathid_map );
    
    bool success = false;
    AlignPosList il, dl;
    if ( overlap >= 0 ) {
        int init = lstr.length()-1*overlap;
        int tgap = rstr.length()-1*overlap;
        //if ( param.verbose ) std::cout << "Rstr:" << rstr.size() << "\tOverlap:" << overlap << "\tInit:" << init << "\tTgap:" << tgap << "\n";
        if ( il.size() > 0 || dl.size() > 0 ) {
            if ( param.verbose ) std::cout << "Gappy joining\n";
        }
        lspath->join(rspath, init, 0, tgap, il, dl, param.kmer_size, iindex, bstrs, param);
        if ( param.verbose ) std::cout << "New Seq:";
        KmerId *rrkmers = lspath->getKmers();
        size_t rrnkmer = lspath->getKmerCount();
        std::string rrstr = biostr::getSequenceString(rrkmers, rrnkmer, param.kmer_size);
        if ( param.verbose ) std::cout << rrstr << "\n";
        success = true;
    }
    else {
        //bool success = false;
        if ( param.verbose ) {
            std::cout << "MINUS OFFSET\n";
            std::cout << "Overlap:" << overlap << "\n";
        }
        for ( size_t i = 0; i < good_reads.size(); i++ ) {
            std::string aread = bstrs[good_reads[i]].toString();
            int s = lstr.length()-read_inits[i];
            if ( param.verbose ) std::cout << good_reads[i] << "\tsize:" << aread.size() << "\tlstr:" << lstr.length() << "\t" << read_inits[i] << "\t" << s << "\n";
            if ( s -1*overlap > (int)aread.size() ) continue;
            
            std::string mid = aread.substr(s, -1*overlap);
            if ( param.verbose ) std::cout << "MID:" << mid << "\n";
            
            int init = lstr.length()-1*overlap;
            int tgap = rstr.length();
            if ( il.size() > 0 || dl.size() > 0 ) {
                if ( param.verbose ) std::cout << "Gappy joining\n";
            }
            lspath->join(rspath, mid, init, 0,  tgap, il, dl, param.kmer_size, iindex, bstrs, param);
            success = true;
            break;
        }
        //if ( !success ) return false;
    }

    if (!success) {
        // reset to original mapping 
        addPathIds( left,  lkmers, lnkmer, pathid_map );
        addPathIds( right, rkmers, rnkmer, pathid_map );
        return false;
    }

    addPathIds( left, lspath->getKmers(), lspath->getKmerCount() , pathid_map );
    
//     size_t qnkmer = rspath->getKmerCount();
//     KmerId *qkmer = rspath->getKmers();
//     for ( size_t i = 0; i < qnkmer; i++ ) {
//         if ( pathid_map.find(qkmer[i]) == pathid_map.end() ) continue;
//         pathid_map[qkmer[i]].erase(right);
//     }

    merged_paths.insert(right);
    success_paths.insert(left);
    
//     for ( size_t i = 0; i < qnkmer; i++ ) 
//         pathid_map[qkmer[i]].insert(left);

    /* add bridging reads */
    lspath->addReads( &good_reads[0], &read_inits[0], good_reads.size() );
    //    if ( param.verbose ) {
//         std::cout << "Good reads\n";
//         for ( size_t i = 0; i < good_reads.size(); i++ )
//             std::cout << good_reads[i] << " ";
//         std::cout << "\n";
//     }
    
    error = 0;
    SpaPath *old_lpath = new SpaPath( *lspath );
    //SpaPath *old_rpath = new SpaPath( *rspath );
    SpaPath *new_lpath = lspath;
	MSA msa(*new_lpath, left, used_reads, bstrs, strands, pairs, PILE|FLAT|TRIM|CODON, error, param);
    if ( error ) {
        if ( param.verbose ) std::cout << "MSA fail:\tcode:" << error << "\n";
        //dropPathIds( left, new_lpath->getKmers(), new_lpath->getKmerCount() , pathid_map );
        //path2aln_map[left]->resetReads(used_reads, NOT_PATH);
        //path2aln_map.erase(left);
        ReadId *rreads = rspath->getReads();
        for ( size_t i = 0; i < rspath->getReadCount(); i++ )
            if ( used_reads[rreads[i]] != right ) used_reads[rreads[i]] = right;
        delete new_lpath;
        path2aln_map[left] = old_lpath;
        return false;
    }
    else delete old_lpath;

    if ( param.verbose ) {
        std::cout << "\nAlignment\n";
        msa.printAlignment(std::cout, *new_lpath, 100);
        std::cout << "\nProfile\n";
        msa.printProfile(std::cout);
    }

    return true;
}


void assem::remapLatchReads( PathId left, 
                  PathId right, 
                  PathPairMap &pair_map,
                  PathReadMap &latch_left,
                  PathReadMap &latch_right,
                  ReadIdArray &good_reads,
                  int lstr_size,
                  int &overlap )
{
    if ( overlap < 0 ) overlap = 0;

    //for ( ReadPosList::iterator it = latch_left[right].begin(); it != latch_left[right].end(); ++it ) {
    for ( ReadPosList::iterator it = latch_right[right].begin(); it != latch_right[right].end(); ++it ) {
        it->rpos += ( lstr_size-overlap );

        latch_right[left].push_back(*it);
    }

    pair_map[left] = pair_map[right];
    pair_map.erase(right);
}

bool assem::joinReadsLeftOver( PathId *used_reads,
                               BitString *bstrs,
                               ReadIdArray inspected,
                               ReadIdArray path_reads, /* sorted already */
                               PathId pid,
                               PathToAlnMap &path2aln_map,
                               Param &param)
{
    KmerId *nlkmers = path2aln_map[pid]->getKmers();
    size_t nlskmer  = path2aln_map[pid]->getKmerCount();
    std::string consensus = biostr::getSequenceString( nlkmers, nlskmer, param.kmer_size );
    
    ReadIdArray more_reads;
    std::vector<int> more_inits;
    
    std::sort(inspected.begin(), inspected.end()); 
    ReadIdArray nope_reads = Set::Difference<ReadId>(inspected, path_reads, true);
    if ( param.verbose ) std::cout << "# to be recruited:" << nope_reads.size() << "\n";
    for ( size_t i = 0; i < nope_reads.size(); i++ ) {
        if ( used_reads[nope_reads[i]] != NOT_PATH ) continue;
        std::string query = bstrs[nope_reads[i]].toString();
        
        if ( param.verbose ) std::cout << "\t" << i << "\t" << nope_reads[i] << "\t" << query << "\n";
        size_t found = consensus.find(query);            
        if ( found != std::string::npos ) {
            more_reads.push_back(nope_reads[i]);
            more_inits.push_back(found);
            used_reads[nope_reads[i]] = pid;
            if ( param.verbose ) std::cout << "\tSubstring match at " << found << "\n";
        } else {
            GlobalAlignPair paln = GlobalAlignPair(consensus, query);
            AlignSummary summary = paln.getSummary();
                
            if ( param.verbose ) {
                std::cout << paln.getAlignment();
                summary.print(std::cout);
//                 std::cout << "\t#ins:" << summary.ins.size() << "\t";
//                 std::cout << "\t#del:" << summary.del.size() << "\t";
//                 std::cout << "\t#mms:" << summary.mismatch << "\n";
//                 std::cout << "\t#out:" << summary.outer.first << "\t" << summary.outer.second << "\n";
//                 std::cout << "\t#in:" << summary.range.first << "\t" << summary.range.second << "\n";
            }
            if ( summary.lgap.first > 0 ) {
                if (param.verbose) std::cout << "Wrong alignment\n"; continue;
            }
            if ( summary.positive < param.merge_score ) {
                if (param.verbose) std::cout << "Weak alignment\n"; continue;
            }
            if ( summary.ins.size() > 0 || summary.del.size() > 0) {
                if ( param.verbose ) std::cout << "No indel here\n"; continue;
            }
            
            more_reads.push_back(nope_reads[i]);
            more_inits.push_back(summary.lgap.second);
            used_reads[nope_reads[i]] = pid;
        }
    }
    
    if ( more_reads.size() > 0 ) {
        path2aln_map[pid]->addReads( &more_reads[0], &more_inits[0], more_reads.size() );
        return true;
    }
    return false;
}

void assem::inspectLatchReads( PathId left,
                               PathId curr,
                               ReadIdArray &good_reads,
                               std::vector<int> &read_inits,
                               int &overlap, 
                               std::map<int, PathId> &commpid_map,
                               std::tr1::unordered_map<PathId, ReadIdArray> &commset_map,
                               BitString             *bstrs,
                               PathReadMap           &latch_left,
                               PathReadMap           &latch_right,
                               PathToAlnMap          &path2aln_map,
                               PathId *used_reads,                                
                               Param                 &param )
    
{
    KmerId *lkmers = path2aln_map[left]->getKmers();
    size_t lnkmers = path2aln_map[left]->getKmerCount();
    std::string lstr = biostr::getSequenceString(lkmers, lnkmers, param.kmer_size);

    KmerId *rkmers = path2aln_map[curr]->getKmers();
    size_t rnkmers = path2aln_map[curr]->getKmerCount();
    std::string rstr = biostr::getSequenceString(rkmers, rnkmers, param.kmer_size);
    if ( param.verbose ) {
        std::cout << "Curr:" << curr << "\n";
        std::cout << "Right:" << rstr << "\n";
        std::cout << "Length:" << rstr.size() << "\n";
    }
    ReadIdArray comm_reads = commset_map[curr];
    
    validLatchReads(used_reads, bstrs, lstr, rstr, left, curr, latch_left, latch_right, comm_reads, good_reads, read_inits, overlap, param );
}


/*
 paired end? check direction conflict - ignore
 */
void assem::findMaxExtendiblePath( PathId                 left,
                                   std::map<int, PathId> &commpid_map,
                                   std::tr1::unordered_map<PathId, ReadIdArray> &commset_map, 
                                   PathId                &mpath,
                                   ReadIdArray           &mrids,
                                   std::vector<int>      &mbegs,
                                   int                   &mover,
                                   PathId                *used_reads,
                                   BitString             *bstrs,
                                   PathReadMap           &latch_left,
                                   PathReadMap           &latch_right,
                                   PathToAlnMap          &path2aln_map,
                                   Param                 &param
                                   )
{
    std::map<int, PathId>::reverse_iterator pit;
    for ( pit =  commpid_map.rbegin(); pit != commpid_map.rend(); ++pit ) {
        /* already found best one */
        if ( (int)mrids.size() >= pit->first ) return;
        
        PathId curr = pit->second;
        ReadIdArray good_reads;
        std::vector<int> read_inits;
        int overlap = 0;
        inspectLatchReads(left, curr, good_reads, read_inits, overlap, commpid_map, commset_map, bstrs, latch_left, latch_right, path2aln_map, used_reads, param );
        if ( good_reads.size() > mrids.size() ) {
            mpath = curr;
            mover = overlap;
            mrids = good_reads;
            mbegs = read_inits;
        }
    }
}

void assem::getPairedReadCounts( std::multimap<int, PathId> &pcount_map,
                                 PathId left,
                                 PathIdSet &pids,
                                 ReadId *pairs,
                                 PathId *used_reads,
                                 char *strands,
                                 PathToAlnMap &path2aln_map,
                                 PathIdSet &merged_paths,
                                 Param &param )
{
    for ( PathIdSet::iterator it = pids.begin(); it != pids.end(); ++it ) {
        if ( merged_paths.find(*it) != merged_paths.end() ) continue;
        
        ReadIdArray qreads, sreads;
        setPairedReads( qreads, sreads, left, *it, used_reads, pairs, path2aln_map );
        if ( param.verbose ) std::cout << left << "\t" << *it << "\tPair reads:" << qreads.size() << "\n";
        
        int count = (int)qreads.size();

        IntPair strs = determineStrands ( sreads, pairs, strands, param );
        if ( param.verbose ) {
            std::cout << "Strands:" << strs.first << "/" << strs.second << "\n";
        }
        if ( strs.first == strs.second ) {
            if ( param.verbose ) std::cout << "** Same strands -> skip\n"; 
            continue;
        }
        if ( strs.first == '-' ) {
            if ( param.verbose ) std::cout << "** Illegal strands -> skip\n"; 
            continue;
        }
        pcount_map.insert(std::pair<int, PathId>(count, *it));
    }
}

void assem::latchPair( BitString *bstrs,
                       char *strands,
                       ReadId *pairs,
                       PathId *used_reads,
                       PathPairMap &pair_map,
                       //std::tr1::unordered_map<PathId, std::list<std::pair<PathId, size_t> > > &pairedreads_counts,
                       PathReadMap &latch_left,
                       PathReadMap &latch_right,
                       PathToAlnMap &path2aln_map,
                       KmerToPathMap &pathid_map,
                       PathIdSet &merged_paths,
                       PathId left,
                       InvertedIndex &iindex,
                       PathIdSet &successPaths,
                       PathIdSet &dirty_paths,
                       Param &param )
{
    KmerId *lkmers = path2aln_map[left]->getKmers();
    size_t lnkmers = path2aln_map[left]->getKmerCount();
    std::string lstr = biostr::getSequenceString(lkmers, lnkmers, param.kmer_size);
    
    if ( param.verbose ) {
        std::cout << "Left:" << lstr << "\n";
        std::cout << "Length:" << lstr.size() << "\n";
    }

    /* common read sets by path id */
    std::map<int, PathId> commpid_map;
    std::tr1::unordered_map<PathId, ReadIdArray> commset_map;
    commonReadSet( left, pair_map, latch_left, latch_right, commpid_map, commset_map, param );    
    if ( commpid_map.size() == 0 ) return;

    /* greedy best path to extend by by reads supports in single end
       # pairs by paired end reads
     */

    PathId           best;
    ReadIdArray      good_reads;
    int              overlap;
    std::vector<int> read_inits;
    ReadIdArray      pair_reads; // reads from left path that paired to reads in right path.

    std::multimap<int, PathId> pcount_map;
    if ( param.pair_flag ) {
        getPairedReadCounts(pcount_map, left, pair_map[left], pairs, used_reads, strands, path2aln_map, merged_paths, param );
        if ( pcount_map.size() == 0 ) {
            if ( param.verbose ) std::cout << "No paired paths found\n"; 
            return;
        }
//         if ( pairedreads_counts.find(left) == pairedreads_counts.end() ) {
//             if ( param.verbose ) std::cout << "No paired paths found\n"; 
//             return;
//         }
        
//         std::multimap<size_t, PathId> pcount_map;    
//         std::list<std::pair<PathId, size_t> >::iterator pt;
//         for ( pt = pairedreads_counts[left].begin(); pt != pairedreads_counts[left].end(); ++pt )
//             pcount_map.insert(std::pair<size_t, PathId>(pt->second, pt->first) );
        
        bool found = false;
        std::multimap<int, PathId>::reverse_iterator it = pcount_map.rbegin();
        for ( ; it != pcount_map.rend(); ++it ) {
            if ( param.verbose ) std::cout << "# pair-end reads:" << it->first << "\n";
            if ( commset_map.find(it->second) == commset_map.end() ) continue;
            if ( (int)commset_map[it->second].size() < param.latch_support ) continue;
            
            /***
             * It may be better use mininum no. of pair-end read instead of latch_support
             */
            if ( it->first < param.latch_support ) {
                if (param.verbose) std::cout << "Weak latch supports\n";
                //return;
                continue;
            }

            found = true;
            best = it->second;
            inspectLatchReads(left, best, good_reads, read_inits, overlap, commpid_map, commset_map, bstrs, latch_left, latch_right, path2aln_map, used_reads, param );    
        }
        if ( !found ) return;
    }
    else {
        findMaxExtendiblePath(left, commpid_map, commset_map, best, good_reads, read_inits, overlap, used_reads, bstrs, latch_left, latch_right, path2aln_map, param);
        if ( (int)good_reads.size() < param.latch_support ) {
            if ( param.verbose ) std::cout << "Weak supported latch:\n";
            return;
        }
    }

    bool success = latch( bstrs, strands, pairs, path2aln_map, pathid_map, merged_paths, successPaths, left, best, used_reads, good_reads, read_inits, overlap, iindex, param);
    
    if ( success ) {
        if ( param.verbose ) std::cout << "Latch Success\n";
    
        std::sort(good_reads.begin(), good_reads.end());

        /* Handle unrecruited reads */
        if ( joinReadsLeftOver( used_reads, bstrs, getReadIds(latch_left[best]),  good_reads, left, path2aln_map, param ) )
            dirty_paths.insert(left);
        if ( joinReadsLeftOver( used_reads, bstrs, getReadIds(latch_right[left]), good_reads, left, path2aln_map, param ) )
            dirty_paths.insert(left);

        // Best right path has been joined to the left.
        // So, drop it 
        latch_left.erase(best); 

        /* No more paths to be joined */
        if ( pair_map.find(best) == pair_map.end() ) {
            latch_right.erase(left);
            return;
        }

        /* remap the read ownership from right to left */
        latch_right[left].clear();
        remapLatchReads( left, best, pair_map, latch_left, latch_right, good_reads, lstr.length(), overlap );

        if ( param.verbose ) std::cout << "Recursive Latch\n";
        //latchPair( bstrs, strands, pairs, used_reads, pair_map, pairedreads_counts, latch_left, latch_right, path2aln_map, pathid_map, merged_paths, left, iindex, successPaths, dirty_paths, param );
        latchPair( bstrs, strands, pairs, used_reads, pair_map, latch_left, latch_right, path2aln_map, pathid_map, merged_paths, left, iindex, successPaths, dirty_paths, param );
        
        
//         KmerId *kmers = path2aln_map[left]->getKmers();
//         size_t  nkmer = path2aln_map[left]->getKmerCount();
//         //         for ( size_t i = 0; i < nkmer; i++ ) {
//         //             pathid_map[kmers[i]].erase(best);
//         //             pathid_map[kmers[i]].insert(left);
//         //         }
//         dropPathIds( best, kmers, nkmer, pathid_map );
//         addPathIds( left,  kmers, nkmer, pathid_map );
        
    }
}
        

//         KmerId *nlkmers = path2aln_map[left]->getKmers();
//         size_t nlskmer  = path2aln_map[left]->getKmerCount();
//         //std::string consensus = biostr::stripGap( biostr::getSequenceString( nlkmers, nlskmer, param.kmer_size ) );
//         std::string consensus = biostr::getSequenceString( nlkmers, nlskmer, param.kmer_size );
                
//         ReadIdArray more_reads;
//         std::vector<int> more_inits;

//         ReadIdArray lreads = getReadIds( latch_left[best] );
//         std::sort(lreads.begin(), lreads.end());
//         ReadIdArray lrdiff = Set::Difference<ReadId>(lreads, good_reads, true);
//         if ( param.verbose ) std::cout << "Left over left:" << lrdiff.size() << "\n";
//         for ( size_t i = 0; i < lrdiff.size(); i++ ) {
//             if ( used_reads[lrdiff[i]] != NOT_PATH ) continue;
//             std::string query = bstrs[lrdiff[i]].toString();

//             if ( param.verbose ) std::cout << "\t" << i << "\t" << lrdiff[i] << "\t" << query << "\n";
//             size_t found = consensus.find(query);            
//             if ( found != std::string::npos ) {
//                 more_reads.push_back(lrdiff[i]);
//                 more_inits.push_back(found);
//                 used_reads[lrdiff[i]] = left;
//                 if ( param.verbose ) std::cout << "\tFound:" << found << "\n";
//             } else {
//                 GlobalAlignPair paln = GlobalAlignPair(consensus, query);
//                 AlignSummary summary = paln.getSummary();
                
//                 if ( param.verbose ) {
//                     std::cout << paln.getAlignment();
                    
//                     std::cout << "\t#ins:" << summary.ins.size() << "\t";
//                     std::cout << "\t#del:" << summary.del.size() << "\t";
//                     std::cout << "\t#mms:" << summary.mismatch << "\n";
//                     std::cout << "\t#out:" << summary.outer.first << "\t" << summary.outer.second << "\n";
//                     std::cout << "\t#in:" << summary.range.first << "\t" << summary.range.second << "\n";
//                 }
//                 if ( summary.lgap.first > 0 ) {
//                     if (param.verbose) std::cout << "Wrong alignment\n"; continue;
//                 }
//                 if ( summary.positive < param.merge_score ) {
//                     if (param.verbose) std::cout << "Weak alignment\n"; continue;
//                 }
//                 if ( summary.ins.size() > 0 || summary.del.size() > 0) {
//                     if ( param.verbose ) std::cout << "No indel here\n"; continue;
//                 }

//                 more_reads.push_back(lrdiff[i]);
//                 more_inits.push_back(summary.lgap.second);
//                 used_reads[lrdiff[i]] = left;
//             }
//         }

//         path2aln_map[left]->addReads( &more_reads[0], &more_inits[0], more_reads.size() );

//         more_reads.clear();
//         more_inits.clear();

//         ReadIdArray rreads = getReadIds( latch_right[left] );
//         std::sort(rreads.begin(), rreads.end());
//         ReadIdArray rrdiff = Set::Difference<ReadId>(rreads, good_reads, true);
//         // no dups
//         rrdiff = Set::Difference<ReadId>(rrdiff, lrdiff, true);
//         if (param.verbose) std::cout << "Left over right:" << rrdiff.size() << "\n";
//         for ( size_t i = 0; i < rrdiff.size(); i++ ) {
//             if ( used_reads[rrdiff[i]] != NOT_PATH ) continue;
//             std::string query = bstrs[rrdiff[i]].toString();
//             if (param.verbose) std::cout << "\t" << i << "\t" << rrdiff[i] << "\t" << query << "\n";
//             size_t found = consensus.find(query);            
//             if ( found != std::string::npos ) {
//                 more_reads.push_back(rrdiff[i]);
//                 more_inits.push_back(found);
//                 used_reads[rrdiff[i]] = left;
//                 if (param.verbose) std::cout << "\tFound:" << found << "\n";
//             } else {
//                 GlobalAlignPair paln = GlobalAlignPair(consensus, query);
//                 AlignSummary summary = paln.getSummary();
//                 if ( param.verbose ) {
//                     std::cout << paln.getAlignment();
//                     std::cout << "\t#ins:" << summary.ins.size() << "\t";
//                     std::cout << "\t#del:" << summary.del.size() << "\t";
//                     std::cout << "\t#mms:" << summary.mismatch << "\n";
//                     std::cout << "\t#out:" << summary.outer.first << "\t" << summary.outer.second << "\n";
//                     std::cout << "\t#in:" << summary.range.first << "\t" << summary.range.second << "\n";
//                 }
//                 if ( summary.lgap.first > 0 ) {
//                     if (param.verbose) std::cout << "Wrong alignment\n"; continue;
//                 }
//                 if ( summary.positive < param.merge_score ) {
//                     if (param.verbose) std::cout << "Weak alignment\n"; continue;
//                 }
//                 if ( summary.ins.size() > 0 || summary.del.size() > 0) {
//                     if (param.verbose) std::cout << "No indel here\n"; continue;
//                 }

//                 more_reads.push_back(rrdiff[i]);
//                 more_inits.push_back(summary.lgap.second);
//                 used_reads[rrdiff[i]] = left;
//             }
//         }



//         path2aln_map[left]->addReads( &more_reads[0], &more_inits[0], more_reads.size() );

        
//         latch_left.erase(best);
//         if ( pair_map.find(best) == pair_map.end() ) {
//             latch_right.erase(left);
//             return;
//         }
        
//         latch_right[left].clear();


//         KmerId *kmers = path2aln_map[left]->getKmers();
//         size_t  nkmer = path2aln_map[left]->getKmerCount();
// //         for ( size_t i = 0; i < nkmer; i++ ) {
// //             pathid_map[kmers[i]].erase(best);
// //             pathid_map[kmers[i]].insert(left);
// //         }
//         dropPathIds( best, kmers, nkmer, pathid_map );
//         addPathIds( left,  kmers, nkmer, pathid_map );

//         if ( param.verbose ) std::cout << "Recursive Latch\n";
//         remapPathId( left, best, pair_map, latch_left, latch_right, good_reads, lstr.length(), overlap );
//         latchPair( bstrs, strands, pairs, used_reads, pair_map, latch_left, latch_right, path2aln_map, pathid_map, merged_paths, left, iindex, successPaths, param );
//     }
// }

void assem::latchPairs(std::list<LatchRead> &latch_reads,
                       PathId *used_reads, 
                       PathToAlnMap &path2aln_map,
                       BitString* bstrs,
                       char *strands,
                       ReadId *pairs,
                       KmerToPathMap &pathid_map,
                       PathIdSet &merged_paths,
                       PathIdSet &successPaths,
                       PathIdSet &dirty_paths,
                       InvertedIndex &iindex,
                       PathReadMap &latch_left, 
                       PathReadMap &latch_right,
                       PathPairMap &pair_map,
                       //std::tr1::unordered_map<PathId, std::list<std::pair<PathId, size_t> > > &pairedreads_counts,
                       Param &param )
{
    if ( latch_right.size() == 0 ) return;

    /* Ordeded by # of latchable reads to left paths */
    std::multimap<int, PathId> rmap = util::sortBySize<PathId, ReadPosList>(latch_right);
    std::multimap<int, PathId>::reverse_iterator ct = rmap.rbegin();
    for ( ; ct != rmap.rend(); ++ct ) {
        PathId left = ct->second;
        if ( merged_paths.find(left) != merged_paths.end() ) continue;

        if ( param.verbose ) std::cout << "\nPath:" << ct->second << "\tSize:" << ct->first << "\n";
        if ( pair_map.find(left) == pair_map.end() ) {
            KmerId *lkmers = path2aln_map[left]->getKmers();
            size_t lnkmers = path2aln_map[left]->getKmerCount();
            std::string lstr = biostr::getSequenceString(lkmers, lnkmers, param.kmer_size);
            if ( param.verbose ) {
                std::cout << "No matching pair\n"; 
                std::cout << "Lstr:" << lstr << "\n";
            }
            continue; 
        }
        //latchPair( bstrs, strands, pairs, used_reads, pair_map, pairedreads_counts, latch_left, latch_right, path2aln_map, pathid_map, merged_paths, left, iindex, successPaths, dirty_paths, param );
        latchPair( bstrs, strands, pairs, used_reads, pair_map, latch_left, latch_right, path2aln_map, pathid_map, merged_paths, left, iindex, successPaths, dirty_paths, param );
    }
}


bool assem::recruitToPath( ReadId rid,
                    PathId path,
                    PathId *used_reads, 
                    BitString* bstrs,
                    PathToAlnMap &path2aln_map,
                    KmerToPathMap &pathid_map,
                    PathIdSet &merged_paths,
                    PathIdSet &successPaths,
					Param &param
                    )
{
    std::string rstr = bstrs[rid].toString();
    if ( param.verbose ) std::cout << "Read:" << rid << "\t" << rstr << "\n";

    KmerArray kmers = biostr::getKmers( rstr, param.kmer_size );
    //int mink = getMinKmerCount( kmers.size()+param.kmer_size-1, param );
    //int mink = filter::minSameKmerCount( param.latch_length, param.kmer_size, 1-param.latch_score );
    //if ( mink < 1 ) mink = 1;
    int mink = 1;

    PosPathPairList pos_paths = searchSimilarPaths(PathIdSet(), &kmers[0], kmers.size(), path2aln_map, pathid_map, merged_paths, mink, true, param);
    if ( pos_paths.size() == 0 ) {
        if ( param.verbose ) std::cout << "0 similar path\n";
        return false;
    }
    std::vector<int> pos_vec;
    PosPathPairList::iterator pt;
    for ( pt = pos_paths.begin(); pt != pos_paths.end(); ++pt ) {
        if ( pt->second == path ) pos_vec = pt->first;
    }
    if ( pos_vec.size() == 0 ) {
        if ( param.verbose ) std::cout << "No similar kmers to the path\n";
        return false;
    }
    
    IntPair rrange = getMatchRange(pos_vec, param);
    KmerId *rkmers = path2aln_map[path]->getKmers();
    size_t rnkmers = path2aln_map[path]->getKmerCount();

    
    if ( rrange.first == -1 || rrange.second == -1 ) return false;
    IntPair qrange; 
    qrange.first = pos_vec[rrange.first];
    qrange.second = pos_vec[rrange.second];
    
    //int type = getLinkType( qrange, rrange, &kmers[0], rkmers, kmers.size(), rnkmers, param.latch_length, param.read_spur, param );
    int type = getPathPairType( pos_vec, qrange, rrange, &kmers[0], rkmers, kmers.size(), rnkmers, param.latch_length, param.read_spur, param );
    
    if ( type == NOTYPE ) return false;
    
    if ( !ordered(pos_vec, rrange) ) { if ( param.verbose ) std::cout << "Not in order\n"; return false; }
    if ( indeled(pos_vec, rrange) ) { if ( param.verbose ) std::cout << "Indeled\n"; return false; }
    
    
    int qnkmers = kmers.size();
    if ( type == SPUR || type == WEAK_SPUR || type == ROPE || type == WEAK_ROPE ) 
        adjustRanges(qrange, rrange, qnkmers, rnkmers, param);
    else if ( type == LATCH ) adjustLatchRange(qrange, rrange, qnkmers, rnkmers, param );

    std::string query = biostr::getSequenceString( &kmers[qrange.first], qrange.second-qrange.first+1, param.kmer_size);
    std::string sbjct = biostr::getSequenceString(&rkmers[rrange.first], rrange.second-rrange.first+1, param.kmer_size);
    
    if ( param.verbose ) {
        std::cout << "Query:" << query << "\n";
        std::cout << "Sbjct:" << sbjct << "\n";
    }
    
    AlignSummary summary = compareBySubstitution(query, sbjct, qnkmers, rnkmers, qrange.first, qrange.second, rrange.first, rrange.second, param);
    if ( type == LATCH ) {
        if ( summary.posrate < param.latch_score ) { if ( param.verbose ) std::cout << "Weak score\n"; return false; }
    } else {
        if ( summary.posrate < param.merge_score ) { if ( param.verbose ) std::cout << "Weak score\n"; return false; }
    }

    int aln_pos = rrange.first;
    if ( type != LATCH ) {
        if ( type != LATCH && qrange.first > 0 ) aln_pos = ( rrange.first - qrange.first );
    } else {
        if ( (double)query.size() < param.latch_ratio*rstr.size() ) {
            if ( param.verbose ) std::cout << "Short read latch\n"; 
            return false;
        }
        if ( qrange.first > 0 ) aln_pos = -1*qrange.first;
    }

    SpaPath *spath = path2aln_map[path];
    spath->addRead( rid, aln_pos );
    used_reads[rid] = path;

    return true;
}

int assem::recruitReadsToPath( ReadIdSet &rset,
                        PathId path,
                        PathId *used_reads, 
                        BitString* bstrs,
                        PathToAlnMap &path2aln_map,
                        KmerToPathMap &pathid_map,
                        PathIdSet &merged_paths,
                        PathIdSet &successPaths,
						Param &param
                        )
{
    int count = 0;
    
    for ( ReadIdSet::iterator it = rset.begin(); it != rset.end(); ++it ) {
        if ( used_reads[*it] != NOT_PATH ) 
            continue;

        bool succ = recruitToPath(*it, path, used_reads, bstrs, path2aln_map, pathid_map, merged_paths, successPaths, param);
        if ( succ ) count++;
    }
    return count;
}

void assem::latchReadsToPath( PathReadMap &latch_reads,
                       PathId *used_reads, 
                       char *strands,
                       ReadId *pairs, 
                       BitString* bstrs,
                       PathToAlnMap &path2aln_map,
                       KmerToPathMap &pathid_map,
                       PathIdSet &merged_paths,
                       PathIdSet &successPaths,
                       int direction,
					   Param &param )
{
    if ( param.verbose ) {
        if ( direction == RIGHT ) 
            std::cout << "Latch reads to right\n";
        else 
            std::cout << "Latch reads to left\n";
        std::cout << latch_reads.size() << "\n";
    }
    for ( PathReadMap::iterator it = latch_reads.begin(); it != latch_reads.end(); ++it ) {
        PathId path = it->first;
        if ( merged_paths.find(path) != merged_paths.end() ) continue;
        if ( param.verbose ) std::cout << "PathId:" << path << "\n";

        ReadIdSet rset;
        for ( ReadPosList::iterator rt = (it->second).begin(); rt != (it->second).end(); ++rt )
            rset.insert(rt->rid);

        ReadIdArray old_reads = ReadIdArray( path2aln_map[path]->getReads(), path2aln_map[path]->getReads() + path2aln_map[path]->getReadCount() );

        int count = recruitReadsToPath(rset, path, used_reads, bstrs, path2aln_map, pathid_map, merged_paths, successPaths, param);

        if (param.verbose) std::cout << "# Reads recruited:" << count << "\n";
        if ( count == 0 ) continue;

        int error = 0;
        SpaPath *old_path = new SpaPath( *path2aln_map[path] );
        SpaPath *new_path = path2aln_map[path];
		MSA msa(*new_path, path, used_reads, bstrs, strands, pairs, PILE|FLAT|TRIM|CODON, error, param);
        if ( error ) {
            if ( param.verbose ) std::cout << "MSA fail:\tcode:" << error << "\n";
            ReadIdArray new_reads = ReadIdArray( new_path->getReads(), new_path->getReads() + new_path->getReadCount() );            
            ReadIdArray invalids = Set::Difference<ReadId>( new_reads, old_reads, false );
            for ( size_t i = 0; i < invalids.size(); i++ ) used_reads[invalids[i]] = NOT_PATH;
            //new_path->resetReads(used_reads, NOT_PATH);
            //path2aln_map.erase(path);
            delete new_path;
            path2aln_map[path] = old_path;
            continue;
        } else delete old_path;

        if ( param.verbose ) {
            std::cout << "\nAlignment\n";
            msa.printAlignment(std::cout, *new_path, 100);
            std::cout << "\nconsensus:";
            std::cout << msa.getConsensus() << "\n"; 
        }

        KmerId *kmers = new_path->getKmers();
        size_t  nkmer = new_path->getKmerCount();
        for ( size_t i = 0; i < nkmer; i++ ) {
            pathid_map[kmers[i]].insert(path);
        }
        successPaths.insert(path);
    }
}

void assem::inspectPathPairs( PathPairMap &path_pairs,
                              //std::tr1::unordered_map<PathId, std::list<std::pair<PathId, size_t> > > &pairedreads_counts,
                              PathToAlnMap &path2aln_map,
                              PathId *used_reads,
                              char *strands,
                              ReadId *pairs, 
                              Param &param
                              )
{


    PathPairMap::iterator it;
    for ( it = path_pairs.begin(); it != path_pairs.end(); ++it ) {
        PathId lpath = it->first;
        ReadIdArray lreads = ReadIdArray( path2aln_map[lpath]->getKmers(), path2aln_map[lpath]->getKmers() + path2aln_map[lpath]->getKmerCount() );
        for ( PathIdSet::iterator jt = it->second.begin(); jt != it->second.end(); ) {
            PathId rpath = *jt;
            ReadIdArray qreads = lreads;
            ReadIdArray sreads;
            if ( param.verbose ) std::cout << "Path pair:" << lpath << "\t" << rpath << "\n";
            setPairedReads( qreads, sreads, lpath, rpath, used_reads, pairs, path2aln_map );

            if ( qreads.size() == 0 ) {
                if ( param.verbose ) std::cout << "Zero supported latch\n";
            }

            IntPair strs = determineStrands ( sreads, pairs, strands, param );            
            if ( qreads.size() && strs.first == strs.second ) {
                if ( param.verbose ) 
                    std::cout << "Strand conflict:" << strs.first << "\t" << strs.second << "\n";
                it->second.erase(jt++);
            } 
            else {
                //pairedreads_counts[lpath].push_back( std::pair<PathId, size_t>( *jt, qreads.size() ) );
                ++jt;
            }
        }
    }

    for ( it = path_pairs.begin(); it != path_pairs.end(); ) {
        if ( it->second.size() == 0 ) path_pairs.erase(it++);
        else ++it;
    }
}

void assem::connectPairedPathsByBridgingReads( PathId *used_reads, 
                                              PathToAlnMap &path2aln_map,
                                              BitString* bstrs,
                                              char *strands,
                                              ReadId *pairs,
                                              InvertedIndex &iindex,
                                              Param &param  )
{
    double t0 = mytime();
    
    setbuf(stdout, NULL); // no buffering
    //if ( param.verbose ) 
    std::cout << "\nConnecting paths by bridging reads ...\n";
    if ( param.verbose ) std::cout << "#Total paths:" << path2aln_map.size() << "\n";

    PathIdSet merged_paths, successPaths, dirty_paths;

    KmerToPathMap pathid_map;
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) 
        addPathIds( it->first, path2aln_map[it->first]->getKmers(), path2aln_map[it->first]->getKmerCount(), pathid_map );

    Map<PathIdPair, bool>::Type scanned_pairs;
    
    /* auxilary variables for progress log */
    double pratio = 1.0;
    int    nproc  = 0;
    double lt0    = mytime();
    size_t npath  = path2aln_map.size();

    PathLengthMap plen_map = getPathLengths( path2aln_map );
    for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {

        /* Print progress status at every percent of path proceeded */
        nproc++;
        if ( nproc/(double)npath*100 >= pratio ) {
            fprintf( stdout, "%d/%d: %.2f sec\n", nproc, int(npath), mytime()-lt0);
            pratio += 1;
        }


        PathId sbjct_pid = it->second;

        if ( merged_paths.find(sbjct_pid) != merged_paths.end() ) continue;

        /* Recursive latching */
        while(true) {
            if ( param.verbose) 
                std::cout << "\nSbjct:" << sbjct_pid << "\t" << path2aln_map[sbjct_pid]->getConsensusString() << "\n";
            
            PathIdArray pids = findPairedPaths( sbjct_pid, used_reads, pairs, path2aln_map[sbjct_pid]->getReads(), path2aln_map[sbjct_pid]->getReadCount(), param );
            if ( pids.size() == 0 ) break;

            PathIdList reduced;
            for ( PathIdArray::iterator pt = pids.begin(); pt != pids.end(); ++pt ) {
                PathIdPair ppair = sbjct_pid > *pt ? PathIdPair(*pt, sbjct_pid) : PathIdPair(sbjct_pid, *pt);
                if ( scanned_pairs.find( ppair ) == scanned_pairs.end() && 
                     merged_paths.find( *pt ) == merged_paths.end() )
                    reduced.push_back( *pt );
            }
            
            pids = PathIdArray( reduced.begin(), reduced.end() );
            if ( pids.size() == 0 ) break;
            
            std::list<PairedPath> ppaths = makePairedPaths( pids, sbjct_pid, pairs, used_reads, strands, path2aln_map, param );            
            if ( ppaths.size() == 0 ) break;

            
            std::multimap<int, PairedPath> orders;
            orderPathsBySupport( orders, ppaths, param );

            bool success = false;
            for ( std::multimap<int, PairedPath>::reverse_iterator jt = orders.rbegin(); jt != orders.rend(); ++jt ) {
                PairedPath apath = jt->second;
                if ( param.verbose ) std::cout << "Query:" << apath.pid << "\n";

                 PathIdPair ppair = sbjct_pid > apath.pid ? PathIdPair(apath.pid, sbjct_pid) : PathIdPair(sbjct_pid, apath.pid);
//                 if ( scanned_pairs.find( ppair ) != scanned_pairs.end() ) {
//                     if ( param.verbose ) std::cout << "Scanned early - skip\n";
//                     continue;
//                 }
                if ( (int)apath.match_reads.size() < param.latch_support ) {
                    if ( param.verbose ) std::cout << "Weak support\n"; 
                    scanned_pairs.insert( std::pair<PathIdPair, bool>( ppair, true) );
                    //avoids.insert( apath.pid );
                    //break;
                    continue;
                }

                success = latchPathsWithReads( sbjct_pid, apath, pathid_map, path2aln_map, bstrs, strands, pairs, used_reads, merged_paths, successPaths, iindex, param.kmer_size, param );
                if ( success ) {
                    scanned_pairs.insert( std::pair<PathIdPair, bool>( ppair, true) );
                    //avoids.insert( apath.pid );
                    //dropPathIds( apath.pid, qkmers, qnkmer, pathid_map );
                    break;
                }
            }
            if ( !success ) break;
        }
    }
    dropMergedPaths(merged_paths, path2aln_map);

    std::cout << "# Merged paths:" << merged_paths.size() << "\n";
    std::cout << "Path count:" << path2aln_map.size() << "\n";
    std::cout << "Conneting paths with reads:" << mytime()-t0 << " sec\n";

    if ( param.verbose ) std::cout << "Path connected by bridging reads:" << mytime()-t0 << " sec\n";
}

void assem::extendPathsByBridgingReads( std::list<LatchRead> &latch_reads,
                                        PathId *used_reads, 
                                        PathToAlnMap &path2aln_map,
                                        BitString* bstrs,
                                        char *strands,
                                        ReadId *pairs,
                                        KmerToPathMap &pathid_map,
                                        PathIdSet &merged_paths,
                                        PathIdSet &successPaths,
                                        PathIdSet &dirty_paths,
                                        InvertedIndex &iindex,
                                        Param &param  )
{
    double time1 = mytime();
    
    //if ( param.verbose ) 
        std::cout << "Connecting paths by bridging reads ...\n";

    //std::tr1::unordered_map<PathId, std::list<std::pair<PathId, size_t> > > pairedreads_counts;
    PathReadMap latch_left, latch_right; // right paths(read latched to the left), left paths(read latched to the right)
    PathPairMap pair_map;
    PathCountMap count_map;

    pairUpPaths( pair_map, latch_left, latch_right, latch_reads, param);
    //if ( param.pair_flag ) inspectPathPairs(pair_map, pairedreads_counts, path2aln_map, used_reads, strands, pairs, param );
    //latchPairs(latch_reads, used_reads, path2aln_map, bstrs, strands, pairs, pathid_map, merged_paths, successPaths, dirty_paths, iindex, latch_left, latch_right, pair_map, pairedreads_counts, param);
    latchPairs(latch_reads, used_reads, path2aln_map, bstrs, strands, pairs, pathid_map, merged_paths, successPaths, dirty_paths, iindex, latch_left, latch_right, pair_map, param);
    
    if ( param.extend_read_flag ) {
        latchReadsToPath(latch_left,  used_reads, strands, pairs, bstrs, path2aln_map, pathid_map, merged_paths, successPaths, LEFT, param );
        latchReadsToPath(latch_right, used_reads, strands, pairs, bstrs, path2aln_map, pathid_map, merged_paths, successPaths, RIGHT, param );
    }

    //if ( param.verbose ) {
    std::cout << "Paths connected by bridging reads:" << mytime()-time1 << " sec\n";
    //}
}


PathIdArray assem::findPairedPaths( PathId query_path,
                                    PathId *used_reads,
                                    ReadId *pairs, 
                                    ReadId *reads, 
                                    size_t nreads,
                                    Param &param )
{
    if ( nreads == 0 ) return PathIdArray();

    PathIdSet pids;
    for ( size_t i = 0; i < nreads; i++ ) {
        ReadId pair_read = pairs[reads[i]];
        if ( pair_read == NOT_PAIR ) continue;
      
        PathId pid = used_reads[pair_read];
        if ( pid == NOT_PATH || pid == query_path ) continue;
        
        pids.insert(pid);
    }
    if (param.verbose) std::cout << "# paired paths:" << pids.size() << "\n";

    return PathIdArray( pids.begin(), pids.end() );
}

ReadIdArray assem::getPairedReads( ReadIdArray &this_reads,
                                   PathId &this_path,
                                   PathId &that_path,
                                   PathId *used_reads,
                                   ReadId *pairs, 
                                   PathToAlnMap &path2aln_map )
{
    std::tr1::unordered_map<ReadId, bool> read_map;
    for ( size_t i = 0; i < this_reads.size(); i++ )
        read_map.insert( std::pair<ReadId, bool>( this_reads[i], true) );

    ReadIdList preads;
    ReadId *reads = path2aln_map[that_path]->getReads();
    size_t  nread = path2aln_map[that_path]->getReadCount();
    for ( size_t i = 0; i < nread; i++ ) {
        ReadId pair_read = pairs[reads[i]];
        if ( pair_read == NOT_PAIR ) continue;
        if ( used_reads[pair_read] != this_path ) continue;
        if ( read_map.find(pair_read) == read_map.end() ) continue;
        preads.push_back( reads[i] );
    }
    return ReadIdArray( preads.begin(), preads.end() );
}

// ReadIdArray assem::getPairedReads( PathId &this_path,
//                                    PathId &that_path,
//                                    PathId *used_reads,
//                                    ReadId *pairs, 
//                                    PathToAlnMap &path2aln_map )
// {
//     ReadIdList preads;

//     ReadId *reads = path2aln_map[that_path]->getReads();
//     size_t  nread = path2aln_map[that_path]->getReadCount();

//     for ( size_t i = 0; i < nread; i++ ) {
//         ReadId pair_read = pairs[reads[i]];
//         if ( pair_read == NOT_PAIR ) continue;
//         PathId pid = used_reads[pair_read];
        
//         if ( pid == this_path ) {
//             preads.push_back( reads[i] );
//         }
//     }
//     return ReadIdArray( preads.begin(), preads.end() );
// }

void assem::setPairedReads ( ReadIdArray &qreads,
                             ReadIdArray &sreads,
                             PathId query_path,
                             PathId sbjct_path,
                             PathId *used_reads,
                             ReadId *pairs,
                             PathToAlnMap &path2aln_map )
{
    ReadIdList qlist, slist;

    ReadId *reads = path2aln_map[sbjct_path]->getReads();
    size_t  nread = path2aln_map[sbjct_path]->getReadCount();
    for ( size_t i = 0; i < nread; i++ ) {
        ReadId pair_read = pairs[reads[i]];
        if ( pair_read == NOT_PAIR ) continue;
        PathId pid = used_reads[pair_read];
        
        if ( pid == query_path ) {
            qlist.push_back( pair_read );
            slist.push_back( reads[i] );
        }
    }
    qreads = ReadIdArray( qlist.begin(), qlist.end() );
    sreads = ReadIdArray( slist.begin(), slist.end() );
}

/**
 * Determine strands of sequence pair.
 * Each strand is determined by majority vote.
 */
IntPair assem::determineStrands( ReadIdArray &preads, // subject reads
                                 ReadId *pairs,
                                 char *strands, Param &param )
{
    IntPair strs;
    int qpos, qneg, spos, sneg;
    qpos = qneg = spos = sneg = 0;
    for ( size_t i = 0; i < preads.size(); i++ ) {
        char sstr = strands[preads[i]];
        char qstr = strands[pairs[preads[i]]];
        
        sstr == '+' ? spos++ : sneg++;
        qstr == '+' ? qpos++ : qneg++;
    }
    if ( param.verbose ) std::cout << "Strands:\tquery+:" << qpos << " query-:" << qneg << "\tsbjct+:" << spos << " sbjct-:" << sneg << "\n";
    
    qpos > qneg ? strs.first  = POSITIVE : strs.first  = NEGATIVE;
	spos > sneg ? strs.second = POSITIVE : strs.second = NEGATIVE;
	return strs;
}

/**
 * If strands of two paths are identical, report invalid strands
 * Otherwise, return the strand of a reference sequence.
 */
int assem::getDirection(ReadIdArray &preads,
                        ReadId *pairs,
                        char *strands, 
                        Param &param )
{
    IntPair strs = determineStrands ( preads, pairs, strands, param );
    if ( strs.first == strs.second ) return NONSENSE;
    //return strs.first;
    return strs.second;
}

// IntPair assem::getSimpleLatchLeftRange( std::vector<int>& pvec, int qnkmer )
// {
//     int s = 0;
//     if ( pvec[s] == -1 ) return IntPair(-1,-1);

//     int prev = pvec[s];
//     int curr = pvec[s];
//     int i = s+1;
//     for ( ; i < (int)pvec.size(); i++ ) {
//         curr = pvec[i];
//         if ( curr != prev+1 ) break;
//         prev = curr;
//     }
//     int e = i-1;
//     if ( pvec[e] != qnkmer-1 ) return IntPair(-1,-1);

//     return IntPair(s,e);
// }

// IntPair assem::getSimpleLatchRightRange( std::vector<int>& pvec )
// {
//     int e = pvec.size()-1;
//     if ( pvec[e] == -1 ) return IntPair(-1,-1);

//     int prev = pvec[e];
//     int curr = pvec[e];
//     int i = e-1;
//     for ( ; i >= 0; i-- ) {
//         curr = pvec[i];
//         if ( curr != prev-1 ) break;
//         prev = curr;
//     }
//     int s = i+1;
//     if ( pvec[s] != 0 ) return IntPair(-1, -1); 
//     return IntPair(s,e);
// }



// IntPair assem::checkLongOverlap(PathId qpid, PathId spid, PathToAlnMap &path2aln_map, int min_length, int direction, Param &param  )
// {
//     KmerId *qkmers = path2aln_map[qpid]->getKmers();
//     KmerId *skmers = path2aln_map[spid]->getKmers();
//     size_t  qnkmer = path2aln_map[qpid]->getKmerCount();
//     size_t  snkmer = path2aln_map[spid]->getKmerCount();

//     std::vector<int> pvec = makePosFlags(qnkmer, snkmer, qkmers, skmers, min_length, param);

    
//     intPair srange, qrange;
//     if ( direction == POSITIVE )
//         srange = getSimpleLatchRightRange(pvec);
//     else
//         srange = getSimpleLatchLeftRange(pvec, qnkmer);

//     if ( srange.first != -1 && srange.second != -1 && param.kmer_size + (srange.second - srange.first) >= min_length ) {
//         qrange.first  = pvec[srange.first];
//         qrange.second = pvec[srange.second];        
 
//         std::string query = biostr::getSequenceString(&qkmers[qrange.first], qrange.second-qrange.first+1, param.kmer_size);
//         std::string sbjct = biostr::getSequenceString(&skmers[srange.first], srange.second-srange.first+1, param.kmer_size);

//         if ( param.verbose) std::cout << "Simple Latch Long overlap:" << query.length() << "\n";

//         int npos = scoring::countPositive(query, sbjct, BLOSUM62);
//         double score = (double)npos/query.length();
//         if ( param.verbose ) std::cout << "Score:" << score << "\n";
        
//         if ( score == 1 ) {
//             if ( srange.first == 0 ) 
//                 return IntPair( LEFT, query.length() );
//             else                 
//                 return IntPair( RIGHT, query.length() );
//         }

//     } 

//     srange = getMatchRange(pvec, param);
//     std::string region = biostr::getSequenceString(&skmers[srange.first], srange.second-srange.first+1, param.kmer_size);

//     qrange.first  = pvec[srange.first];
//     qrange.second = pvec[srange.second];
    
//     //int type = getLinkType( qrange, srange, qkmers, skmers, qnkmer, snkmer, param.kmer_size, param.read_spur, param );
//     int type = getPathPairType( pvec, qrange, srange, qkmers, skmers, qnkmer, snkmer, min_length, param.read_spur, param );
//     if ( type == LATCH ) {
//         adjustLatchRange(qrange, srange, qnkmer, snkmer, param );
//         std::string query = biostr::getSequenceString(&qkmers[qrange.first], qrange.second-qrange.first+1, param.kmer_size);
//         std::string sbjct = biostr::getSequenceString(&skmers[srange.first], srange.second-srange.first+1, param.kmer_size);

//         /* Short */
//         if ( query.length() < min_length ) {
//             if ( param.verbose) std::cout << "Short:" << query.length() << "\n";
//             return IntPair( -1, NOT_OVERLAP );
//         }
//         else if ( param.verbose) std::cout << "Long overlap:" << query.length() << "\n";

//         int npos = scoring::countPositive(query, sbjct, BLOSUM62);
//         double score = (double)npos/query.length();
//         if ( param.verbose ) std::cout << "Score:" << score << "\n";
        
//         if ( score == 1 ) {
//             if ( srange.first == 0 ) 
//                 return IntPair( LEFT, query.length() );
//             else                 
//                 return IntPair( RIGHT, query.length() );
//         }
//     } else {
//         //if ( param.verbose ) std::cout << "\tNot a latch\n";
//     }
//     return IntPair( -1, NOT_OVERLAP );
// }



IntPair assem::checkLongOverlap(PathId qpid, PathId spid, PathToAlnMap &path2aln_map, int min_length, int direction, Param &param  )
{
    KmerId *qkmers = path2aln_map[qpid]->getKmers();
    KmerId *skmers = path2aln_map[spid]->getKmers();
    size_t  qnkmer = path2aln_map[qpid]->getKmerCount();
    size_t  snkmer = path2aln_map[spid]->getKmerCount();

    std::vector<int> pvec = makePosFlags(qnkmer, snkmer, qkmers, skmers, min_length, param);
    IntPair srange = getMatchRange(pvec, param);
    if ( srange.first == -1 || srange.second == -1 ) return IntPair( -1, NOT_OVERLAP );

    std::string region = biostr::getSequenceString(&skmers[srange.first], srange.second-srange.first+1, param.kmer_size);
//     if ( param.verbose ) {
//         std::cout << "\tLong overlap\n";
//     }

    IntPair qrange;
    qrange.first  = pvec[srange.first];
    qrange.second = pvec[srange.second];
    
//     if ( param.verbose ) 
//         std::cout << "\trrange:" << srange.first << ":" << srange.second 
//                   << "\tqrange:" << qrange.first << ":" << qrange.second << "\n";
    
    //int type = getLinkType( qrange, srange, qkmers, skmers, qnkmer, snkmer, param.kmer_size, param.read_spur, param );
    int type = getPathPairType( pvec, qrange, srange, qkmers, skmers, qnkmer, snkmer, min_length, param.read_spur, param );
    if ( type == LATCH ) {
        //adjustLatchRange(qrange, srange, qnkmer, snkmer, param );
        adjustRanges(qrange, srange, qnkmer, snkmer, param );
        std::string query = biostr::getSequenceString(&qkmers[qrange.first], qrange.second-qrange.first+1, param.kmer_size);
        std::string sbjct = biostr::getSequenceString(&skmers[srange.first], srange.second-srange.first+1, param.kmer_size);

        /* Short */
        if ( (int)query.length() < min_length ) {
            if ( param.verbose) std::cout << "Short:" << query.length() << "\n";
            return IntPair( -1, NOT_OVERLAP );
        }
        else if ( param.verbose) std::cout << "Long overlap:" << query.length() << "\n";

        int npos = scoring::countPositive(query, sbjct, BLOSUM62);
        double score = (double)npos/query.length();
        if ( param.verbose ) std::cout << "Score:" << score << "\n";
        
//         if ( score == 1 ) {
//             if ( srange.first == 0 ) 
//                 return IntPair( LEFT, query.length() );
//             else                 
//                 return IntPair( RIGHT, query.length() );
//         }
        if ( score >= param.latch_score ) {
            if ( srange.first > qrange.first ) 
                return IntPair( RIGHT, query.length() );
            else
                return IntPair( LEFT, query.length() );
                
//             if (direction == NEGATIVE)
//                 return IntPair( LEFT, query.length() );
//             else return IntPair( RIGHT, query.length() );
        }
    } else {
        //if ( param.verbose ) std::cout << "\tNot a latch\n";
    }
    return IntPair( -1, NOT_OVERLAP );
}

// IntPair assem::checkShortOverlap( PathId qpid, PathId spid, 
//                            PathToAlnMap &path2aln_map,
//                            int min_len,
// 						   Param &param )
// {
//     KmerId *qkmers = path2aln_map[qpid]->getKmers();
//     KmerId *skmers = path2aln_map[spid]->getKmers();
//     size_t  qnkmer = path2aln_map[qpid]->getKmerCount();
//     size_t  snkmer = path2aln_map[spid]->getKmerCount();
    
//     std::string query = biostr::stripGap( biostr::getSequenceString(qkmers, qnkmer, param.kmer_size) );
//     std::string sbjct = biostr::stripGap( biostr::getSequenceString(skmers, snkmer, param.kmer_size) );
//     for ( int l = param.kmer_size-1; l > min_len; l-- ) {
//         if ( query.substr( query.size()-l, l ) == sbjct.substr(0, l) ) 
//             return IntPair( LEFT, l);
//         if ( query.substr( 0, l ) == sbjct.substr(sbjct.size()-l, l) )
//             return IntPair( RIGHT, l);
//     }
//     return IntPair( -1, NOT_OVERLAP );
// }

IntPair assem::checkShortOverlap( PathId qpid, PathId spid, 
                                  PathToAlnMap &path2aln_map,
                                  int direction, 
                                  Param &param )
{
    KmerId *qkmers = path2aln_map[qpid]->getKmers();
    KmerId *skmers = path2aln_map[spid]->getKmers();
    size_t  qnkmer = path2aln_map[qpid]->getKmerCount();
    size_t  snkmer = path2aln_map[spid]->getKmerCount();
    
//     std::string query = biostr::stripGap( biostr::getSequenceString(qkmers, qnkmer, param.kmer_size) );
//     std::string sbjct = biostr::stripGap( biostr::getSequenceString(skmers, snkmer, param.kmer_size) );
    std::string query = biostr::getSequenceString(qkmers, qnkmer, param.kmer_size);
    std::string sbjct = biostr::getSequenceString(skmers, snkmer, param.kmer_size);
    
    if ( direction == NEGATIVE ) 
        util::swap<std::string>(query, sbjct);

    //for ( int l = 2*param.kmer_size; l >= min_len; l-- ) {
    //int l = 2*param.kmer_size;
    std::vector<double> sizes;
    sizes.push_back(2*param.kmer_size);
    sizes.push_back(query.size());
    sizes.push_back(sbjct.size());
    int min = (int)math::min(&sizes[0], 3);
    for ( int l = min-1; l >= param.pairend_overlap; l-- ) {
        std::string ssub = sbjct.substr( sbjct.size()-l, l );
        std::string qsub = query.substr(0, l);
        int count = scoring::countPositive(qsub, ssub, BLOSUM62);    

        if ( (double)count/qsub.length() >= param.latch_score ) {
            if (direction == NEGATIVE)
                return IntPair( LEFT, l );
            else 
                return IntPair( RIGHT, l );
        }
//         //std::cout << qsub << "\t" << ssub << "\t" << count << "\n";
//         if ( count == l ) return IntPair( LEFT, l );

//         qsub = query.substr(0, l);
//         ssub = sbjct.substr(sbjct.size()-l, l);
//         count = scoring::countPositive(qsub, ssub, BLOSUM62);    
//         //std::cout << qsub << "\t" << ssub << "\t" << count << "\n";
//         if ( count == l ) return IntPair( RIGHT, l );
    }
    return IntPair( -1, NOT_OVERLAP );
}

IntPair assem::findOverlap( PathId qpid,
                            PathId spid,
                            PathToAlnMap &path2aln_map,
                            int direction, 
                            Param &param )
{
    //if ( param.verbose ) std::cout << "\tOverlap checking\tquery:" << qpid << "\tsbjct:" << spid << "\n";

    //---------------------------------
    // Check long overlapping sequences
    // based on k-mer filtering
    //---------------------------------
    IntPair found = checkLongOverlap(qpid, spid, path2aln_map, param.kmer_size, direction, param );
    //std::cout << "Long overlap:?" << found.second << "\n";
    
    //-------------------------------
    // Check short substring matching
    //-------------------------------
    if ( found.second < 0 ) 
        found = checkShortOverlap(qpid, spid, path2aln_map, direction, param);
        //found = checkShortOverlap(qpid, spid, path2aln_map, 1, direction, param);

    return found;
}
                   

KmerId assem::getCenterKmer( PathId query_path, 
                      PathId sbjct_path, 
                      int overlap,
                      int direction, 
                      PathToAlnMap &path2aln_map,
					  Param &param )
{
//     std::string seq1 = biostr::stripGap( biostr::getSequenceString( path2aln_map[query_path]->getKmers(),path2aln_map[query_path]->getKmerCount(), param.kmer_size) );
//     std::string seq2 = biostr::stripGap( biostr::getSequenceString( path2aln_map[sbjct_path]->getKmers(),path2aln_map[sbjct_path]->getKmerCount(), param.kmer_size) );
    std::string query = biostr::getSequenceString( path2aln_map[query_path]->getKmers(),path2aln_map[query_path]->getKmerCount(), param.kmer_size);
    std::string sbjct = biostr::getSequenceString( path2aln_map[sbjct_path]->getKmers(),path2aln_map[sbjct_path]->getKmerCount(), param.kmer_size);

    std::string center, qsub, ssub;
    if ( overlap >= param.kmer_size ) {
        assert( (int)query.size() >= param.kmer_size && (int)sbjct.size() >= param.kmer_size );
        //direction == NEGATIVE ? center = sbjct.substr(0, param.kmer_size) : center = query.substr(0, param.kmer_size);
        int i = 0;
        while (true) {
            if ( direction == NEGATIVE ) {
                if ( i >= (int)sbjct.size()-param.kmer_size ) break;
                if ( i > overlap-param.kmer_size ) break;
                ssub = sbjct.substr(i, param.kmer_size); 
                qsub = query.substr(query.size()-overlap+i, param.kmer_size);
                if ( ssub == qsub ) return alpha::AminoAcidToInteger<KmerId>(ssub);
            } else {
                if ( i >= (int)query.size()-param.kmer_size ) break;
                if ( i > overlap-param.kmer_size ) break;
                qsub = query.substr(i, param.kmer_size); 
                ssub = sbjct.substr(sbjct.size()-overlap+i, param.kmer_size);
                if ( ssub == qsub ) return alpha::AminoAcidToInteger<KmerId>(ssub);
            }
            i++;
        }
        // in case of not exact match found, just set first kmer
        direction == NEGATIVE ? center = sbjct.substr(i, param.kmer_size) : center = query.substr(i, param.kmer_size);
    } else {
        assert( (int)query.size() >= overlap && (int)sbjct.size() >= overlap );
        direction == NEGATIVE ? center = sbjct.substr(0, overlap) : center = query.substr(0, overlap);
        int i = overlap;
        while (true) {
            assert( (int)query.size() > i && (int)sbjct.size() > i );
            if ( direction == NEGATIVE ) {
                center = query.substr( query.size()-1-i, 1 ) + center;
                if ( (int)center.size() == param.kmer_size ) break;
                center += sbjct.substr( i, 1 );
                if ( (int)center.size() == param.kmer_size ) break;
            } else {
                center = sbjct.substr( sbjct.size()-1-i, 1 ) + center;
                if ( (int)center.size() == param.kmer_size ) break;
                center += query.substr( i, 1 );
                if ( (int)center.size() == param.kmer_size ) break;
            }
            i++;
        }
    }
    return alpha::AminoAcidToInteger<KmerId>(center);
}

bool assem::__checkLatchableRead( int &start,
                                 std::vector<KmerId> &kmers, 
                                 std::vector<int> &pos_vec, 
                                 PathId &path, 
                                 PathToAlnMap &path2aln_map, 
                                 Param &param, 
                                 int direction)
{
    //bool latch = false;
    
    IntPair srange;
    direction == LEFT ? srange = getMatchRangeLeft(pos_vec, param) :
        srange = getMatchRangeRight(pos_vec, param) ;

    if ( srange.first != -1 || srange.second != -1 ) return false;
    
    IntPair qrange; 
    qrange.first  = pos_vec[srange.first];
    qrange.second = pos_vec[srange.second];

    int type = getMatchType( qrange, srange, &kmers[0], path2aln_map[path]->getKmers(), kmers.size(), path2aln_map[path]->getKmerCount(), param );
    
    if ( direction == LEFT ) {
        if ( type != LATCH_LEFT_EASY && type != LATCH_LEFT_DIFF ) 
            return false;
        else {
            start = srange.first;
            return true;
        }
    } else {
        if ( type != LATCH_RIGHT_EASY && type != LATCH_RIGHT_DIFF ) 
            return false;
        return true;
    }

    return false;
}

void assem::getLatchReads( //ReadIdArray &latch_reads,
                           //std::vector<int> &read_inits,
                           BitString *bstrs,
                           PathId *used_reads,
                           PathId sbjct_path, 
                           PairedPath &pair_path,
                           //PathId sbjct_path,
                           KmerToPathMap &pathid_map,
                           PathToAlnMap &path2aln_map, 
                           PathIdSet &merged_paths,
                           InvertedIndex &iindex,
                           size_t latch_offset, Param &param)
{
    KmerId kid = getCenterKmer( pair_path.pid, sbjct_path, pair_path.overlap, pair_path.direction, path2aln_map, param );
    if ( param.verbose ) std::cout << "Center Kmer:" << alpha::IntegerToAminoAcid(kid, param.kmer_size) << "\n";
    if ( ! iindex.has(kid) ) { 
        if ( param.verbose ) std::cout << "Index not exist\n"; 
        return;
    }
    
    PathId query_path = pair_path.pid;
//     if ( pair_path.direction == POSITIVE ) {
//         if ( param.verbose ) std::cout << "Swapped\n";
//         util::swap<PathId>(query_path, sbjct_path);
//     }

//     std::string qstr = biostr::stripGap( biostr::getSequenceString( path2aln_map[query_path]->getKmers(), path2aln_map[query_path]->getKmerCount(), param.kmer_size ) );
//     std::string sstr = biostr::stripGap( biostr::getSequenceString( path2aln_map[sbjct_path]->getKmers(), path2aln_map[sbjct_path]->getKmerCount(), param.kmer_size ) );
    std::string qstr = biostr::getSequenceString( path2aln_map[query_path]->getKmers(), path2aln_map[query_path]->getKmerCount(), param.kmer_size );
    std::string sstr = biostr::getSequenceString( path2aln_map[sbjct_path]->getKmers(), path2aln_map[sbjct_path]->getKmerCount(), param.kmer_size );

    ReadIdArray reads = ReadIdArray( iindex.getValue(kid)->rid, iindex.getValue(kid)->rid + iindex.getValue(kid)->size);
    
    ReadIdList latch;
    std::list<int> inits;
    for ( size_t i = 0; i < reads.size(); i++ ) {
        //--------------------------------------
        // ignore if it is already in other path.
        //--------------------------------------
        if ( used_reads[reads[i]] != NOT_PATH ) continue;
             
        if (param.verbose) std::cout << "Read:" << reads[i] << "\t" << bstrs[reads[i]].toString() << "\n";

        std::string rstr = bstrs[reads[i]].toString();
        KmerArray kmers = biostr::getKmers( rstr, param.kmer_size );
        PosPathPairList pos_paths = searchSimilarPaths(PathIdSet(), &kmers[0], kmers.size(), path2aln_map, pathid_map, merged_paths, 1, true, param);
        if ( param.verbose ) std::cout << "\tLatch reads:" << pos_paths.size() << "\n";
        bool qfound = false;
        bool mfound = false;
        int start;
        for ( PosPathPairList::iterator pt = pos_paths.begin(); pt != pos_paths.end(); ++pt ) {
            if ( param.verbose ) std::cout << "\tPath:" << pt->second << "\n";

            if ( pt->second == query_path ) {
                qfound = __checkLatchableRead(start, kmers, pt->first, query_path, path2aln_map, param, RIGHT);
            } 
            else if ( pt->second == pair_path.pid ) {
                mfound = __checkLatchableRead(start, kmers, pt->first, pair_path.pid, path2aln_map, param, LEFT);
            }
            else continue;


//             IntPair rrange = getMatchRange(pt->first, param);
//             if ( rrange.first == -1 || rrange.second == -1 ) continue;
//             IntPair qrange; 
//             qrange.first = (*pt).first[rrange.first];
//             qrange.second = (*pt).first[rrange.second];

//             //if ( LATCH != getLinkType( qrange, rrange, &kmers[0], path2aln_map[pt->second]->getKmers(), kmers.size(), path2aln_map[pt->second]->getKmerCount(), latch_offset, param.path_spur, param ) )
//             int type = getPathPairType( pt->first, qrange, rrange, &kmers[0], path2aln_map[pt->second]->getKmers(), kmers.size(), path2aln_map[pt->second]->getKmerCount(), latch_offset, param.path_spur, param );
//             if ( param.verbose ) std::cout << "Type:" << type << "\n";
//             if ( type != LATCH ) {
//                 std::cout << "Not a latch\n";
//                 continue;
//             }
//             if ( pt->second == query_path ) {
//                 qfound = true;
//                 start = rrange.first;
//                 std::cout << "Start:" << start << "\n";
//             }
//             if ( pt->second == sbjct_path ) {
//                 mfound = true;
//             }
        }

        if ( !qfound && !mfound ) continue;

        // if k-mer based filtering is failed check substr match
        if ( !qfound ) {
            for ( size_t k = param.kmer_size-1; k > 0; k-- ) {
                assert( qstr.size() > k && rstr.size() >= k );
                if ( qstr.substr( qstr.size()-k, k ) == rstr.substr(0,k) ) {
                    if (param.verbose) std::cout << "Qfound:Left:" << k << "\n";
                    start = qstr.size()-k;
                    qfound = true; break;
                }
            }
        }
        if ( !mfound ) {
            for ( size_t k = param.kmer_size-1; k > 0; k-- ) {
                assert( sstr.size() >= k && rstr.size() > k );
                if ( sstr.substr( 0, k ) == rstr.substr(rstr.size()-k,k) ) {
                    if (param.verbose) std::cout << "Mfound:Left:" << k << "\n";
                    mfound = true; break;
                }
            }
        }

        if ( qfound && mfound ) {
            if (param.verbose) {
                std::cout << "Good:" << bstrs[reads[i]].toString() << "\n";
                std::cout << "Start:" << start << "\n";
            }
            latch.push_back(reads[i]);
            inits.push_back(start);
        }
    }
    if (param.verbose) std::cout << "# Good reads:" << latch.size() << "\n";
    pair_path.latch_reads = ReadIdArray( latch.begin(), latch.end() );
    pair_path.inits = std::vector<int>( inits.begin(), inits.end() );
}

ReadIdArray assem::getBridgeReads( PathToAlnMap &path2aln_map, 
                                   PathId query_path, 
                                   PathId sbjct_path, 
                                   InvertedIndex &iindex, 
                                   PathId *used_reads,
                                   int direction, 
                                   int length,
                                   Param &param)
{
    KmerId *qkmers = path2aln_map[query_path]->getKmers();
    size_t  qnkmer = path2aln_map[query_path]->getKmerCount();
    KmerId *skmers = path2aln_map[sbjct_path]->getKmers();
    size_t  snkmer = path2aln_map[sbjct_path]->getKmerCount();

    if ( (int)qnkmer+param.kmer_size-1 < length || (int)snkmer+param.kmer_size-1 < length ) return ReadIdArray();
    std::string ssub, qsub;
    std::tr1::unordered_map<ReadId, bool> sreads, qreads;
    size_t sbeg, send, qbeg, qend;
    if ( direction == POSITIVE ) {
        sbeg = (snkmer+param.kmer_size-1) - length; 
        qbeg = 0;
        send = snkmer;
        qend = length-param.kmer_size+1;
    } else {
        sbeg = 0; 
        qbeg = (qnkmer+param.kmer_size-1) - length; 
        send = length-param.kmer_size+1;
        qend = qnkmer;
    }
    
    for ( size_t i = sbeg; i < send; i++ ) {
        KmerId kid = skmers[i];
        if ( iindex.has(kid) ) {
            ReadIdArray reads = ReadIdArray( iindex.getValue(kid)->rid, iindex.getValue(kid)->rid + iindex.getValue(kid)->size);            
            for ( ReadIdArray::iterator it = reads.begin(); it != reads.end(); ++it )
                sreads.insert( std::pair<ReadId, bool>(*it, true) );
        }
    }
    
    for ( size_t i = qbeg; i < qend; i++ ) {
        KmerId kid = qkmers[i];
        if ( iindex.has(kid) ) {
            ReadIdArray reads = ReadIdArray( iindex.getValue(kid)->rid, iindex.getValue(kid)->rid + iindex.getValue(kid)->size);            
            for ( ReadIdArray::iterator it = reads.begin(); it != reads.end(); ++it )
                qreads.insert( std::pair<ReadId, bool>(*it, true) );
        }
    }

    ReadIdArray comm_reads;
    for ( std::tr1::unordered_map<ReadId, bool>::iterator it = sreads.begin(); it != sreads.end(); ++it ) {
        if ( qreads.find(it->first) != qreads.end() ) {
            if ( used_reads[it->first] != NOT_PATH ) continue;
            comm_reads.push_back(it->first);
        }
    }

    if ( param.verbose ) std::cout << "direction:" << direction << "\tcommon reads:" << comm_reads.size() << "\n";
    return comm_reads;
}

void assem::makePathReadMap( ReadIdArray &comm_reads,
                             PathId pid,
                             BitString *bstrs,
                             PathReadMap &latch,
                             PathToAlnMap &path2aln_map,
                             Param &param,
                             int direction )
{
    for ( ReadIdArray::iterator it = comm_reads.begin(); it != comm_reads.end(); ++it ) {
        std::string str = bstrs[*it].toString();    
        if ( param.verbose ) std::cout << "\nReadId:" << *it << "\tSequence:" << str << "\n";
        
        KmerArray kmers = biostr::getKmers( str, param.kmer_size );
        std::vector<int> pvec = makePosFlags(kmers.size(), path2aln_map[pid]->getKmerCount(), &kmers[0], path2aln_map[pid]->getKmers(), param.latch_length, param);

        IntPair qrange, rrange;
        int type;
        AlignSummary summary;
        if ( __goodReadToJoin( pvec, qrange, rrange, type, summary, &kmers[0], kmers.size(), path2aln_map[pid]->getKmers(), path2aln_map[pid]->getKmerCount(), param, ANCHOR_CENTER ) ) {
            if ( param.verbose ) std::cout << "direction:" << direction << "\ttype:" << type << "\n";
            if ( ( direction == RIGHT && type == LATCH_LEFT_EASY ) ||
                 ( direction == LEFT && type == LATCH_RIGHT_EASY ) ) {
                latch[pid].push_back( ReadPos(*it, qrange.first, rrange.first) );                               
            } else {
                if ( param.verbose ) std::cout << "Invalid latch read\n";
                continue;
            }
        }
    }
}

bool assem::latchPathsWithReads( PathId &sbjct_path, 
                                 PairedPath &pair_path,
                                 KmerToPathMap &pathid_map,
                                 PathToAlnMap &path2aln_map, 
                                 BitString *bstrs,
                                 char *strands,
                                 ReadId *pairs,
                                 PathId *used_reads,
                                 PathIdSet &merged_paths,
                                 PathIdSet &success_paths,
                                 InvertedIndex &iindex,
                                 int kmer_size,
                                 Param &param)
{
    PathId query_path = pair_path.pid;
    std::string qstr = biostr::getSequenceString( path2aln_map[query_path]->getKmers(), path2aln_map[query_path]->getKmerCount(), param.kmer_size );
    std::string sstr = biostr::getSequenceString( path2aln_map[sbjct_path]->getKmers(), path2aln_map[sbjct_path]->getKmerCount(), param.kmer_size );

    //size_t PATH_OVERLAP_REGION = 30;
    if ( qstr.size() < PATH_OVERLAP_REGION || sstr.size() < PATH_OVERLAP_REGION ) return false;
    
    ReadIdArray comm_reads = getBridgeReads( path2aln_map, query_path, sbjct_path, iindex, used_reads, pair_path.direction, PATH_OVERLAP_REGION, param );
    if ( (int)comm_reads.size() < param.latch_support ) {
        if ( param.verbose ) std::cout << "Weak pair-end read support\n";
        return false;
    }

    std::string lstr, rstr;
    PathId left, right;
    if ( pair_path.direction == POSITIVE ) {
        lstr = sstr; rstr = qstr;
        left = sbjct_path; right = query_path;
    } else {
        lstr = qstr; rstr = sstr;
        left = query_path; right = sbjct_path;
    }

    PathReadMap latch_left;
    PathReadMap latch_right;
    makePathReadMap( comm_reads, left, bstrs, latch_right, path2aln_map, param, LEFT );
    makePathReadMap( comm_reads, right, bstrs, latch_left, path2aln_map, param, RIGHT );
    
    ReadIdArray good_reads;
    std::vector<int> read_inits;
    int overlap;
        
    validLatchReads(used_reads, bstrs, lstr, rstr, left, right, latch_left, latch_right, comm_reads, good_reads, read_inits, overlap, param );

    if ( (int)good_reads.size() < param.latch_support ) {
        if ( param.verbose ) std::cout << "Weak supported latch:\n";
        return false;
    }
    
    //bool success = latch( bstrs, strands, pairs, path2aln_map, pathid_map, merged_paths, success_paths, left, right, used_reads, good_reads, read_inits, overlap, iindex, param);
    bool success = connectPathPair(bstrs, strands, pairs, path2aln_map, pathid_map, merged_paths, success_paths, sbjct_path, query_path, used_reads, good_reads, read_inits, overlap, pair_path.direction, iindex, param );
    
    if ( success ) {
        if ( param.verbose ) std::cout << "Latch Success\n";
        //if ( left != sbjct_path ) sbjct_path = left; // for recursive latch
        return true;
    }
    return false;
}

//====================================================================
// Determine the best path to be latched based on the following
// a. # paired reads
// b. # compatible paired reads
// c. Determine direction
// d. Find overlap
//====================================================================
std::list<PairedPath> assem::makePairedPaths( PathIdArray &pids,
                                              PathId &sbjct_pid,
                                              ReadId *pairs,
                                              PathId *used_reads,
                                              char *strands,
                                              PathToAlnMap &path2aln_map,
                                              Param &param )
{
    assert( path2aln_map.find(sbjct_pid) != path2aln_map.end() );
    std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();
    //std::string sbjct = biostr::getSequenceString( path2aln_map[sbjct_pid]->getKmers(), path2aln_map[sbjct_pid]->getKmerCount(), param.kmer_size );
    std::list<PairedPath> ppaths;
    for ( size_t i = 0; i < pids.size(); i++ ) {
        if ( path2aln_map.find(pids[i]) == path2aln_map.end() ) {
            //std::cout << "[WARN]:path " << pids[i] << " does not exist\n";
            continue;
        }

        //assert( path2aln_map.find(pids[i]) != path2aln_map.end() );
        std::string query = path2aln_map[pids[i]]->getConsensusString();
        //if ( query.size() > sbjct.size() ) continue;

        if ( param.verbose ) {
            std::cout << i << ". Paired PID:" << pids[i] << "\n";
            std::cout << "Query:" << query << "\n";
            //std::cout << "Query:" << biostr::getSequenceString( path2aln_map[pids[i]]->getKmers(), path2aln_map[pids[i]]->getKmerCount(), param.kmer_size )  << "\n";
        }
        ReadIdArray qreads, sreads;
        setPairedReads( qreads, sreads, pids[i], sbjct_pid, used_reads, pairs, path2aln_map );
        if ( param.verbose ) std::cout << "Pair reads:" << qreads.size() << "\n";
        int direction = getDirection( sreads, pairs, strands, param );
        if ( param.verbose ) {
            std::cout << "Direction:" <<  direction << "\n";
        }
        if ( direction == NONSENSE ) {
            if ( param.verbose ) std::cout << "** Same strands -> skip\n"; continue;
        }
        IntPair overlap = findOverlap( pids[i], sbjct_pid, path2aln_map, direction, param );
        if ( param.verbose ) std::cout << "Overlap:" << overlap.first << "\t" << overlap.second << "\n";
        //if ( overlap.second <= 0 ) continue;
        if ( //overlap.second > 0 &&
             ( ( direction == POSITIVE && (size_t)overlap.first == LEFT  ) ) || 
             ( ( direction == NEGATIVE && (size_t)overlap.first == RIGHT ) ) ) { 
            if ( param.verbose ) std::cout << "\tDirectionality conflicts\n"; 
            continue;
        }
        
        std::string lstr;
        if ( direction == POSITIVE ) lstr = path2aln_map[sbjct_pid]->getConsensus();
        else lstr = path2aln_map[pids[i]]->getConsensus();
        if ( lstr[lstr.length()-1] == '*' ) {
            if ( param.verbose ) std::cout << "\tStop codon. No latch here\n"; 
            continue;
        }

        ppaths.push_back(PairedPath( pids[i], direction, overlap.second, qreads, sreads, ReadIdArray(), std::vector<int>() ) ); 
    }
    return ppaths;
}

void assem::orderPathsByOverlap( std::multimap<int, PairedPath> &orders, std::list<PairedPath> &ppaths, Param &param  )
{
    if ( param.verbose ) std::cout << "# Merge candidate paths:" << ppaths.size() << "\n";
    PairedPath max_path;
    //std::multimap<int, PairedPath> orders;
    
    for ( std::list<PairedPath>::iterator it = ppaths.begin(); it != ppaths.end(); ++it ) {
        orders.insert( std::pair<int, PairedPath>( it->overlap, *it ) );
    }
}

void assem::orderPathsBySupport( std::multimap<int, PairedPath> &orders, std::list<PairedPath> &ppaths, Param &param  )
{
    if ( param.verbose ) std::cout << "# Merge candidate paths:" << ppaths.size() << "\n";
    PairedPath max_path;
    //std::multimap<int, PairedPath> orders;
    
    for ( std::list<PairedPath>::iterator it = ppaths.begin(); it != ppaths.end(); ++it ) {
        orders.insert( std::pair<int, PairedPath>( (*it).match_reads.size(), *it ) );
    }
}

PairedPath assem::determineBestPair( std::list<PairedPath> &ppaths, Param &param  )
{
    if ( param.verbose ) std::cout << "# Merge candidate paths:" << ppaths.size() << "\n";
    PairedPath max_path;
    for ( std::list<PairedPath>::iterator it = ppaths.begin(); it != ppaths.end(); ++it ) {
        if ( max_path.pid == NOT_PATH ) {
            max_path = *it;
            continue;
        }

        if ( (int)(*it).match_reads.size() < param.overlap_support ) continue;
        
        if ( (*it).overlap > max_path.overlap ) {
            max_path = *it;
        }
    }
    return max_path;
}

std::vector<int> assem::extractInits( ReadIdArray &rids, 
                               SpaPath *spath )
{
    std::list<int> starts;
    
    int    *inits = spath->getInits();
    size_t  nread = spath->getReadCount();
    ReadId *reads = spath->getReads();
    for ( size_t i = 0; i < rids.size(); i++ ) {
        for ( size_t j = 0; j < nread; j++ ) {
            if ( rids[i] == reads[j] ) {
                starts.push_back( inits[j] );
                break;
            }
        }
    }
    return std::vector<int>( starts.begin(), starts.end() );
}


// void assem::stitchPair( PathId &query_path, 
//                  PairedPath &pair_info,
//                  KmerToPathMap &pathid_map,
//                  PathToAlnMap &path2aln_map,
//                  BitString *bstrs,
//                  char *strands,
//                  ReadId *pairs,
//                  PathId *used_reads,
//                  PathIdSet &merged_paths,
//                  PathIdSet &success_paths,
//                  InvertedIndex &iindex,
// 				 Param &param )
// {
//     if ( param.verbose ) {
// //         std::cout << "Query:" << query_path << "\t" << biostr::stripGap( biostr::getSequenceString( path2aln_map[query_path]->getKmers(),path2aln_map[query_path]->getKmerCount(), param.kmer_size) ) << "\n";
// //         std::cout << "Sbjct:" << pair_info.pid << "\t" << biostr::stripGap( biostr::getSequenceString( path2aln_map[pair_info.pid]->getKmers(),path2aln_map[pair_info.pid]->getKmerCount(), param.kmer_size) ) << "\n";            
//         std::cout << "Query:" << query_path << "\t" << biostr::getSequenceString( path2aln_map[query_path]->getKmers(),path2aln_map[query_path]->getKmerCount(), param.kmer_size) << "\n";
//         std::cout << "Sbjct:" << pair_info.pid << "\t" << biostr::getSequenceString( path2aln_map[pair_info.pid]->getKmers(),path2aln_map[pair_info.pid]->getKmerCount(), param.kmer_size) << "\n";            
        
//         std::cout << "Overlap:" << pair_info.overlap << "\n";
//         std::cout << "Paired:" << pair_info.match_reads.size() << "\n";
//         std::cout << "Direction:" << pair_info.direction << "\n";
//         std::cout << "Latch:" << pair_info.latch_reads.size() << "\n";
//     }

//     PathId sbjct_path = pair_info.pid;
//     int direction = pair_info.direction;
//     ReadIdArray links = pair_info.latch_reads;
//     std::vector<int>  inits = pair_info.inits;
    
//     if ( direction == NEGATIVE ) util::swap<PathId>( query_path, sbjct_path );

    
//     latch( bstrs, strands, pairs, path2aln_map, pathid_map, merged_paths, success_paths, sbjct_path, query_path, used_reads, links, inits, pair_info.overlap, iindex, param );
//     query_path = sbjct_path;

//     //if (param.verbose) std::cout << "NewSeq:" << biostr::stripGap( biostr::getSequenceString( path2aln_map[query_path]->getKmers(),path2aln_map[query_path]->getKmerCount(), param.kmer_size) ) << "\n";
// if (param.verbose) std::cout << "NewSeq:" << biostr::getSequenceString( path2aln_map[query_path]->getKmers(),path2aln_map[query_path]->getKmerCount(), param.kmer_size) << "\n";
// }

bool assem::linkPairedPath( PathId &sbjct_path,
                            PairedPath &max_pair,
                            KmerToPathMap &pathid_map,
                            PathToAlnMap &path2aln_map, 
                            BitString *bstrs,
                            char *strands,
                            ReadId *pairs, 
                            PathId *used_reads,
                            int nreads, 
                            PathIdSet &merged_paths,
                            PathIdSet &success_paths,
                            InvertedIndex &iindex,
                            size_t latch_offset,
                            Param &param )
{
    if ( max_pair.overlap <= 0 ) return false;

    getLatchReads( bstrs, used_reads, sbjct_path, max_pair, pathid_map, path2aln_map, merged_paths, iindex, latch_offset, param );        
    //bool success = latch( bstrs, strands, pairs, path2aln_map, pathid_map, merged_paths, success_paths, query_path, max_pair.pid, used_reads, max_pair.latch_reads, max_pair.inits, max_pair.overlap, iindex, param );
    bool success = connectPathPair(bstrs, strands, pairs, path2aln_map, pathid_map, merged_paths, success_paths, sbjct_path, max_pair.pid, used_reads, max_pair.latch_reads, max_pair.inits, max_pair.overlap,max_pair.direction, iindex, param );
    return success;
}


// bool assem::latchPathsWithReads( PathId &sbjct_path,
//                                  PairedPath &max_pair,
//                                  KmerToPathMap &pathid_map,
//                                  PathToAlnMap &path2aln_map, 
//                                  BitString *bstrs,
//                                  char *strands,
//                                  ReadId *pairs, 
//                                  PathId *used_reads,
//                                  int nreads, 
//                                  PathIdSet &merged_paths,
//                                  PathIdSet &success_paths,
//                                  InvertedIndex &iindex,
//                                  size_t latch_offset,
//                                  Param &param )
// {
//     if ( max_pair.overlap <= 0 ) return false;

//     getBridgeReads( bstrs, used_reads, sbjct_path, max_pair, pathid_map, path2aln_map, merged_paths, iindex, latch_offset, param );        
//     //bool success = latch( bstrs, strands, pairs, path2aln_map, pathid_map, merged_paths, success_paths, query_path, max_pair.pid, used_reads, max_pair.latch_reads, max_pair.inits, max_pair.overlap, iindex, param );
//     bool success = connectPathPair(bstrs, strands, pairs, path2aln_map, pathid_map, merged_paths, success_paths, sbjct_path, max_pair.pid, used_reads, max_pair.latch_reads, max_pair.inits, max_pair.overlap,max_pair.direction, iindex, param );
//     return success;
// }

// bool assem::linkPairedPath( PathId &query_path,
//                             PairedPath &max_pair,
//                             KmerToPathMap &pathid_map,
//                             PathToAlnMap &path2aln_map, 
//                             BitString *bstrs,
//                             char *strands,
//                             ReadId *pairs, 
//                             PathId *used_reads,
//                             int nreads, 
//                             PathIdSet &merged_paths,
//                             PathIdSet &success_paths,
//                             InvertedIndex &iindex,
//                             size_t latch_offset,
//                             Param &param )
// {
//     if ( max_pair.overlap <= 0 ) return false;


//     bool swapped = false;
//     //PathId sbjct_pid = max_pair.pid;
//     //if ( path2aln_map[query_pid]->getKmerCount() < path2aln_map[max_pair.pid]->getKmercount() ) {

//     /* Make query goes left */
//     if ( max_pair.direction == POSITIVE ) {
//         //std::cout << "Query:" << query_path << "\tSbjct:" << max_pair.pid << "\n";
//         if ( param.verbose ) std::cout << "Swapped\n";
//         util::swap<PathId>(query_path, max_pair.pid);
//         //std::cout << "Query:" << query_path << "\tSbjct:" << max_pair.pid << "\n";
//         max_pair.direction = NEGATIVE;
//         swapped = true;
//     }
//     //ReadIdArray latch_rids;
//     //std::vector<int> inits;
//     //getLatchReads( latch_rids, inits, bstrs, used_reads, query_path, sbjct_path, pathid_map, path2aln_map, merged_paths, iindex, latch_offset, param );            
        
//     getLatchReads( bstrs, used_reads, query_path, max_pair, pathid_map, path2aln_map, merged_paths, iindex, latch_offset, param );        
// //         std::string qstr = biostr::stripGap( path2aln_map[query_path]->getConsensusString() );
// //         std::string sstr = biostr::stripGap( path2aln_map[max_pair.pid]->getConsensusString() );
// //         std::string qstr = path2aln_map[query_path]->getConsensusString();
// //         std::string sstr = path2aln_map[max_pair.pid]->getConsensusString();
        
// //         if ( param.verbose ) 
// //             std::cout << "Merging:query:" << query_path << "\tsbject:" << max_pair.pid 
// //                       << "\tqlen:" << qstr.length() << "\tslen:" << sstr.length() 
// //                       <<"\toverlap:" << max_pair.overlap 
// //                       << "\tqreads:" << path2aln_map[query_path]->getReadCount() << "\tsreads:" << path2aln_map[max_pair.pid]->getReadCount() 
// //                       << "\tpairs:" << max_pair.query_reads.size() << "\tlatch:" << max_pair.latch_reads.size() << "\n";
//     //stitchPair( query_path, max_pair, pathid_map, path2aln_map, bstrs, strands, pairs, used_reads, merged_paths, success_paths, iindex, param );
//     //return true;
//         //}
//     bool success = latch( bstrs, strands, pairs, path2aln_map, pathid_map, merged_paths, success_paths, query_path, max_pair.pid, used_reads, max_pair.latch_reads, max_pair.inits, max_pair.overlap, iindex, param );
//     if ( success ) {
//         // no need to update query
//         //query_path = sbjct_path;
//         return true;
//     } else {
//         if (swapped) util::swap<PathId>(query_path, max_pair.pid);
//     }
//     return false;
// }


// //====================================================================
// // Latch a pair of paths when sufficient paired end reads exists, even
// // if there is no supporting latch reads between two paths.
// // If two sequences are overlapped, latch the pair.
// // Otherwise, just print out the possible latching.
// //====================================================================
// void assem::extendShortOverlapPaths( PathToAlnMap &path2aln_map, 
//                                      BitString *bstrs,
//                                      char *strands,
//                                      ReadId *pairs, 
//                                      PathId *used_reads,
//                                      int nreads, 
//                                      InvertedIndex &iindex,
//                                      Param &param )
// {
//     if ( !param.extend_flag || !param.pair_flag ) return;

//     double t0 = mytime();

//     std::cout << "\nSTAGE 4:\nMerging short ovelapping paths ....\n";
    
//     PathIdSet merged_paths; // Keep track of merged paths
//     PathIdSet successPaths; // Keep track of new extended pathso

//     //-------------------------
//     // k-mer to path id mapping
//     //-------------------------
//     KmerToPathMap pathid_map;
//     for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) 
//         addPathIds( it->first, path2aln_map[it->first]->getKmers(), path2aln_map[it->first]->getKmerCount(), pathid_map );

//     /*********************************************************************
//       1. For each path, find paths with paired end reads supports
//       2. Determine the best path to be latched based on the following
//       	 a. # paired reads
//          b. # compatible paired reads
//          c. Determine direction
//          d. Find overlap
//       3. Latch the pair.
//       4. Fill gaps if exist.
//      *********************************************************************/
//     double pratio = 1.0;
//     int nproc = 0;
//     int npath = path2aln_map.size();
//     double lt0 = mytime();

//     PathLengthMap plen_map = getPathLengths( path2aln_map );
//     for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
//         PathId sbjct_pid = it->second;

//         //for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
//         nproc++;
//         if ( nproc/(double)npath*100 >= pratio ) {
//             fprintf( stdout, "%d/%d: %.2f sec\n", nproc, npath, mytime()-lt0);
//             pratio += 1;
//         }

//         if ( merged_paths.find(sbjct_pid) != merged_paths.end() ) continue;
        
//         //PathId query_path = it->first;

//         /* Recursive latching */
//         while(true) {
//             //if ( param.verbose) std::cout << "\nQuery:" << query_path << "\t" << biostr::stripGap( biostr::getSequenceString( path2aln_map[sbjct_pid]->getKmers(),path2aln_map[sbjct_pid]->getKmerCount(), param.kmer_size) ) << "\n";            
//             if ( param.verbose) std::cout << "\nSbjct:" << sbjct_pid << "\t" << biostr::getSequenceString( path2aln_map[sbjct_pid]->getKmers(),path2aln_map[sbjct_pid]->getKmerCount(), param.kmer_size) << "\n";            
            
//             PathIdArray pids = findPairedPaths( sbjct_pid, used_reads, pairs, path2aln_map[sbjct_pid]->getReads(), path2aln_map[sbjct_pid]->getReadCount(), param );

//             //if ( pids.size() == 0 ) break;

//             std::list<PairedPath> ppaths = makePairedPaths( pids, sbjct_pid, pairs, used_reads, strands, path2aln_map, param );
//             if ( ppaths.size() == 0 ) break;

            
//             PairedPath best_pair = determineBestPair( ppaths, param );
//             if (param.verbose) {
//                 //std::cout << "Max Sbjct:" << best_pair.pid << "\t" << biostr::stripGap( biostr::getSequenceString( path2aln_map[best_pair.pid]->getKmers(),path2aln_map[best_pair.pid]->getKmerCount(), param.kmer_size) ) << "\n";            
//                 std::cout << "Max Query:" << best_pair.pid << "\t" << biostr::getSequenceString( path2aln_map[best_pair.pid]->getKmers(),path2aln_map[best_pair.pid]->getKmerCount(), param.kmer_size) << "\n";            
//                 std::cout << "Max Overlap:" << best_pair.overlap << "\n";
//                 std::cout << "Max Support:" << best_pair.match_reads.size() << "\n";
//                 std::cout << "Max Links:" << best_pair.latch_reads.size() << "\n";
//                 std::cout << "Max Direction:" << best_pair.direction << "\n";
//             }

//             if ( best_pair.pid == NOT_PATH ) {
//                 if ( param.verbose ) std::cout << "Not a path\n";
//                 break;
//             }
//             if ( best_pair.overlap < param.pairend_overlap ) {
//                 if ( param.verbose ) std::cout << "Weak overlap\n"; 
//                 break;
//             }

//             if ( (int)best_pair.match_reads.size() < param.overlap_support ) {
//                 if ( param.verbose ) std::cout << "Weak support\n"; 
//                 break;
//             }
            
//             bool success = linkPairedPath( sbjct_pid, best_pair, pathid_map, path2aln_map, bstrs, strands, pairs, used_reads, nreads, merged_paths, successPaths, iindex, param.kmer_size, param );
//             if ( !success ) break;
//         }
        
// //         /***/
// //         dropPathIds( sbjct_pid, path2aln_map[sbjct_pid]->getKmers(), path2aln_map[sbjct_pid]->getKmerCount(), pathid_map );
//     }

//     dropMergedPaths(merged_paths, path2aln_map);


//     //--------------------------------------------------------------------
//     // Try to merge or latch similar pairs from read supported latch paths
//     //--------------------------------------------------------------------
//     double lt = mytime();
//     linkPaths( successPaths, path2aln_map, iindex, used_reads, bstrs, strands, pairs, nreads, param, MERGELATCH );
//     std::cout << "Paths merged/extended with newly extended paths:" << mytime()-lt << " sec\n";


//     if ( param.verbose ) 
//         printPaths(path2aln_map, param);
//     std::cout << "Merged:" << merged_paths.size() << "\n";
//     std::cout << "Path count:" << path2aln_map.size() << "\n";
//     std::cout << "Merging short overlapping paths:" << mytime()-t0 << " sec\n";
// }



void assem::scaffoldPaths(PathToAlnMap &path2aln_map, 
                          BitString *bstrs,
                          char *strands,
                          ReadId *pairs, 
                          PathId *used_reads,
                          Param &param )
{
    if ( !param.scaffold_flag || !param.pair_flag ) return;
    double t0 = mytime();

    setbuf(stdout, NULL); // no buffering
    std::cout << "\nStage 5:\nScaffolding paths ....\n";
    if ( param.verbose ) std::cout << "#Total paths:" << path2aln_map.size() << "\n";

    ScaffoldList scaffolds;

    PathIdSet merged_paths; // Keep track of merged paths
    Map<PathIdPair, bool>::Type scanned_pairs;

    
    /* auxilary variables for progress log */
    double pratio = 1.0;
    int    nproc  = 0;
    double lt0    = mytime();
    size_t npath  = path2aln_map.size();

    /* Make a path length map */
    PathLengthMap plen_map = getPathLengths( path2aln_map );
    for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        /* Print progress status at every percent of path proceeded */
        nproc++;
        if ( nproc/(double)npath*100 >= pratio ) {
            fprintf( stdout, "%d/%d: %.2f sec\n", nproc, int(npath), mytime()-lt0);
            pratio += 1;
        }

        PathId sbjct_pid = it->second;
        std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();
        if ( merged_paths.find( sbjct_pid ) != merged_paths.end() ) continue;

        if ( param.verbose ) {
            std::cout << "\nCenter:" << sbjct_pid << "\tLength:" << it->first << "\n";
            std::cout << sbjct << "\n";
        }

        ReadIdArray neg_rids, pos_rids;
        splitReadsByStrands(neg_rids, pos_rids, sbjct_pid, path2aln_map, strands);

        if ( param.verbose ) std::cout << "-reads:" << neg_rids.size() << "\t+reads:" << pos_rids.size() << "\n";
        neg_rids = dropSelfPairs( neg_rids, sbjct_pid, used_reads, pairs );
        pos_rids = dropSelfPairs( pos_rids, sbjct_pid, used_reads, pairs );
        if ( param.verbose ) std::cout << "-reads:" << neg_rids.size() << "\t+reads:" << pos_rids.size() << "\n";

//         if ( param.verbose ) std::cout << "-ReadId\tStrand\n";
//         for ( size_t i = 0; i < neg_rids.size(); i++ )
//             std::cout << neg_rids[i] << "\t" << strands[neg_rids[i]] << "\n";
//         if ( param.verbose ) std::cout << "+ReadId\tStrand\n";
//         for ( size_t i = 0; i < pos_rids.size(); i++ )
//             std::cout << pos_rids[i] << "\t" << strands[pos_rids[i]] << "\n";
                                                        
        Scaffold scaffold;
        scaffold.paths.push_back( sbjct_pid );
        scaffold.dists.push_back( 0 );

        __scaffold( scaffold, neg_rids, sbjct_pid, path2aln_map, bstrs, pairs, strands, used_reads, merged_paths, scanned_pairs, param, LEFT );
        __scaffold( scaffold, pos_rids, sbjct_pid, path2aln_map, bstrs, pairs, strands, used_reads, merged_paths, scanned_pairs, param, RIGHT );

        scaffolds.push_back(scaffold);
     }

    dropMergedPaths(merged_paths, path2aln_map);

    std::cout << "# Merged paths:" << merged_paths.size() << "\n";
    std::cout << "Path count:" << path2aln_map.size() << "\n";
    std::cout << "\nScaffolding paths:" << mytime()-t0 << " sec\n";
}

void assem::splitReadsByStrands( ReadIdArray &neg_rids,
                                 ReadIdArray &pos_rids,
                                 PathId pid,
                                 PathToAlnMap &path2aln_map,
                                 char *strands )
{
    ReadId *reads = path2aln_map[pid]->getReads();
    size_t nreads = path2aln_map[pid]->getReadCount();
    for ( size_t i = 0; i < nreads; i++ ) {
        if ( strands[reads[i]] == '+' )
            pos_rids.push_back( reads[i] );
        else if ( strands[reads[i]] == '-' )
            neg_rids.push_back( reads[i] );
    }
}

void assem::__scaffold( Scaffold &scaffold,
                        ReadIdArray &rids,
                        PathId sbjct_pid,
                        PathToAlnMap &path2aln_map,
                        BitString *bstrs,
                        ReadId *pairs,
                        char *strands, 
                        PathId *used_reads, 
                        PathIdSet &merged_paths,
                        Map<PathIdPair, bool>::Type &scanned_pairs,
                        Param &param,
                        int direction )
{
    if ( rids.size() == 0 ) return;

    if ( param.verbose ) std::cout << "Direction:" << direction << "\n";

    std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();
    if ( sbjct.size() == 0 ) return;
    if ( direction == RIGHT && sbjct[sbjct.size()-1] == '*' ) {
        if ( param.verbose ) std::cout << "Stop codon. No scaffold right\n";
        return;
    }

    //rids = ( direction==LEFT ) ? dropConflictReads(rids, strands, '-' ) : dropConflictReads(rids, strands, '+' );
    //  if ( rids.size() == 0 ) return;

//     /* Recursive latching */
//     while(true) {
    
    PathIdArray pids = findPairedPaths( sbjct_pid, used_reads, pairs, &rids[0], rids.size(), param );
    if ( pids.size() == 0 ) return;
    

    PathIdList candidates;
    for ( PathIdArray::iterator pt = pids.begin(); pt != pids.end(); ++pt ) {
        PathIdPair ppair = sbjct_pid > *pt ? PathIdPair(*pt, sbjct_pid) : PathIdPair(sbjct_pid, *pt);
        if ( scanned_pairs.find( ppair ) == scanned_pairs.end() && 
            merged_paths.find( *pt ) == merged_paths.end() )
                candidates.push_back( *pt );
    }
    
    //pids = PathIdArray( candidates.begin(), candidates.end() );
    if ( candidates.size() == 0 ) return;
        
    std::multimap<size_t, PathId> orders;    
    std::tr1::unordered_map<PathId, ReadIdArray> supports;
    std::tr1::unordered_map<PathId, int> overlaps;
    std::tr1::unordered_map<PathId, int> ngaps;
    for ( PathIdList::iterator it = candidates.begin(); it != candidates.end(); ++it ) {
        std::string query = path2aln_map[*it]->getConsensusString();
        if ( param.verbose ) std::cout << "Candidate:" << *it << "\n";
        if ( param.verbose ) std::cout << query << "\n";
        if ( query.size() == 0 ) continue;

        if ( direction == LEFT && query[query.size()-1] == '*' ) {
            if ( param.verbose ) std::cout << "Stop codon. Invalid candidate\n";
            continue;
        }

        ReadIdArray qpreads = getPairedReads( rids, sbjct_pid, *it, used_reads, pairs, path2aln_map );
        if ( param.verbose ) std::cout << "Pair reads:" << qpreads.size() << "\t";

        qpreads = ( direction==LEFT ) ? dropConflictReads(qpreads, strands, '+' ) : dropConflictReads(qpreads, strands, '-' );

        if ( param.verbose ) std::cout << "Pair reads:" << qpreads.size() << "\n";

        if ( (int)qpreads.size() < param.latch_support ) {
            PathIdPair ppair = sbjct_pid > *it ? PathIdPair(*it, sbjct_pid) : PathIdPair(sbjct_pid, *it);
            scanned_pairs.insert( std::pair<PathIdPair, bool>( ppair, true) );
            continue;
        }
        IntPair overlap = findOverlap( *it, sbjct_pid, path2aln_map, direction, param );
        if ( param.verbose ) std::cout << "Overlap:" << overlap.first << "\t" << overlap.second << "\n";

        int ngap = estimateGap( sbjct_pid, *it, qpreads, path2aln_map, bstrs, strands, pairs, direction, param );
        if ( param.verbose ) std::cout << "Gap estimate:" << ngap << "\n";

        orders.insert( std::pair<size_t, PathId>( qpreads.size(), *it ) );
        supports.insert( std::pair<PathId, ReadIdArray>( *it, qpreads ) );
        overlaps.insert( std::pair<PathId, int>( *it, overlap.second ) );
        ngaps.insert( std::pair<PathId, int>( *it, ngap ) );
    }

    if ( orders.size() == 0 ) return;

    if ( param.verbose ) std::cout << "# good candidates:" << orders.size() << "\n";
    for ( std::multimap<size_t, PathId>::reverse_iterator ot = orders.rbegin(); ot != orders.rend(); ++ot ) {
        if ( param.verbose ) std::cout << "Query:" << ot->second << "\tSupport:" << ot->first << "\tGap:" << ngaps[ot->second] << "\n";
    }

    bool success = false;
    for ( std::multimap<size_t, PathId>::reverse_iterator ot = orders.rbegin(); ot != orders.rend(); ++ot ) {

        PathId query = ot->second;
        size_t support = ot->first;

        PathIdPair ppair = sbjct_pid > query ? PathIdPair(query, sbjct_pid) : PathIdPair(sbjct_pid, query);
        scanned_pairs.insert( std::pair<PathIdPair, bool>( ppair, true) );
        
        if (param.verbose) {
            std::cout << "Query:" << query << "\t" << path2aln_map[query]->getConsensusString() << "\n";
            std::cout << "Overlap:" <<  overlaps[query] << "\n";
            std::cout << "Support:" << support << "\n";
        }
            
        if ( overlaps[query] < param.pairend_overlap ) {
            if ( param.verbose ) std::cout << "Weak overlap\n"; 
            //continue;
        }

        //sbjct_pid = query;
        ReadIdArray neg_rids, pos_rids;
        splitReadsByStrands(neg_rids, pos_rids, query, path2aln_map, strands);

        neg_rids = dropSelfPairs( neg_rids, query, used_reads, pairs );
        pos_rids = dropSelfPairs( pos_rids, query, used_reads, pairs );

        //rids = (direction==LEFT) ? neg_rids : pos_rids;


        //__scaffold( rids, sbjct_pid, path2aln_map, bstrs, pairs, strands, used_reads, merged_paths, scanned_pairs, param, direction );
        if ( direction == LEFT )
            __scaffold( scaffold, neg_rids, query, path2aln_map, bstrs, pairs, strands, used_reads, merged_paths, scanned_pairs, param, direction );
        else 
            __scaffold( scaffold, pos_rids, query, path2aln_map, bstrs, pairs, strands, used_reads, merged_paths, scanned_pairs, param, direction );
            
        //             success = linkPairedPath( sbjct_pid, best_pair, pathid_map, path2aln_map, bstrs, strands, pairs, used_reads, nreads, merged_paths, successPaths, iindex, param.kmer_size, param );
        //             if ( success ) {
        //                 scanned_pairs.insert( std::pair<PathIdPair, bool>( ppair, true) );
        //                 if ( param.verbose ) std::cout << "Search recursive latch for path:" << sbjct_pid << "\n";
        //                 break;
        //             }
    }
    if ( !success ) {
        if ( param.verbose ) std::cout << "Stop scaffolding\n";
        return;
    }
}

ReadIdArray assem::dropConflictReads( ReadIdArray &qreads, 
                                      char *strands, 
                                      char strand )
{
    ReadIdList good;
    for ( size_t i = 0; i < qreads.size(); i++ ) {
        if ( strands[qreads[i]] == strand ) 
            good.push_back( qreads[i] );
    }
    return ReadIdArray( good.begin(), good.end() );
}

int assem::estimateGap( PathId sbjct_pid, 
                        PathId query_pid, 
                        ReadIdArray &qpreads, 
                        PathToAlnMap &path2aln_map,
                        BitString *bstrs,
                        char *strands,
                        ReadId *pairs,
                        int direction,
                        Param &param)
{
    // min/max library insert size
    // one std. dev: 90%
    //double min_insert = (param.insert_size - param.insert_sd)/3;
    //double max_insert = (param.insert_size + param.insert_sd)/3;
    
    std::string query = path2aln_map[query_pid]->getConsensusString();
    std::string sbjct = path2aln_map[sbjct_pid]->getConsensusString();

    if ( param.verbose ) std::cout << "Lengths\tquery:" << query.size() << "\tsbjct:" << sbjct.size() << "\n";
    std::tr1::unordered_map<ReadId, int> query_starts = getReadInitPosMap( path2aln_map[query_pid] );
    std::tr1::unordered_map<ReadId, int> sbjct_starts = getReadInitPosMap( path2aln_map[sbjct_pid] );

    std::vector<double> min_gaps, max_gaps, avg_gaps;
    for ( size_t i = 0; i < qpreads.size(); i++ ) {
        ReadId query_rid = qpreads[i];
        ReadId sbjct_rid = pairs[query_rid];

        int qinit = query_starts[query_rid];
        int sinit = sbjct_starts[sbjct_rid];

        std::string qread = bstrs[query_rid].toString();
        std::string sread = bstrs[sbjct_rid].toString();

        //size_t qlen = qread.size();
        //size_t slen = sread.size();

        int lend, rend;
        int min_gap, max_gap;
        double avg_gap;
        if ( direction == LEFT ) {
            //lend = (query.size()-1) - qinit;
            rend = sinit;
            lend = query.size() - (qinit+qread.size());
        }
        else {
            //lend = (sbjct.size()-1) - sinit;
            rend = qinit;            
            lend = sbjct.size() - (sinit+sread.size());
        }
        if ( param.verbose ) std::cout << "query:" << query_rid << "\tsbjct:" << sbjct_rid 
                                       << "\tstarts:" << qinit << "," << sinit 
                                       << "\tstrands:" << strands[query_rid] << "," << strands[sbjct_rid] 
                                       << "\trange:" << lend << "," << rend << "\n";
//         min_gap = (min_insert+qlen+slen) - (lend+rend);
//         max_gap = (max_insert+qlen+slen) - (lend+rend);
        min_gap = param.insert_size/3 - param.insert_sd/3 - lend - rend;
        max_gap = param.insert_size/3 + param.insert_sd/3 - lend - rend;
        avg_gap = (min_gap+max_gap)/2;
        min_gaps.push_back(min_gap);
        max_gaps.push_back(max_gap);
        avg_gaps.push_back(avg_gap);
    }
    
    double min_avg = math::mean(&min_gaps[0], min_gaps.size());
    double max_avg = math::mean(&max_gaps[0], max_gaps.size());
    double avg_avg = math::mean(&avg_gaps[0], avg_gaps.size());
    if ( param.verbose ) std::cout << "Min avg gap:" << min_avg << "\tMax avg gap:" << max_avg << "\tAvg of avg:" << avg_avg << "\n";
    int avg2 = (int)(min_avg + max_avg)/2;
    if ( param.verbose ) std::cout << "Mean:" << avg2 << "\n";

    return avg2;
}

std::tr1::unordered_map<ReadId, int> assem::getReadInitPosMap( SpaPath *path )
{
    std::tr1::unordered_map<ReadId, int> spos_map;
    ReadId *reads = path->getReads();
    size_t nreads = path->getReadCount();
    int    *inits = path->getInits();

    for ( size_t i = 0; i < nreads; i++ ) {
        //std::cout << reads[i] << "\t" << inits[i] << "\n";
        spos_map.insert( std::pair<ReadId, int>( reads[i], inits[i] ) );
    }
    //std::cout << "\n";
    
    return spos_map;
}

ReadIdArray assem::dropSelfPairs( ReadIdArray &orids,
                           PathId pid,
                           PathId *used_reads,
                           ReadId *pairs )
{
    ReadIdList nrids;
    for ( size_t i = 0; i < orids.size(); i++ ) {
        ReadId pair_read = pairs[orids[i]];
        if ( pair_read == NOT_PAIR ) continue; 
        if ( used_reads[pair_read] == pid ) continue;
        nrids.push_back( orids[i] );
    }
    return ReadIdArray( nrids.begin(), nrids.end() );
}

void assem::writeConsensus( std::fstream &out, SpaPath *spath , int count )
{
//     // debug 2
//     std::string orf = spath->getConsensus();
//     //std::string orf = spath->getConsensusString(); //??? this causing problem??? 
//     assert( orf.size() > 0 );

//     orf = biostr::stripGap(orf);

    // debug 1 //
    char *orf = spath->getConsensus();
    assert( strlen(orf) > 0 ) ;
    
    out << ">" << count << "\n";
    out << orf << "\n";
    out.flush();
}


void assem::writePlacement( std::fstream &out, SpaPath *spath, int count )
{
    size_t nreads  = spath->getReadCount();
    ReadId *reads  = spath->getReads();
    int    *inits  = spath->getInits();
    Mismatch *inss = spath->getInsertions();
    Mismatch *dels = spath->getDeletions();
    unsigned cins  = spath->countInsertions();
    unsigned cdel  = spath->countDeletions();

    std::tr1::unordered_map<ReadId, std::set<int> > iPos, dPos;
    for ( unsigned i = 0; i < cins; i++ ) {
        ReadId rid = inss[i].read;
        int   spos = inss[i].spos;
        iPos[rid].insert(spos);
    }
    for ( unsigned i = 0; i < cdel; i++ ) {
        ReadId rid = dels[i].read;
        int   spos = dels[i].spos;
        dPos[rid].insert(spos);
    }
    out << "ID:" << count << "\n";
    out << "Raw:" << spath->getConsensus() << "\n";
    out << "Count:" << nreads << "\n";
    for ( size_t i = 0; i < nreads; i++ ) {
        out << i << "\t" << reads[i] << "\t" << inits[i] << "\tI:";
        for ( std::set<int>::iterator st = iPos[reads[i]].begin(); st != iPos[reads[i]].end(); ++st )
            out << *st << ",";
        out << "\tD:";
        for ( std::set<int>::iterator st = dPos[reads[i]].begin(); st != dPos[reads[i]].end(); ++st )
            out << *st << ",";
        out << "\n";
    }
    out << "//\n";
    out.flush();
}


// void writeStatistic( std::fstream &out, MSA &msa, int count )
// {
//     out << "ID:" << count << "\n";
//     std::string raw = msa.getConsensus();
//     std::string orf = biostr::stripGap(raw);
//     out << "Length:\t" << "Raw:" << raw.length() << "\tClean:" << orf.length() << "\n";
    
//     Profile profile = msa.getProfile();
//     ScoreSummary s = assem::computeEntropy( profile );
//     out << "Entropy:\t"<< "Sum:" << s.sum << "\tMin:" << s.min << "\tMax:" << s.max << "\tAvg:" << s.avg << "\tMed:" << s.med << "\n";
//     s = assem::computeDepth(profile);
//     out << "Coverage:\t"<< "Sum:" << s.sum << "\tMin:" << s.min << "\tMax:" << s.max << "\tAvg:" << s.avg << "\tMed:" << s.med << "\n";
//     out << "//\n";
// }

void assem::writeStatistic( std::fstream &out, SpaPath *spath, int count )
{
    out << "ID:" << count << "\n";

//     // debug 2
//     std::string raw = spath->getConsensus();
//     std::string orf = biostr::stripGap(raw);
//     out << "Length:\t" << "Raw:" << raw.length() << "\tClean:" << orf.length() << "\n";

    // debug 1
    char *orf = spath->getConsensus();
    out << "Length:" << strlen(orf) << "\n";

    ProfileVector *vprof = spath->getProfileVector();
    //Profile profile;
    //vprof->convert(profile);
    Profile profile = vprof->convert();
    ScoreSummary s = assem::computeEntropy( profile );
    out << "Entropy:\t"<< "Sum:" << s.sum << "\tMin:" << s.min << "\tMax:" << s.max << "\tAvg:" << s.avg << "\tMed:" << s.med << "\n";
    s = assem::computeDepth(profile);
    out << "Coverage:\t"<< "Sum:" << s.sum << "\tMin:" << s.min << "\tMax:" << s.max << "\tAvg:" << s.avg << "\tMed:" << s.med << "\n";
    out << "//\n";
    out.flush();

    profile.clean();
}

// void writeAlignment( std::fstream &out, MSA &msa, SpaPath *spath, int count )
// {
//     out << "ID:" << count << "\n";
//     msa.printAlignment( out, *spath, 100 );
//     out << "//\n";
// }
void assem::writeAlignment( std::fstream &out, SpaPath *spath, int count, BitString *bstrs )
{
    out << "ID:" << count << "\n";
    spath->printAlignment( out, 100, bstrs );
    out << "//\n";
    out.flush();
}

void assem::writeProfile( std::fstream &out, SpaPath *spath, int count )
{
    out << "ID:" << count << "\n";
    ProfileVector *vprof = spath->getProfileVector();
    //Profile profile;
    //vprof->convert(profile);
    Profile profile = vprof->convert();
    profile.print(out);
    out << "//\n";
    profile.clean();
    out.flush();
}
