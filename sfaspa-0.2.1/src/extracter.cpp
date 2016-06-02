#include "extracter.h"

Extracter::Extracter()
{
    preads = NULL;
    status = false;

    read_starts = NULL;
}

Extracter::Extracter(Loader *l)
{
    init(l);
}

Extracter::~Extracter()
{
    clear();
}

void Extracter::clear()
{
    release();

    path_map.clear();
    
    if ( preads != NULL ) {
        delete[] preads;
        preads = NULL;
    }

    if ( read_starts != NULL ) {
        delete[] read_starts;
        read_starts = NULL;
    }

    // if ( read_ends != NULL ) {
    //     delete[] read_ends;
    //     read_ends = NULL;
    // }

    status = false;

    if ( Param::verbose ) std::cout << "Extracter cleaned\n";
}

//====================================================================
// Initialization
//====================================================================
void Extracter::init(Loader *l)
{
    //-------------------------------
    // Set pointers to loader objects
    //------------------------------- 
    loader = l;
    seqs   = loader->getReads();
    nreads = loader->getCount();
    gsa    = loader->getGSA();

    read_starts = new ReadStartCount[Param::nparts];
    //read_ends = new ReadEndCount[Param::nparts];

    //----------------------------
    // Set sequence pointer to GSA
    //----------------------------
    for ( int i = 0; i < Param::nparts; i++ ) {
        //---------------------------------------------
        // Find the first sequence of each suffix array
        //---------------------------------------------
        int s = i*int(nreads/Param::nparts);

        //---------------------------------------
        // No. of sequences for each suffix array
        //---------------------------------------
        int n = (i < Param::nparts-1) ? 
            nreads/Param::nparts : 
            nreads - int(nreads/Param::nparts)*(Param::nparts-1);
        if (Param::verbose) std::cout << "#sequences to load:" << n << "\n";

        //------------------------------------------
        // Set the correct starting sequence pointer
        //------------------------------------------ 
        gsa[i].setSequences(seqs+s);

        size_t g = gsa[i].getSize();
        double tic = mytime();
        std::cout << "Generating read start count array ...\n";
        read_starts[i].init(gsa+i, g); 
        std::cout << "Read start count array generated:" << mytime()-tic << " sec\n";

        if ( Param::verbose ) gsa[i].setVerbosity(true);
        // read_ends[i].init(gsa+i, g, seqs+s); 
        // if ( Param::verbose > 1 ) {
        //     std::cout << "Read end counts\n";
        //     read_ends[i].print(gsa+i, std::cout);
        // }
    }

    //-----------------------
    // No. of assembled paths
    //-----------------------
    npaths = 0;

    //-------------------------------
    // Set assembled reads (NOT_PATH)
    //-------------------------------
    preads = new PathId[nreads];
    memset( preads, NOT_PATH, nreads*sizeof(PathId));

    //--------------------------------
    // Mark if any bad reads are given
    //--------------------------------
    bad_reads = loader->getBadReads();
    for ( size_t i = 0; i < bad_reads->size(); i++ )
        preads[(*bad_reads)[i]] = BAD_READ;

    //----------------------------
    // Construct a de Bruijn graph
    //----------------------------
    buildGraph();

    srch_log.ct_reads = nreads;
}

//====================================================================
// Construct a de Bruijn graph
//====================================================================
void Extracter::buildGraph()
{
    if ( Param::verbose ) std::cout << "\nBuilding graph ...\n";

    //----------------------------------------------------------------
    // This graph input file (gin) must exist in the working directory
    //----------------------------------------------------------------
    std::string graph_file = Param::out_dir + "/gin";

    //-------------------------------------------
    // Construct a graph by calling graph builder
    //-------------------------------------------
    graph::build(graph, kmer_coverage, vertex_map, Param::kmer_size, graph_file, Param::verbose);
    printElapsed( INIT_TIME, mytime(), "Graph built" );
    printMemoryUsage();

    //----------------------------------------
    // Display graph summary (nodes and edges)
    //----------------------------------------
    graph::graphSummary(graph, std::cout);

    //----------------------------
    // Display coverage statistics
    //----------------------------
    graph::coverageSummary(graph, kmer_coverage, vertex_map, std::cout);
}

//====================================================================
// Trim low supported nodes and edges
//====================================================================
void Extracter::trim()
{
    if ( ! Param::trim_flag ) return;

    if ( Param::verbose) std::cout << "\nTrimming graph ...\n";
    trimmed_nodes = graph::trimGraph(graph, kmer_coverage, vertex_map, Param::min_depth, Param::verbose );
    printElapsed(INIT_TIME, mytime(), "Graph trimmed");
	if ( Param::verbose ) graph::graphSummary(graph, std::cout);    

    //-------------------------------------------------
    // Now delete coverage information of trimmed nodes
    //-------------------------------------------------
    trimMap();

    graph::graphSummary(graph, std::cout);
    graph::coverageSummary(graph, kmer_coverage, vertex_map, std::cout);
}

//====================================================================
// Delete coverage information of deleted nodes
//====================================================================
void Extracter::trimMap( )
{
    //NodeSet::iterator it;
    //for ( it = dropped.begin(); it != dropped.end(); ++it ) {
    for ( auto item : trimmed_nodes ) {
        Vertex v = item.first;
        assert( vertex_map.find(v) != vertex_map.end() );
        KmerId kid = vertex_map[v];
        if ( Param::verbose ) std::cout << "dropped: " << alpha::IntegerToAminoAcid(kid,Param::kmer_size ) << "\n";         
        kmer_coverage.erase( kmer_coverage.find(vertex_map[v]) );
        trimmed_kmers[kid] = true;
    }
}

//====================================================================
// Run extractor
//====================================================================
void Extracter::run()
{
    //---------------------------------
    // Display read sequence statistics
    //---------------------------------
    stat();

    //------------------------
    // Trim graph if requested
    //------------------------
    trim();

    //--------------
    // Extract paths
    //--------------
    extract();

    //---------------------------------------
    // Release memories that no further needs
    //---------------------------------------
    release();

    //-----------------
    // Extractor ran OK
    //-----------------
    status = true;
}

//====================================================================
// Display read statistics
//====================================================================
void Extracter::stat()
{
    std::vector<double> lengths;
    lengths.reserve(nreads);
    for ( int i = 0; i < nreads; i++ ) {
        if ( preads[i] == BAD_READ ) continue;
        std::string str = seqs[i];
        lengths.push_back((double)str.size());
    }

    //-------------------------------------------
    // Sort read lengths for median and quantiles
    //-------------------------------------------
    std::sort( lengths.begin(), lengths.end() );

    size_t n = lengths.size();
    double max = math::max( &lengths[0], n );
    double min = math::min( &lengths[0], n );
    double avg = math::mean( &lengths[0], n );
    double med = math::median( &lengths[0], n, true );

    std::cout << "\nRead length summary\n";
    printf( "Max:%d ", (int)max );
    printf( "Min:%d ", (int)min );
    printf( "Mean:%.2f ", avg );
    printf( "Median:%.2f\n", med );
    printf( "Quantiles: ( ");
    for ( int i = 10; i<=100; i+=10 ) {
        printf("%d%%:%.2f ", i, math::quantile(&lengths[0], n, i/100.0, true) );
    }
    std::cout << ")\n\n";
            
}

//====================================================================
// Extract all paths
//====================================================================
void Extracter::extract()
{
    //-------------
    // No buffering
    //------------- 
    setbuf(stdout, NULL);

    //-------------------------------------------------------------
    // Set minimum depth of seed kmers with given percentile option
    //-------------------------------------------------------------
    setStartMinCoverage();

    //-------------------------------------------
    // Ancestor mapping
    // de Bruijn graph is directed graph.
    // Maintain ancestors map for quick retrieval
    //-------------------------------------------
    ancestors = graph::predecessorMap(graph);

    //-----------------------------------------
    // Make search order based on support value
    //-----------------------------------------
    initPathSearchPriority();

    //--------------
    // Extract paths
    //--------------
    double t = mytime();
    extractPaths(t);

    //---------------------------
    // Display extraction summary
    //---------------------------
    displaySummary();

    if ( Param::skip_fail && Param::try_later ) {
        makeFailedSeedsPriority();

        srch_failed.clear();
        read_failed.clear();
        save_failed.clear();

        extractPaths(t);

        displaySummary();
    }

    //--------------------
    // Write path if asked
    //--------------------
    writePath();
}

//====================================================================
// Find the minimum seed coverage if percentile option is given.
//====================================================================
void Extracter::setStartMinCoverage()
{
    if ( !Param::percent_flag ) return;
    size_t n = kmer_coverage.size();
    double *coverage = new double[n];

    int i = 0;
    CoverageMap::iterator it;
    for (  it = kmer_coverage.begin(); it != kmer_coverage.end(); ++it )
        coverage[i++] = (double)it->second;
    math::sort(coverage, n);  

    //-----------------------------------
    // Get seed value of given percentile
    //-----------------------------------
    double value = math::quantile(coverage, n, (100-Param::kmer_percentile)/100.0, 1);

    //------------------------
    // Set min seed value here
    //------------------------
    Param::min_seed = value;
    std::cout << "Min. coverage of seed kmers:" << Param::min_seed << "\n";

    //-------------
    // Release this
    //-------------
    delete[] coverage;
}

//====================================================================
// Set priorities of all seed kmers.
//====================================================================
void Extracter::initPathSearchPriority()
{
    Vertex_iter vc, ve;
    for ( boost::tie(vc,ve) = vertices(graph); vc != ve; ++vc ) {
        updatePathSearchPriority(*vc);
    }
    printElapsed( INIT_TIME, mytime(), "Seed kmers extracted" );
    printMemoryUsage();
}

void Extracter::makeFailedSeedsPriority()
{
    double t0 = mytime();

    //-------------------------------------
    // Drop all seed nodes for initial runs
    //-------------------------------------
    seed_nodes.clear();

    //----------------------------
    // Now collect new seeds nodes
    //----------------------------
    NodeFlags new_seeds;
    for ( auto item : srch_failed ) {
        if ( ! item.second ) continue;
        if ( used_seeds.find(item.first) != used_seeds.end() ) continue;
        new_seeds.insert(item);
    }
    for ( auto item : read_failed ) {
        if ( ! item.second ) continue;
        if ( used_seeds.find(item.first) != used_seeds.end() ) continue;
        new_seeds.insert(item);
    }
    for ( auto item : save_failed ) {
        if ( ! item.second ) continue;
        if ( used_seeds.find(item.first) != used_seeds.end() ) continue;
        new_seeds.insert(item);
    }

    //-------------------------------------
    // Generate seed priority for new seeds
    //-------------------------------------
    for ( auto item : new_seeds ) {
        Vertex v = item.first;
        if ( deleted_nodes.find(v) != deleted_nodes.end() ) continue;
        // if ( ! Param::use_all_fail ) {
        //     if ( ancestors[v].size() <= 1 && graph::successors(graph, v).size() <=1 )
        //         continue;
        // }
        updatePathSearchPriority(v);
    }

    //-------------------------------------------
    // Now clean previously used seed information
    //-------------------------------------------
    used_seeds.clear();

    std::cout << "Seeds regenerated:" << mytime()-t0 << " sec\n";
}

//====================================================================
// Set the priority of given seed node
//====================================================================
void Extracter::updatePathSearchPriority(Vertex v)
{    
    size_t sd = ancestors[v].size() +  graph::successors(graph, v).size();
    if ( sd == 0 ) return;

    CoverageType count = kmer_coverage[ vertex_map[v] ];
    if ( (int)count < Param::min_seed ) return;

    double weight = count / exp((double)sd);

    // 
    weight = round(weight*10)/10;

    //------------------------------------------------
    // Insert a seed node and its weight to seed table
    //------------------------------------------------
    seed_nodes.insert( v, (BinType)weight );
}

//====================================================================
// Extract all paths from seed nodes
//====================================================================
void Extracter::extractPaths(double lt0)
{
    double t0 = mytime();
    //double lt0 = mytime(), 
    double tic = mytime(), tbeg=mytime(), tend=mytime(), tpat=mytime();

    std::cout << "# Seed k-mers:" << seed_nodes.getSize() << "\n";
    std::cout << "# Seed bins:" << seed_nodes.getBinSize() << "\n";

    int nseeds = seed_nodes.getSize();
    progress = Progress( 1.0, 0, nseeds, mytime() );
    
    bool valid = false;
    while ( seed_nodes.getSize() ) {
        tbeg = mytime();
        srch_log.et_lapse += (tbeg-tend);
        srch_log.et_total = (mytime()-lt0);

        GraphPath gpat;
        ReadAligner raln;
        ReadAlignLog rlog;

        //---------------------------
        // Reset for a new extraction
        //---------------------------
        check_pid.clear();
        check_pos.clear();
        //check_starts.clear();
        //check_ends.clear();
        rhist_map.clear();
        //lhist_map.clear();
        
        valid = false;

        //------------------------------------------
        // Get a seed from seed table by popping one
        //------------------------------------------
        tic = mytime();
        Vertex seed = seed_nodes.pop();
        srch_log.et_init_seed += (mytime()-tic);
        if ( seed == NULL ) break;

        if ( Param::trim_flag && deleted_nodes.find(seed) != deleted_nodes.end() )
            continue;

        //-------------------------------------------
        // Get successors and predecessor of the seed
        //-------------------------------------------
        tic = mytime();
        NodeArray succs, preds;
        //assert( deleted_nodes.find(seed) == deleted_nodes.end());
        succs = graph::successors(graph, seed);
        preds = ancestors[seed];
        srch_log.et_init_graph += (mytime()-tic);


        tic = mytime();
        //----------------------------------
        // Get progress ratio before updated
        //----------------------------------
        double prev = progress.ratio;

        //-------------------
        // Set progress count
        //-------------------
        progress.count = nseeds - (int)seed_nodes.getSize();

        //--------------------------------------------------
        // Display progress when progress ratio is increased
        //--------------------------------------------------
        progress.showProgress();
        if ( progress.ratio > prev && Param::verbose >= 1 ) {
            srch_log.ct_assem = countAssembledReads();
            srch_log.printProgress(read_log);
            printf("# Saved GSA search entries (start reads ranges:%zu, end read ranges:%zu)\n", used_starts.size(), used_ends.size() );
            printMemoryUsage();
        }
        srch_log.et_init_prog += (mytime()-tic);
        srch_log.ct_init++;
        srch_log.et_init += (mytime()-tbeg);

        //--------------------------------------
        // Check whether it is valid search node
        //--------------------------------------
        tic = mytime();
        valid = validStartNode(seed);
        srch_log.et_seed += ( mytime() - tic );

        NodeFlags::iterator st = srch_failed.find(seed);
        NodeFlags::iterator rt = read_failed.find(seed);        
        NodeFlags::iterator ut = save_failed.find(seed);

        
        tic = mytime();
        if ( !valid ) {
            srch_log.ct_skip_seed++;
            srch_log.et_skip_seed += (mytime()-tic);
            srch_log.et_skip += (mytime()-tic);
            tend = mytime();
            continue;
        }
        srch_log.ct_good_seed++;    

        if ( Param::skip_fail ) used_seeds.insert( std::pair<Vertex, bool>(seed, 1 ) );

        //---------------------------------------
        // Now extract a path from the given seed
        //---------------------------------------
        tpat = mytime();
        valid = extractOnePath(seed, gpat);
        srch_log.et_path += ( mytime()-tic );

        tic = mytime();
        if ( !valid ) {
            if ( Param::verbose ) std::cout << "Extraction fail\n";
            srch_log.et_skip += (mytime()-tic);

            if ( ut != save_failed.end() ) srch_log.ct_path_fail_save++;
            if ( st != srch_failed.end() ) srch_log.ct_path_fail_srch++;
            if ( rt != read_failed.end() ) srch_log.ct_path_fail_read++;

            if ( st != srch_failed.end() || rt != read_failed.end() || ut != save_failed.end() ) {
                ( succs.size() <= 1 && preds.size() <= 1 ) ? 
                    srch_log.ct_path_sim_node_fail++ :
                    srch_log.ct_path_mul_node_fail++;
            }
            
            addFlags(gpat.path, srch_failed);

            size_t plen = gpat.path.size() ? gpat.path.size()+Param::kmer_size-1 : 0;
            if ( srch_log.srch_fail_min > plen ) srch_log.srch_fail_min = plen;
            if ( srch_log.srch_fail_max < plen ) srch_log.srch_fail_max = plen;
            srch_log.srch_fail_sum += plen;
            srch_log.srch_fail++;
            
            tend = mytime();
            
            srch_log.et_srch_fail += (mytime()-tpat);
            continue;
        }

        if ( ut != save_failed.end() ) srch_log.ct_path_succ_save++;
        if ( st != srch_failed.end() ) srch_log.ct_path_succ_srch++;
        if ( rt != read_failed.end() ) srch_log.ct_path_succ_read++;

        if ( st != srch_failed.end() || rt != read_failed.end() || ut != save_failed.end() ) {
            ( succs.size() <= 1 && preds.size() <= 1 ) ? 
                srch_log.ct_path_sim_node_succ++ :
                srch_log.ct_path_mul_node_succ++;
        }

        //-------------------------------------------------
        // Check whether the extracted path is valid or not
        //-------------------------------------------------
        valid = checkPath(gpat, raln, rlog);
        srch_log.et_align += ( mytime()-tic );
        srch_log.ct_align++;

        tic = mytime();
        if ( valid ) {
            srch_log.ct_good_align++;
            if ( ut != save_failed.end() ) srch_log.ct_align_succ_save++;
            if ( st != srch_failed.end() ) srch_log.ct_align_succ_srch++;
            if ( rt != read_failed.end() ) srch_log.ct_align_succ_read++;

            if ( st != srch_failed.end() || rt != read_failed.end() || ut != save_failed.end() ) {
                ( succs.size() <= 1 && preds.size() <= 1 ) ? 
                    srch_log.ct_align_sim_node_succ++ :
                    srch_log.ct_align_mul_node_succ++;
            }

        }
        else {
            if ( Param::verbose ) std::cout << "Path check fail\n";
            srch_log.et_skip += (mytime()-tic);

            if ( ut != save_failed.end() ) srch_log.ct_align_fail_save++;
            if ( st != srch_failed.end() ) srch_log.ct_align_fail_srch++;
            if ( rt != read_failed.end() ) srch_log.ct_align_fail_read++;

            if ( st != srch_failed.end() || rt != read_failed.end() || ut != save_failed.end() ) {
                ( succs.size() <= 1 && preds.size() <= 1 ) ? 
                    srch_log.ct_align_sim_node_fail++ :
                    srch_log.ct_align_mul_node_fail++;
            }

            addFlags(gpat.path, read_failed);

            size_t plen = gpat.path.size() ? gpat.path.size()+Param::kmer_size-1 : 0;
            if ( srch_log.align_fail_min > plen ) srch_log.align_fail_min = plen;
            if ( srch_log.align_fail_max < plen ) srch_log.align_fail_max = plen;
            srch_log.align_fail_sum += plen;
            srch_log.align_fail++;

            tend = mytime();

            srch_log.et_align_fail += (mytime()-tpat);
            continue;
        }
        
        //-------------------------------------------------
        // Save the path and update working data structures
        //-------------------------------------------------
        tic = mytime();
        valid = update(gpat, raln);
        if ( !valid && Param::verbose ) std::cout << "Update fail\n";
        srch_log.et_update += (mytime()-tic);
        valid ? srch_log.ct_update_succ++ : srch_log.ct_update_fail++;
        if ( valid ) {
            if ( ut != save_failed.end() ) srch_log.ct_update_succ_save++;
            if ( st != srch_failed.end() ) srch_log.ct_update_succ_srch++;
            if ( rt != read_failed.end() ) srch_log.ct_update_succ_read++;

            if ( st != srch_failed.end() || rt != read_failed.end() || ut != save_failed.end() ) {
                ( succs.size() <= 1 && preds.size() <= 1 ) ? 
                    srch_log.ct_update_sim_node_fail++ :
                    srch_log.ct_update_mul_node_fail++;
            }
            srch_log.et_save_succ += (mytime()-tpat);

            size_t plen = gpat.path.size() ? gpat.path.size()+Param::kmer_size-1 : 0;
            if ( srch_log.save_succ_min > plen ) srch_log.save_succ_min = plen;
            if ( srch_log.save_succ_max < plen ) srch_log.save_succ_max = plen;
            srch_log.save_succ_sum += plen;
            srch_log.save_succ++;

        } else {
            if ( ut != save_failed.end() ) srch_log.ct_update_fail_save++;
            if ( st != srch_failed.end() ) srch_log.ct_update_fail_srch++;
            if ( rt != read_failed.end() ) srch_log.ct_update_fail_read++;

            if ( st != srch_failed.end() || rt != read_failed.end() || ut != save_failed.end() ) {
                ( succs.size() <= 1 && preds.size() <= 1 ) ? 
                    srch_log.ct_update_sim_node_fail++ :
                    srch_log.ct_update_mul_node_fail++;
            }


            addFlags(gpat.path, save_failed);

            size_t plen = gpat.path.size() ? gpat.path.size()+Param::kmer_size-1 : 0;
            if ( srch_log.save_fail_min > plen ) srch_log.save_fail_min = plen;
            if ( srch_log.save_fail_max < plen ) srch_log.save_fail_max = plen;
            srch_log.save_fail_sum += plen;
            srch_log.save_fail++;
            
            srch_log.et_save_fail += (mytime()-tpat);
        }

        if ( Param::verbose ) printf("Update complete:%.4f\n", mytime()-t0);
        
        if ( Param::verbose ) srch_log.printOneLog();

        tend = mytime();
    }
}

//====================================================================
// Check whether the given seed node is valid
//====================================================================
bool Extracter::validStartNode( Vertex &v )
{
    //-------------------------------------
    // Case: only unique seeds are allowed.
    //-------------------------------------
    //if ( Param::uniq_seed ) {
    if ( ! Param::seed_reuse ) {
        if ( discovered.find(v) != discovered.end() ) {
            return false;
        }
    }

    if ( Param::trim_flag ) {
        if ( deleted_nodes.find(v) != deleted_nodes.end() ) 
            return false;
    }

    // island
    if ( out_degree(v, graph) == 0  ) {
        AdjacencyMap::iterator it = ancestors.find(v);
        if ( it == ancestors.end() || it->second.size() == 0 ) 
            return false;
    }

    if ( Param::right_ext_only && ! graph::isSource( v, ancestors ) ) 
        return false;

    //----------------------------------------------------------------------
    // Case: only previously succeeded seeds are allowed if a seed is reused
    //----------------------------------------------------------------------
    if ( Param::skip_fail ) {
        if ( srch_failed.find(v) != srch_failed.end() ) 
            return false;
        if ( read_failed.find(v) != read_failed.end() ) 
            return false;
        if ( save_failed.find(v) != save_failed.end() ) 
            return false;
    }

    //------------------------------------------
    // Shouldn't occur (No vertex mapping exist)
    //------------------------------------------
    if ( vertex_map.find(v) == vertex_map.end() ) { 
        std::cerr << "[Error] Vertex map for :" << v << " does not exit\n"; 
        exit(EXIT_FAILURE);
    }
    
    //--------------------
    // Current pivot k-mer
    //--------------------
    KmerId kid = vertex_map[v];
    
    //----------------------
    // Invalid start kmer
    // Start with stop codon
    //----------------------
    if ( alpha::getFirstAminoAcid<KmerId>(kid, Param::kmer_size) == '*') {
        if ( Param::verbose )  std::cout << "** Invalid Start Kmer:" << alpha::IntegerToAminoAcid(kid, Param::kmer_size ) << "\tskip\n"; 
        return false;
    }
    
    //----------------------------
    // Shouldn't occur
    // Coverage value is not found
    //----------------------------
    if ( kmer_coverage.find(kid) == kmer_coverage.end() ) {
        if ( Param::verbose ) std::cerr << "** Coverage does not exist: " << alpha::IntegerToAminoAcid(kid,Param::kmer_size ) << "\n"; 
        exit(EXIT_FAILURE);
    }

    //------------------
    // Low coverage seed
    //------------------
    CoverageType size = kmer_coverage[kid];
    if ( (int)size < Param::min_seed ) {
        if ( Param::verbose ) std::cout << "** Low read coverage:" << alpha::IntegerToAminoAcid(kid,Param::kmer_size ) << "\tcoverage:" << size << "\tskip\n";
        return false;
    }
    
    //-------------------------------------
    // Shouldn't occur (No predecessor map)
    //-------------------------------------
    if ( ancestors.find(v) == ancestors.end() ) {
        if ( Param::verbose ) std::cerr << "** Predecessor map not defined: " << alpha::IntegerToAminoAcid(kid,Param::kmer_size ) << "\n"; 
        exit(EXIT_FAILURE);
    }

    //--------------------
    // Should be ok by now
    //--------------------
    return true;
}

//====================================================================
// Extract one path starting from the given node v
//====================================================================
bool Extracter::extractOnePath( Vertex v, GraphPath &gpath )
{
    double t0 = mytime();
    double tic;

    //-----------------------
    // Reset trace string map
    //-----------------------
    traces_map.clear();

    //--------------------
    // Current pivot k-mer
    //--------------------
    assert( vertex_map.find(v) != vertex_map.end() );
    KmerId kid = vertex_map[v];
        
    //--------------------
    // Depth of start kmer
    //--------------------
    assert( kmer_coverage.find(kid) != kmer_coverage.end() );
    size_t start_depth = kmer_coverage[kid];

    //-----------------------------
    // Integer to string conversion
    //-----------------------------
    std::string seed_str = alpha::IntegerToAminoAcid(kid,Param::kmer_size );
    if ( Param::verbose ) {
        std::cout << "\n*Start:" << seed_str 
                  << "\tcoverage:" << start_depth << "\t"
                  << "#indeg:" << ancestors[v].size() << "\t" 
                  << "#outdeg:" << graph::successors(graph, v).size() << "\n";
    }

    /*===========================================*/
    /* Discover path with the current start node */
    /*===========================================*/
    tic = mytime();
    
    //----------------------------
    // Set a seed of a path object
    //----------------------------
    gpath.seed = kid;

    //-----------------------------------------------------------------------------
    // Find the seed string from suffix arrays. 
    // Save suffix array(s) result(s) which is left-most and right-most boundaries
    //-----------------------------------------------------------------------------
    BoundArray bounds(Param::nparts, BoundType());
    for ( int i = 0; i < Param::nparts; i++ ) {
        double tic = mytime();
        BoundType srch = gsa[i].search((SfaChar*)seed_str.c_str(), seed_str.size());
        srch_log.et_array += (mytime()-tic);
        srch_log.ct_array++;
        
        bounds[i] = srch;
        if ( Param::verbose ) {
            std::cout << "(" << srch.first << "\t" << srch.second << ")\n";
        }   
    }
    // srch_log.et_array += (mytime()-t0);
    // srch_log.ct_array++;
    
    //-------------------------
    // Add suffix array results
    //------------------------- 
    gpath.bounds.push_back(bounds);

    //---------------------------------------------
    // Add the trace string size
    // Here, the trace string is seed string itself
    //---------------------------------------------
    gpath.traces.push_back(seed_str.size());

    //----------------------------------------
    // Insert seed string to traces string map
    //----------------------------------------
    traces_map.insert( std::pair<std::string, bool>(seed_str, true) );
    srch_log.et_path_init += (mytime()-t0);

    //-----------------------------------------------
    // Start greedy path extension from the seed node
    //-----------------------------------------------
    tic = mytime();
    bool good = extractGreedyBestPath(v, gpath);
    srch_log.et_path_path += (mytime()-tic);

    tic = mytime();
    /*=================*/
    /* Path Validation */
    /*=================*/
    if ( !good ) {
        srch_log.ct_poor_path++;
        srch_log.et_poor_path += (mytime()-t0);
        srch_log.et_path_post += (mytime()-tic);
        return false;
    }
    
    //----------------------------------------------
    // Check whether the path is too short to report
    //----------------------------------------------
    if ( shortPath(gpath.path) ) {
        if ( Param::strict_length ||
             ( gpath.lstop != PATHEND || gpath.rstop != PATHEND) ) {
            srch_log.ct_short_path++;
            srch_log.ct_poor_path++;
            srch_log.et_poor_path += (mytime()-t0);
            srch_log.et_path_post += (mytime()-tic);
            return false;
        }
    }
    
    //--------------------------------------
    // Now get numeric kmers from path nodes
    //--------------------------------------
    setPathKmers(gpath);
    srch_log.et_path_post += (mytime()-tic);

    srch_log.ct_good_path++;
    srch_log.et_good_path += (mytime()-t0);
    
    return true;
}

//====================================================================
// Check whether a path is valid or not.
//====================================================================
bool Extracter::checkPath( GraphPath &gpath, ReadAligner &raln, ReadAlignLog &rlog )
{
    if ( Param::verbose ) std::cout << "Path:" << gpath.toString(Param::kmer_size) << "\n";

    //--------------------
    // Align reads to path
    //--------------------
    raln = ReadAligner(&gpath, gsa, read_starts, seqs, nreads, preads, &rlog);
    read_log.add(rlog);
    
    //--------------
    // Enough reads?
    //--------------
    if ( raln.getSize() < (size_t)Param::min_share ) {
        srch_log.ct_check_fail_read++;
        if ( Param::verbose ) std::cout << "Small reads placement\n";
        return false;
    }

    //----------------------------
    // Check read alignment status
    //---------------------------- 
    if ( raln.getStatus() == false ) {
        srch_log.ct_check_fail_trim++;
        if ( Param::verbose ) std::cout << "Read alignment fail\n";
        return false;
    }

    //------------------
    // Check path length
    //------------------
    if ( ! Param::short_trim_len ) {
        if ( ! validLengthPath(raln) ) {
            srch_log.ct_check_fail_len++;
            if ( Param::verbose ) std::cout << "Invalid length:" << raln.getReference().size() << "\n";
            return false;
        }
    }

    //---------------------------------------------------------
    // Trim nodes if the path is shortened after read alignment
    //---------------------------------------------------------
    trimPath(gpath,raln);

    if ( Param::verbose ) 
        std::cout << "Path:" << raln.getReference() 
                  << "\tstops:" << raln.getLStop() << " " << raln.getRStop() 
                  << "\ttrims:" << raln.getLTrim() << " " << raln.getRTrim() << "\n";

    return true;
}

//====================================================================
// Make numeric kmers from path nodes and save them to path object
//====================================================================
void Extracter::setPathKmers( GraphPath &p )
{
    KmerArray kids = KmerArray( p.path.size(), 0 );
    int i = 0;
    for ( PathType::iterator pt = p.path.begin(); pt != p.path.end(); ++pt )
        kids[i++] = vertex_map[*pt];
    p.kmers = kids;
}


//====================================================================
// Check whether the path is too short
//====================================================================
bool Extracter::shortPath( PathType &max_path )
{
    std::string path_str = biostr::getSequenceString(max_path, vertex_map, Param::kmer_size);
    if ( (int)path_str.length() < Param::min_length ) {
        if (Param::verbose) std::cout << "\nShort Path\n";
        return true;
    }
    return false;
}

//====================================================================
// Extract a path from a given node by greedy extension 
//====================================================================
bool Extracter::extractGreedyBestPath( Vertex curr,
                                       GraphPath &gpath )
{
    double t0 = mytime();

    assert( ancestors.find(curr) != ancestors.end() );
    NodeArray preds = ancestors[curr];
    NodeArray succs = graph::successors(graph, curr);


    //------------------------------
    // append current node to a path
    //------------------------------
    gpath.path.push_back(curr); 

    //----------------------
    // seed position in path
    //----------------------
    gpath.spos = 0;    
   
    assert( vertex_map.find(curr) != vertex_map.end() );
    KmerId kid = vertex_map[curr];    

    std::string path_str = alpha::IntegerToAminoAcid(kid, Param::kmer_size);
    int reason;
    double tic;

    if ( Param::verbose ) 
        std::cout << "(+) "
                  << alpha::IntegerToAminoAcid(kid, Param::kmer_size) << "\n";
    
    //-------------
    // trace string
    //-------------
    // double tic = mytime();
    //std::string nback_str = alpha::IntegerToAminoAcid(kid, Param::kmer_size);
    // srch_log.et_node_nback += (mytime()-tic);
    // srch_log.et_node_trace += (mytime()-tic);
    
    //------------------
    // Do left extension
    //------------------
    tic = mytime();
    
    //std::string path_str = nback_str;
    
    //extendNode( curr, path_str, nback_str, gpath, LEFT, reason );
    extendNode( curr, path_str, gpath, RIGHT, reason );
    srch_log.et_path_extend_right += (mytime()-tic);
    srch_log.et_path_extend += (mytime()-tic);
    if (Param::verbose) std::cout << "Right Extension time:" << mytime()-tic << "\n";
    
    //-------------------------------
    // Set left extension stop reason
    //-------------------------------
    tic = mytime();
    gpath.rstop = reason;
    
    //----------------------------------------------------------
    // In case of zero length path after duplicated path removal
    //----------------------------------------------------------
    if ( reason == CYCLICPATH ) {
        //if ( Param::drop_dup_path ) 
            if ( gpath.path.size() == 0 ) return false;
    }
    
    
    
    //---------------------------------
    // extension in the other direction
    //---------------------------------
    if ( Param::verbose ) 
        std::cout << "(-) "
                  << alpha::IntegerToAminoAcid(kid, Param::kmer_size) << "\n";
    
    //-----------------------------------------
    // Initialize trace str for right extension
    //-----------------------------------------
    tic = mytime();
    //PathType nback_nodes = getSubPath( gpath.path, RIGHT, Param::back_trace);
    //nback_str  = biostr::getSequenceString(nback_nodes, vertex_map, Param::kmer_size); 
    srch_log.et_node_nback += (mytime()-tic);
    srch_log.et_node_trace += (mytime()-tic);

    check_pid.clear();
    check_pos.clear();
    //check_starts.clear();
    //check_ends.clear();

    // if ( nback_str.size() <= (size_t)Param::back_trace && 
    //      nback_str.size() >= (size_t)Param::kmer_size+2 ) {
    if ( path_str.size() <= (size_t)Param::back_trace && 
         path_str.size() >= (size_t)Param::kmer_size+2 ) {
        //std::string nstr = nback_str.substr(0, Param::kmer_size+2);
        std::string nstr = path_str.substr(0, Param::kmer_size+2);
        findTraceStringFromPaths( nstr, check_pos, check_pid );
    }


    //-------------------
    // Do right extension
    //-------------------
    path_str = biostr::getSequenceString(gpath.path, vertex_map, Param::kmer_size); 
    //extendNode( curr, path_str, nback_str, gpath, RIGHT, reason );
    extendNode( curr, path_str, gpath, LEFT, reason );
    srch_log.et_path_extend_left += (mytime()-tic);
    srch_log.et_path_extend += (mytime()-tic);
    if(Param::verbose) std::cout << "Left Extension time:" << mytime()-tic << "\n";
    tic = mytime();
    gpath.lstop = reason;

    //int spos = gpath.spos;

    if ( reason == CYCLICPATH ) {
        //if ( Param::drop_dup_path ) {
            //-----------------
            // Zero length path
            //-----------------
            if ( gpath.path.size() == 0 ) return false;

            //---------------------------------------------
            // Seed node was lost because of repeat removal
            //---------------------------------------------
            if ( (int)gpath.path.size() <= gpath.spos ) {
                if ( Param::verbose ) {
                    std::cout << "Seed is removed\n";
                    std::cout << "Seed:" << gpath.spos << "\tPath size:" << gpath.path.size() << "\n";
                }
                return false;
            }
            //}
    }
    
    if ( Param::verbose ) 
        std::cout << "\n** Best path search (length:" 
                  << gpath.path.size()+Param::kmer_size-1 
                  << "\ttime:" << mytime() - t0 << " sec)\n";

    return true;
}

//====================================================================
// Check whether the given node is a terminal node
//====================================================================
bool Extracter::isTerminalNode( Vertex curr,
                                GraphPath &gpath,
                                int direction,
                                int &reason )
{
    KmerId cur_kid = vertex_map[curr];

    //----------------------------------
    // Check whether it is a source node
    //----------------------------------
    if ( graph::isSource(curr, ancestors) && direction == LEFT ) {
        reason = PATHEND;
        if ( Param::verbose ) std::cout << "** Path ends\n"; 
        return true;
    }
    
    //--------------------------------
    // Check whether it is a sink node
    //--------------------------------
    if ( graph::isSink(curr, graph) && direction == RIGHT ) {
        reason = PATHEND;
        if ( Param::verbose ) std::cout << "** Path ends\n"; 
        return true;
    }

    //------------------------------------------
    // Check wether it is a node with stop codon
    //------------------------------------------
    char ch = alpha::getLastAminoAcid<KmerId>( cur_kid, Param::kmer_size );
    if ( ch == '*' ) {
        if ( direction == RIGHT ) { 
            if ( Param::verbose ) std::cout << "** Stop codon\n";
            reason = PATHEND;
            return true;
        }
    }
    
    return false;
}

bool Extracter::findMaxNeighbor( Vertex curr,
                                 std::string &path_str,
                                 //std::string &nback_str,
                                 GraphPath &gpath,
                                 int direction,
                                 int &reason,
                                 Vertex &max_node,
                                 BoundArray &max_bound,
                                 std::string &max_trace,
                                 UsedBoundInfoList &max_useds )
{
    double tic = mytime();
    NodeArray vset = getNeighboringNodes(curr, direction);
    srch_log.et_node_ngb += (mytime()-tic);
    srch_log.et_next_ngb += (mytime()-tic);
    if ( Param::verbose ) std::cout << "Neighboring node search:" << mytime()-tic << " sec\n";
    
    
    tic = mytime();
    if ( vset.size() == 0 ) {
        reason = EXTENDFAIL;
        if ( Param::verbose ) std::cout << "All weakly supported neighbors\n"; 
        return false; 
    }
    

    tic = mytime();
    int max_value = 0;
    findGreedyBestNeighbor( gpath, vset, curr, path_str, max_node, max_value, max_bound, max_trace, max_useds, direction);

    direction == LEFT ? 
        srch_log.ct_path_extend_left++ : 
        srch_log.ct_path_extend_right++;

    srch_log.et_node_max += (mytime()-tic);
    srch_log.et_next_max += (mytime()-tic);

    if ( Param::verbose ) std::cout << "Maximal node search:" << mytime()-tic << " sec\n";
        
    tic = mytime();
    if ( max_node == NULL ) { 
        if ( Param::verbose ) std::cout << "** Extension fail\n"; 
        reason = EXTENDFAIL;
        return false; 
    }
    
    tic = mytime();
    char ch = alpha::getFirstAminoAcid<KmerId>(vertex_map[max_node], Param::kmer_size);
    if ( ch == '*') {
        if ( Param::verbose ) std::cout << "Bad extension\n"; 
        reason = EXTENDFAIL;
        return false;
    }

    if ( Param::verbose ) {
        std::cout << "Max k-mer:" 
                  << alpha::IntegerToAminoAcid(vertex_map[max_node], Param::kmer_size)  
                  << "\t# support:" << max_value << "\n";
    }
    
    tic = mytime();
    if ( path_str.size() < (size_t)Param::back_trace ) {
        if ( max_value < Param::min_share ) {
            if ( Param::verbose ) std::cout << "\n** Weak read support - max:" << max_value << "\n";
            reason = WEAKMAX;
            return false; 
        }
    } else {
        if ( max_value == 0 ) {
            if ( Param::verbose ) std::cout << "\n** Weak read support - max:" << max_value << "\n";
            reason = WEAKMAX;
            return false; 
        }
    }


    return true;
}

// bool Extracter::correctableRepeat( Vertex curr,
//                                    std::string &nback_str,
//                                    Vertex max_node,
//                                    GraphPath &gpath,
//                                    CyclicPath &cpath,
//                                    int direction,
//                                    int &reason )
// {
//     double tic = mytime();
//     bool succ = handleRepeat( nback_str, gpath, cpath, curr, max_node, direction );
//     srch_log.et_node_cyc += (mytime()-tic);
//     if ( !succ ) {
//         reason = CYCLICPATH;
//         return false;
//     }
//     if ( Param::verbose ) std::cout << "Repeat check:" << mytime() - tic << " sec\n";
//     return true;
// }

void Extracter::extendNode( Vertex curr,
                            std::string &path_str,
                            //std::string &nback_str,
                            GraphPath &gpath,
                            //CyclicPath &cpath,
                            int direction,
                            int &reason )
    
{
    while(true) {
        double tic = mytime();
        double t0 = mytime();
        //srch_log.ct_see++;
        if ( Param::trim_flag )
            assert( deleted_nodes.find(curr) == deleted_nodes.end());

        if ( Param::verbose ) {
            std::cout << "--------------------------------------------------\n"
                      << "Best:" 
                      << biostr::getSequenceString(gpath.path, vertex_map, Param::kmer_size) 
                      << "\n";
        }
        // srch_log.et_add += mytime()-t0;

        Vertex max_node = NULL;
        BoundArray max_bound;
        std::string max_trace;
        UsedBoundInfoList max_useds;        
        srch_log.et_path_pre += (mytime()-tic);

        // check whether the path is already discovered
        //if ( Param::dupcheck_flag && isRedundantPath( gpath, direction, reason ) ) return;

        // check whether the current node is the terminal node (sink/source/stop codon)
        //if ( isTerminalNode( curr, gpath, direction, reason ) ) return;
        tic = mytime();
        bool term = isTerminalNode( curr, gpath, direction, reason );
        srch_log.et_path_term += (mytime() - tic );
        if ( term ) {
            srch_log.ct_path_end ++;
            return;
        }

        // get the next maximal node
        //if ( !findMaxNeighbor( curr, gpath, direction, reason, max_node, max_bound ) ) return;
        tic = mytime();
        // bool maxn = findMaxNeighbor( curr, path_str, nback_str, gpath, direction, reason, max_node, max_bound, max_trace );
        bool maxn = findMaxNeighbor( curr, path_str, gpath, direction, reason, max_node, max_bound, max_trace, max_useds );
        srch_log.et_path_maxn += (mytime() - tic );
        if ( !maxn ) {
            srch_log.ct_path_weak++;            
            return;
        }
        // clean up the repeat
        //if ( !correctableRepeat( curr, max_node, gpath, cpath, direction, reason ) ) return;
        tic = mytime();
        //bool simple = correctableRepeat( curr, nback_path, nback_str, max_node, gpath, cpath, direction, reason );
        //bool cycle = isCyclic( gpath, nback_str, max_bound, direction );
        bool cycle = false;
        if ( (int)max_trace.size()==Param::back_trace+1 &&
             traces_map.find(max_trace) != traces_map.end() ) {
            if ( Param::verbose ) std::cout << "Cyclic path\n";
            //traces_map.erase(max_trace);

            if ( Param::check_cycle_entry ) {
                if ( fromSameNode(max_trace, gpath.path, direction) ) {
                    if ( Param::verbose ) std::cout << "Same entry cycle\n";
                    cycle = true;
                }
            } 
            else cycle = true;
        }
        srch_log.ct_path_cycle += (mytime()-tic);

        if ( cycle ) {
            reason = CYCLICPATH;
            srch_log.ct_path_cycle ++;
            
            //if ( Param::drop_dup_path ) {
                removeSubPath( gpath, direction );
                if ( gpath.path.size() == 0 ) return;
                
                std::string pstr = biostr::getSequenceString( gpath.path, vertex_map, Param::kmer_size);
                path_str = pstr;
                if ( (int)pstr.size() <= Param::back_trace ) {
                    traces_map.clear(); return;
                }
                TraceStringMap new_traces;
                for ( size_t k = 0; k < pstr.size()-(Param::back_trace+1)+1; k++ ){
                    std::string sstr = pstr.substr(k, Param::back_trace+1);
                    if ( traces_map.find(sstr) != traces_map.end() )
                        new_traces.insert( std::pair<std::string, bool>(sstr, true) );
                }
                traces_map = new_traces;
                //}
            return;
        }
        //if ( !simple ) return;

        // prepare for the next extesion
        tic = mytime();
        size_t trace_size;
        if ( path_str.size() >= (size_t)Param::back_trace+1 ) trace_size = Param::back_trace+1;
        else trace_size = path_str.size();
        if ( direction == LEFT ) {
            gpath.path.push_front(max_node);
            gpath.bounds.push_front(max_bound);
            gpath.traces.push_front(trace_size);
            gpath.used_ends.push_front(max_useds);
        } else {
            gpath.path.push_back(max_node);
            gpath.bounds.push_back(max_bound);
            gpath.traces.push_back(trace_size);
            gpath.used_begs.push_back(max_useds);
            //max_trace.size());
        }
        
        if ( Param::verbose > 1 ) std::cout << "Used begs size:" << gpath.used_begs.size() << "\n";
        if ( Param::verbose > 1 ) std::cout << "Used ends size:" << gpath.used_ends.size() << "\n";

        traces_map.insert( std::pair<std::string, bool>(max_trace, true) );
        

        // assert( nback_str.size() <= (size_t)Param::back_trace );
        // if ( direction == LEFT ) {
        //     // nback_path.push_front(max_node);
        //     if ( nback_str.size() < (size_t)Param::back_trace ) {
        //         nback_str = max_trace;
        //     } else {
        //         //nback_path.pop_back();
        //         nback_str = max_trace.substr(0, Param::back_trace);
        //     }
        // } else {
        //     //nback_path.push_back(max_node);
        //     if ( nback_str.size() < (size_t)Param::back_trace ) {
        //         nback_str = max_trace;
        //     } else {
        //         //nback_path.pop_front();
        //         nback_str = max_trace.substr(1);
        //     }
        // }
        
        // if ( Param::verbose ) std::cout << "New nback-str:" << nback_str << "\n";
        
        

        if ( direction == LEFT ) gpath.spos++;        
        if ( Param::verbose ) std::cout << "One node extension:" << mytime()-t0 << " sec\n";
        
        assert(vertex_map.find(max_node) != vertex_map.end());
        KmerId max_kid = vertex_map[max_node];
        char ch = direction == LEFT ? 
            alpha::getFirstAminoAcid(max_kid, Param::kmer_size):
            alpha::getLastAminoAcid(max_kid, Param::kmer_size);
        
        if ( direction == LEFT )
            path_str = ch + path_str;
        else
            path_str += ch;

        curr = max_node;

        srch_log.et_path_next += (mytime()-tic);
        // srch_log.et_add += (mytime()-tic);
    }
}


bool Extracter::fromSameNode( std::string &max_trace,
                              PathType &path,
                              int direction )
{
    std::string path_str = biostr::getSequenceString( path, vertex_map, Param::kmer_size );
    if ( path_str.size() <= (size_t)Param::back_trace ) return false;

    std::vector<size_t> poss;
    size_t found = path_str.find( max_trace );
    while ( found != std::string::npos ) {
        poss.push_back(found);
        found = path_str.find(max_trace, found+1);
    }
    assert(poss.size()>0);
    //assert( found != std::string::npos );
    //if ( count == 0 ) return false;

    for ( size_t i = 0; i < poss.size(); i++ ) {
        found = poss[i];
        if ( Param::verbose ) std::cout << "found:" << found << "\t#lenth:" << path.size() << "\n";
        if ( direction == LEFT ) {
            if ( found+Param::back_trace+1 == path_str.size() ) {
                if ( Param::verbose ) std::cout << "Sequence end\n";
                return true;
            }
            char prev = path_str[found+Param::back_trace+1];
            char curr = path_str[Param::back_trace];
            if ( Param::verbose ) 
                std::cout << "curr:" << curr << "\tprev:" << prev << "\n";
            if ( curr == prev ) return true;
        }
        else {
            if ( found == 0 ) {
                if ( Param::verbose ) std::cout << "Sequence end\n";
                return true;
            }
            char prev = path_str[found-1];
            char curr = path_str[path_str.size()-(Param::back_trace+1)];
            if ( Param::verbose ) 
                std::cout << "curr:" << curr << "\tprev:" << prev << "\n";
            if ( curr == prev ) return true;
        }   
    }     
    return false;
}

PathType Extracter::getSubPath( PathType &path,
                                int direction,
                                int sub_length )
{
    int length = (int)path.size()+Param::kmer_size-1;
    if ( length <= sub_length ) return path;
    
    int nnode =  sub_length - Param::kmer_size + 1;
    NodeArray nodes = NodeArray( path.begin(), path.end() );
    NodeArray trace = direction==LEFT ? 
        NodeArray(nodes.begin(),  nodes.begin()+nnode) :
        NodeArray(nodes.rbegin(), nodes.rbegin()+nnode);

    return direction==LEFT ? 
        PathType( trace.begin(), trace.end() ) :
        PathType( trace.rbegin(), trace.rend() );
}


NodeArray Extracter::getNeighboringNodes( Vertex &curr,
                                          int direction )
{
    double t0 = mytime();
    NodeArray nexts = direction == LEFT ?  
        ancestors[curr] : graph::successors(graph, curr);

    NodeArray nset = findGoodNeighbors( curr, nexts, direction );
    
    srch_log.ct_extend++;
    srch_log.et_graph += (mytime()-t0);
    srch_log.ct_graph += nset.size();
    
    return nset;
}

NodeArray Extracter::findGoodNeighbors( Vertex curr_node, 
                                        NodeArray &nexts,
                                        int direction )
{
    //double t0 = mytime();

    //NodeSet vset;
    NodeArray vset;
    for ( NodeArray::iterator it = nexts.begin(); it != nexts.end(); ++it )  {
        Vertex v = *it;
        if ( Param::trim_flag ) 
            assert( deleted_nodes.find(v) == deleted_nodes.end() );

        if ( vertex_map.find(v) == vertex_map.end() ) {
            std::cerr << "\n[ERROR] " << v << "\tVertex map does not exist\n";
            exit(EXIT_FAILURE);
        }

        if ( kmer_coverage[vertex_map[v]] < (unsigned)Param::min_share ) {
            if ( Param::verbose ) 
                std::cout << "\t" << alpha::IntegerToAminoAcid(vertex_map[v], Param::kmer_size) << "\t" << kmer_coverage[vertex_map[v]] << "\tlow coverage -> skip\n";
            continue;
        }
        

        std::pair<Edge, bool> ep;
        if ( direction == LEFT ) 
            ep = edge(v, curr_node, graph) ;
        else
            ep = edge(curr_node, v, graph) ;

        if ( ep.second == false ) {
            std::cerr << "\n[ERROR] Edge to " << alpha::IntegerToAminoAcid(vertex_map[v], Param::kmer_size) << " not exist\n";
            exit(EXIT_FAILURE);
        }

        EdgeWeight weight;
        direction == LEFT ? weight = graph[edge(v, curr_node, graph).first].weight : weight = graph[edge(curr_node, v, graph).first].weight;
        if ( weight < Param::min_share ) {
            if ( Param::verbose ) std::cout << "\t" << alpha::IntegerToAminoAcid(vertex_map[v], Param::kmer_size) << " : weak edge support (" << weight << ")\n";
            continue;
        }
        vset.push_back(v);
    }
    return vset;
}


void Extracter::orderNodes( NodeArray &vset,
                            NodeArray &nodes,
                            std::vector<size_t> &covers )
{
    //double t0 = mytime();

    //------------------------------
    // Neighbor node search priority
    //------------------------------
    std::multimap<size_t, Vertex> depths;
    for ( NodeArray::iterator it = vset.begin(); it != vset.end(); ++it )  
        depths.insert( std::pair<size_t, Vertex>( kmer_coverage[vertex_map[*it]], *it ) );


    nodes  = NodeArray(vset.size(), NULL);
    covers = std::vector<size_t>(vset.size(), 0);

    int i;
    std::multimap<size_t, Vertex>::reverse_iterator it;
    for ( it = depths.rbegin(), i=0; it != depths.rend(); ++it, i++ ) {
        Vertex v = it->second;
        nodes[i] = v;
        covers[i] = it->first;
    }
    
    // srch_log.et_ord += mytime()-t0;
}

void Extracter::findGreedyBestNeighbor( GraphPath &gpath,//NodeSet &vset,
                                        NodeArray &vset,
                                        Vertex curr_node,
                                        //PathType &trace,
                                        std::string &path_str,
                                        //std::string &nback_str,
                                        Vertex &max_vertex,
                                        int &max_value,
                                        BoundArray &max_bounds,
                                        std::string &max_trace,
                                        UsedBoundInfoList &max_useds,
                                        int direction)
{
    if ( vset.size() == 0 ) return;

    if ( Param::verbose ) {
        direction == LEFT ? 
            std::cout << "\nPrev k-mers:" : std::cout << "\nNext k-mers:";
        std::cout << vset.size() << "\n";
    }

    double tic = mytime();
    NodeArray nodes;
    std::vector<size_t> covers;
    orderNodes( vset, nodes, covers );
    srch_log.et_node_order += (mytime()-tic);

    
    if ( direction == RIGHT && path_str.size() >= (size_t)Param::back_trace ) {
        std::string curr_str = path_str.substr( path_str.size()-Param::back_trace );
        if ( Param::verbose > 1 ) std::cout << "Curr str:" << curr_str << "\n";
        //BoundArray curr_bounds(Param::nparts, BoundType());
        for ( size_t p = 0; p < (size_t)Param::nparts; p++ ) {
            double tic = mytime();
            BoundType srch = gsa[p].search((SfaChar*)curr_str.c_str(), curr_str.size());
            srch_log.et_array += (mytime()-tic);
            srch_log.ct_array++;

            //curr_bounds[p] = srch;
            if ( Param::verbose > 1 ) 
            std::cout << "Bound: (" << srch.first << ", " << srch.second << ")\n";

            if ( srch.second < srch.first ) continue;

            GsaBound b(p, srch);
            rhist_map.insert( std::pair<std::string, GsaBound>(curr_str, b) );
        }
        //rhist_map.insert( std::pair<std::string, BoundArray>(curr_str, curr_bounds) );
    }

    int k, n = (int)nodes.size();
    std::vector<int> Sums(n,0);
    std::vector<int> Supports(n,0);
    std::vector<std::vector<BoundType> > Bounds( n, std::vector<BoundType>(Param::nparts, BoundType()) );
    
    std::vector<std::string> trace_strs(n, "");
    std::vector<PathIdSet>    Check_Pids(n, check_pid);
    std::vector<PathIdPosMap> Check_Poss(n, check_pos);
    std::vector<UsedBoundInfoList> Check_Used(n, UsedBoundInfoList());
    //std::vector<UsedBoundMap> Check_Ends(n, UsedBoundMap());
    std::vector<SearchHistory> Check_Hists(n, SearchHistory());
    
    int njob = n > Param::ncpus ? Param::ncpus : n;

    bool par_run = ( Param::par_search && Param::ncpus > 1 && njob > 1 ) ? true : false;

    // temporarily set to false
    //par_run = false;

#pragma omp parallel for schedule(dynamic, 1) if (par_run) private(k) num_threads(njob)
    for ( k = 0; k < n; k++ ) {
        Vertex v = nodes[k];
        KmerId kid = vertex_map[v];

        if ( Param::verbose ) { 
            int tid = omp_get_thread_num(); 
#pragma omp critical
            std::cout << "\tt:" << tid << "\t" 
                      << alpha::IntegerToAminoAcid(kid, Param::kmer_size) 
                      << "\t" << "i:" << ancestors[v].size()
                      << " o:"  << graph::successors(graph,v).size() 
                      << " c:"  << kmer_coverage[kid] << "\n";
        }
        
        if ( (int)kmer_coverage[kid] < max_value ) {
            if ( Param::verbose ) {
#pragma omp critical
                std::cout << "\tcoverage < max -> skip\n";
            }
            continue;
        }
        
        double t0 = mytime();
        KmerId nkid = vertex_map[v];
        char ch = direction == LEFT ? 
            alpha::getFirstAminoAcid(nkid, Param::kmer_size):
            alpha::getLastAminoAcid(nkid, Param::kmer_size);
        srch_log.et_node_trace_add += (mytime()-t0);


        BoundArray ba;
        if ( direction == RIGHT && gpath.bounds.size() ) ba = gpath.bounds.back();
        assert( path_str.size() >= (size_t)Param::kmer_size );
        if ( path_str.size() == (size_t)Param::kmer_size ) {
            Supports[k] = getWeight(curr_node, v, path_str, trace_strs[k], ch, ba, Bounds[k], direction);
        } else if ( path_str.size() < (size_t)Param::back_trace ) {
            Supports[k] = getSubstringMatches( path_str, trace_strs[k], ch, ba, Bounds[k], Check_Pids[k], Check_Poss[k], direction);
        } else {
            if ( direction == RIGHT ) 
                Supports[k] = getStartReadsOnPath( path_str, trace_strs[k], ch, Bounds[k], Check_Used[k], Check_Hists[k], direction);
            else
                Supports[k] = getEndReadsOnPath( path_str, trace_strs[k], ch, Bounds[k], Check_Used[k], direction);                
        }

        if ( Supports[k] > max_value ) {
            //#pragma omp critical
            max_value  = Supports[k];
        }   
    }

    
    //-------------------------------
    // Determin maximal neighbor node
    //------------------------------- 
    max_value = 0;
    max_vertex = NULL;
    
    for ( k = 0; k < n; k++ ) {
        if ( Supports[k] > max_value ) {
            max_value  = Supports[k];
            max_vertex = nodes[k];
            max_bounds = Bounds[k];
            max_trace  = trace_strs[k];

            check_pid = Check_Pids[k];
            check_pos = Check_Poss[k];

            // if ( direction == RIGHT ) 
            //     combineMap(check_starts, Check_Starts[k]);
            // //check_starts = Check_Starts[k];
            // else
            //     combineMap(check_ends, Check_Ends[k]);
            //     //check_ends = Check_Ends[k];

            max_useds = Check_Used[k];

            // Reuse previous query during right extension
            if ( direction == RIGHT ) 
                combineHist(Check_Hists[k], direction);
        } 
    }        
    if ( Param::verbose > 1 ) {
        std::cout << "Max Used list count:" << max_useds.size() << "\n";
    }
}

// void Extracter::combineMap( UsedBoundMap &master, UsedBoundMap &slave )
// {
//     for ( auto item : slave ) {
//         master.insert( item );
//     }
// }

void Extracter::combineHist( SearchHistory &new_hist, int direction ) 
{
    for ( auto item : new_hist ) {
        //direction == LEFT ?
        //lhist_map.insert( item ) :
            rhist_map.insert( item ) ;
    }
}


int Extracter::getWeight( Vertex &curr_node, 
                          Vertex &v,
                          std::string &path_str,
                          std::string &trace_str,
                          char ch,
                          BoundArray &prev_bound,
                          BoundArray &curr_bound,
                          int direction )
{
    EdgeWeight weight = direction == LEFT ? 
        graph[edge(v, curr_node, graph).first].weight :
        graph[edge(curr_node, v, graph).first].weight;

    double t0 = mytime();
    trace_str = direction == LEFT ?
        ch + path_str : path_str + ch;
    srch_log.et_node_trace += (mytime()-t0);

    int support = (int)weight;

    t0 = mytime();
    for ( int i = 0; i < Param::nparts; i++ ) {
        BoundType srch;
        if ( direction == RIGHT && prev_bound.size() ) {
            double tic = mytime();
            srch = gsa[i].refine((SfaChar*)trace_str.c_str(), trace_str.size(), prev_bound[i].first, prev_bound[i].second);
            srch_log.et_refine += (mytime()-tic);
            srch_log.ct_refine++;
        } else {
            double tic = mytime();
            srch = gsa[i].search((SfaChar*)trace_str.c_str(), trace_str.size());
            srch_log.et_array += (mytime()-tic);
            srch_log.ct_array++;
        }
        if ( Param::verbose ) {
            //#pragma omp critical
            std::cout << "Trace str:" << trace_str << "\n"
                      << "(" << srch.first << "\t" << srch.second << ")\n";
        }
        
        curr_bound[i] = srch;
        support += srch.second>=srch.first ? (srch.second-srch.first+1) : 0;
    }
    
    return support;
}

int Extracter::getSubstringMatches( std::string &path_str,
                                    std::string &trace_str,
                                    char ch,
                                    BoundArray &prev_bound,
                                    BoundArray &bound,
                                    PathIdSet &Check_Pids,
                                    PathIdPosMap &Check_Poss,
                                    int direction )
{
    int support = 0;

    double t0 = mytime();
    trace_str = direction == LEFT ?
        ch + path_str : path_str + ch;
    srch_log.et_node_trace += (mytime()-t0);

    t0 = mytime();
    for ( int i = 0; i < Param::nparts; i++ ) {
        BoundType srch;
        if ( direction == RIGHT && prev_bound.size() ) {
            double tic = mytime();
            srch = gsa[i].refine((SfaChar*)trace_str.c_str(), trace_str.size(), prev_bound[i].first, prev_bound[i].second );
            srch_log.et_refine += (mytime()-tic);
            srch_log.ct_refine++;
            
            if ( Param::verbose ) {
#pragma omp critical
                std::cout << "Trace str:" << trace_str << "\n"
                          << "Prev: (" << prev_bound[i].first << "\t" << prev_bound[i].second << ")\n"
                          << "Curr: (" << srch.first << "\t" << srch.second << ")\n";
            }
        }
        else {
            double tic = mytime();
            srch = gsa[i].search((SfaChar*)trace_str.c_str(), trace_str.size());
            srch_log.et_array += (mytime()-tic);
            srch_log.ct_array++;
            
            if ( Param::verbose ) {
#pragma omp critical
                std::cout << "Trace str:" << trace_str << "\n"
                          << "(" << srch.first << "\t" << srch.second << ")\n";
            }
        }
        bound[i] = srch;
        support += srch.second>=srch.first ? (srch.second-srch.first+1) : 0;
    }

    int prev = getIntervalCounts(trace_str, Check_Poss, Check_Pids, direction);
    if ( prev > 0 ) srch_log.ct_iback_hit++;
    else srch_log.ct_iback_mis++;

    support -= prev;
    if ( Param::verbose ) std::cout << "Trace string:" << trace_str << "\t#prev-count:" << prev << "\tsupport:" << support << "\n";
// #pragma omp critical
// {
    srch_log.ct_iback++;
    srch_log.et_iback += (mytime()-t0);
    //}
    return support;
}

int Extracter::getSizeOfHistSearch( char ch,
                                    std::string &path_str, 
                                    size_t max_size, 
                                    int direction )
{
    double t0 = mytime();
    std::string test_str;
    if ( direction == RIGHT ) {
        test_str = ( path_str.size() < max_size ) ?
            path_str + ch :
            path_str.substr(path_str.size()-(max_size-1), max_size-1) + ch;
    }
    else {
        test_str = ( path_str.size() < max_size ) ?
            ch + path_str :
            ch + path_str.substr(0,max_size-1);
    }

    if ( Param::verbose > 1 ) std::cout << "Path str:" << test_str << "\n";
    srch_log.et_node_trace += (mytime()-t0);
    int njob = test_str.size() - (Param::back_trace);
    return njob;
}

int Extracter::getStartReadsOnPath( //GraphPath &gpath,
                                    //std::string &curr_str,
                                    //BoundArray &curr_bound,
                                   std::string &path_str,
                                   std::string &trace_str,
                                   char ch,
                                   BoundArray &bound,
                                   UsedBoundInfoList &Check_Starts,
                                   SearchHistory &Check_Hists,
                                   int direction )
{
    assert( direction == RIGHT );

    //-----------------------------------------
    // Initialize as invalid SFA range.
    // This will be used when search is failed.
    //-----------------------------------------
    bound = BoundArray(Param::nparts, BoundType(0,-1)); 
    
    size_t max_size = loader->getMaxLength();

    int njob = getSizeOfHistSearch( ch, path_str, max_size, direction );

    IntArray Used( njob*Param::nparts, 0 );
    IntArray Nums( njob*Param::nparts, 0 );
    
    assert( path_str.size() >= (size_t)Param::back_trace );
    for ( size_t i = Param::back_trace, k=0; i < max_size; i++, k++ ) {
        if ( path_str.size() < i ) break;

        double tic = mytime();
        std::string hist_str = path_str.substr(path_str.size()-i);
        std::string srch_str = hist_str + ch;
        srch_log.et_node_trace += (mytime()-tic);

        //--------------------------------------------------------------
        // Set trace string size.
        // This will be used in ralign to figure out read placment range
        //--------------------------------------------------------------
        if ( srch_str.size() == (size_t)Param::back_trace+1 ) 
            trace_str = srch_str;                

        if ( Param::verbose > 1 ) {
            std::cout << "Search str:" << srch_str << "\n";
            std::cout << "History str:" << hist_str << "\n";
        }

        //------------------------------------------
        // Find previous successful SFA search bound
        //------------------------------------------
        auto range = rhist_map.equal_range(hist_str);
        if ( range.first == range.second ) continue;

        
        tic = mytime();
        BoundType srch;
        BoundArray new_bound;
        for ( size_t p = 0; p < (size_t)Param::nparts; p++ ) {
            BoundType prev;
            bool found = false;
            for ( auto iter = range.first; iter != range.second; ++iter ) {
                if ( iter->second.first == p ) {
                    prev = iter->second.second;
                    found = true;
                    break;
                }
            }

            //---------------------------------------------
            // Require that previously successful SFA range
            //---------------------------------------------
            if ( !found ) continue;

            //-------------------
            // Bounded SFA search
            //-------------------
            double tic = mytime();
            srch = gsa[p].refine((SfaChar*)srch_str.c_str(), srch_str.size(), prev.first, prev.second);
            srch_log.et_refine += (mytime()-tic);
            srch_log.ct_refine++;
            

            if ( Param::verbose > 1 ) {
                std::cout << "Prev: (" << prev.first << ", " << prev.second << ")\n";
                std::cout << "Curr: (" << srch.first << ", " << srch.second << ")\n";
            }

            if ( srch.second < srch.first ) continue;

            //----------------------------------------------------
            // Insert the refined SFA bound of given search string
            //----------------------------------------------------
            Check_Hists.insert(std::pair<std::string,GsaBound>(srch_str,GsaBound(p, srch)));

            //---------------------------------------------------------
            // Update search bound so that later read can be retrieved.
            //---------------------------------------------------------
            if ( srch_str.size() == (size_t)Param::back_trace+1 ) {
                bound[p] = srch;
            }
            
            //--------------------------------------------------------------------------
            // Get the count of starting reads between SFA range using an auxilary array
            //--------------------------------------------------------------------------
            size_t num = read_starts[p].getCount( srch.first, srch.second );
            if ( num == 0 ) continue;
            
            Nums[k] += num; // Sum of starting counts of multiple GSAs


            //---------------------------------------------
            // Get already used starting reads of the range
            //---------------------------------------------
            size_t used = getUsedCounts( &used_starts, srch, srch_str.size(), p, 1 );
            if ( used > num ) used = num;
            Used[k] += used;

            //---------------------------------------------
            // Insert starting read count of this SFA range
            //---------------------------------------------
            UsedBound info( srch.first, srch.second, p );
            //Check_Starts.push_back( UsedBoundInfo( info, num-used ) );
            //-----------------------------------------------------
            // Actual correct count will be computed in ReadAligner
            // after placing reads
            //-----------------------------------------------------
            Check_Starts.push_back( UsedBoundInfo( info, num ) );


            
            if ( Param::verbose ) {
                //#pragma omp critical
                std::cout << "Test str:" << srch_str << "\n"
                          << "(" << srch.first << "\t" << srch.second << ")\t"
                          << "Read start:" << num << "\t"
                          << "Used start:" << used << "\t"
                          << "Read support:" << num-used << "\n";
            }
        }

        //Check_Hists.insert( std::pair<std::string, BoundArray>( srch_str, new_bound ) );
        
        srch_log.et_array += (mytime()-tic);
        srch_log.ct_array++;
        
        // tic = mytime();
        // if ( direction == LEFT ) test_str = test_str.substr(0, test_str.size()-1);
        // else test_str = test_str.substr(1);
        // srch_log.et_node_trace += (mytime()-tic);
        
    }
    int rstart = 0;
    
    for ( size_t i = 0; i < Nums.size(); i++ ) {
        rstart += (Nums[i] - Used[i]);
    }
        
    if ( Param::verbose ) std::cout << "#Read start:" << rstart << "\n";

    return std::max(rstart, 0);
}

//====================================================================
// 1. Generate n-back length string by concatenating current candidate 
//    k-mer and sub-path.
// 2. Find SFA bound of the n-back string
// 3. All the subsequence SFA queries (N-back+1 to maximum length) 
//    refines the SFA query from step 2.
//====================================================================
int Extracter::getEndReadsOnPath( std::string &path_str,
                                  std::string &trace_str,
                                  char ch,
                                  std::vector<BoundType> &bound,
                                  UsedBoundInfoList &Check_Ends,
                                  //SearchHistory &Check_Hists,
                                  int direction )
{
    assert( direction == LEFT );

    //-----------------------------------------
    // Initialize as invalid SFA range.
    // This will be used when search is failed.
    //-----------------------------------------
    bound = BoundArray(Param::nparts, BoundType(0,-1));

    //---------------------------------
    // Initial n-back length SFA search
    //---------------------------------
    assert( path_str.size() >= (size_t)Param::back_trace );
    std::string curr_str = ch + path_str.substr( 0, Param::back_trace-1 );
    if ( Param::verbose > 1 ) std::cout << "Curr str:" << curr_str << "\n";
    //------------------
    // Initial SFA range
    //------------------
    BoundArray curr_bounds(Param::nparts, BoundType());
    for ( size_t p = 0; p < (size_t)Param::nparts; p++ ) {
        double tic = mytime();
        BoundType srch = gsa[p].search((SfaChar*)curr_str.c_str(), curr_str.size());
        srch_log.et_array += (mytime()-tic);
        srch_log.ct_array++;
        
        curr_bounds[p] = srch;
        if ( Param::verbose > 1 ) 
            std::cout << "Bound: (" << srch.first << ", " << srch.second << ")\n";
    }

    size_t max_size = loader->getMaxLength();    
    int njob = getSizeOfHistSearch( ch, path_str, max_size, direction );
    
    IntArray Used(njob*Param::nparts, 0);
    IntArray Nums(njob*Param::nparts, 0);

    for ( size_t i = Param::back_trace, k=0; i < max_size; i++, k++ ) {
        if ( path_str.size() < i ) break;

        double tic = mytime();
        std::string srch_str = ch + path_str.substr(0, i);
        srch_log.et_node_trace += (mytime()-tic);

        //--------------------------------------------------------------
        // Set trace string size.
        // This will be used in ralign to figure out read placment range
        //--------------------------------------------------------------
        if ( srch_str.size() == (size_t)Param::back_trace+1 ) 
            trace_str = srch_str;
        
        if ( Param::verbose > 1 ) std::cout << "\nSearch str:" << srch_str << "\n";
        
        tic = mytime();
        size_t p, srch_len = srch_str.size();
        for ( p = 0; p < (size_t)Param::nparts; p++ ) {
            //-----------------------------
            // Require that valid SFA range
            //-----------------------------
            if ( curr_bounds[p].first > curr_bounds[p].second ) continue;

            //-------------------
            // Bounded SFA search
            //-------------------
            double tic = mytime();
            BoundType srch = gsa[p].refine((SfaChar*)srch_str.c_str(), srch_len, curr_bounds[p].first, curr_bounds[p].second );
            srch_log.et_refine += (mytime()-tic);
            srch_log.ct_refine++;

            //---------------------------------------------------------
            // Update to current SFA range so that the following search
            // (one character longer string) will refine this range.
            //---------------------------------------------------------
            curr_bounds[p] = srch;        

            //--------------
            // Invalid range
            //-------------- 
            if ( srch.second < srch.first ) continue;
            
            //---------------------------------------------------------
            // Update search bound so that later read can be retrieved.
            //---------------------------------------------------------
            if ( srch_str.size() == (size_t)Param::back_trace+1 ) {
                bound[p] = srch;
            }
            
            //=====================================================
            // Now count the number of reads end in this new range
            // 1. Get the inner range of reads end '\0'.
            //=====================================================
            tic = mytime();
            BoundType end_range = gsa[p].getEndBound( (SfaChar*)srch_str.c_str(), srch_len, srch.first, srch.second );
            srch_log.et_endsrch += (mytime()-tic);
            srch_log.ct_endsrch++;

            //-------------------------------
            // Skip if no ending reads found.
            //-------------------------------
            if ( end_range.second < end_range.first ) continue;

            size_t num = end_range.second - end_range.first + 1;
            Nums[k] += num; // Sum of end counts of multiple GSAs


            //---------------------------------------------
            // Get already used starting reads of the range
            //---------------------------------------------
            size_t used = getUsedCounts( &used_ends, end_range, srch_str.size(), p, 0 );
            if ( used > num ) used = num; // Safety
            Used[k] += used;

            //----------------------------------------
            // Insert end read count of this SFA range
            //----------------------------------------
            UsedBound info( end_range.first, end_range.second, (LcpType)p);
            // Check_Ends.push_back( UsedBoundInfo( info, num-used ) );
            //-----------------------------------------------------
            // Actual correct count will be computed in ReadAligner
            // after placing reads
            //-----------------------------------------------------
            Check_Ends.push_back( UsedBoundInfo( info, num ) );
            
            if ( Param::verbose ) {
                std::cout << "Test str:" << srch_str << "\n"
                          << "(" << srch.first << "\t" << srch.second << ")\t"
                          << "End bound:" << " "
                          << "(" << end_range.first << "\t" << end_range.second << ")\t"
                          << "Read end:" << num << "\t"
                          << "Used end:" << used << "\t"
                          << "Read support:" << num-used << "\n";
            }


        }
        
        srch_log.et_array += (mytime()-tic);
        srch_log.ct_array++;
    }

    int rend = 0;

    for ( size_t i = 0; i < Nums.size(); i++ ) {
        rend = rend + Nums[i] - Used[i];
    }

    if ( Param::verbose ) std::cout << "#Read end:" << rend << "\n";
    return std::max(rend, 0);
}

size_t Extracter::getUsedCounts( UsedBoundMap *used_map,
                                 BoundType &srch,
                                 size_t srch_len,
                                 size_t gsa_id,
                                 bool start )
{
    size_t count = 0;
    double t0 = mytime();

    //UsedBound info( srch.first, srch.second, gsa_id, srch_len );
    UsedBound info( srch.first, srch.second, (LcpType)gsa_id );
    //auto ret = used_map->equal_range(srch);
    auto it = used_map->find(info);
    start ? srch_log.et_used_start_srch += (mytime()-t0): 
        srch_log.et_used_end_srch += (mytime()-t0);
    // for ( auto rt = ret.first; rt != ret.second; ++rt) {
    //     if ( rt->second.gsa != gsa_id ) continue;
    //     if ( rt->second.len != srch_len ) continue;
        
    //     t0 = mytime();

    //     size_t used = rt->second.num;
    //     count += used;

    //     if ( Param::verbose >1 ) 
    //         std::cout << "Pre-existing used count:" << used << "\n";
    //     start ? srch_log.ct_used_start_succ++ : srch_log.ct_used_end_succ++;
    // }
            
    // if ( ret.first == ret.second ) {
    //     start ? srch_log.ct_used_start_fail++ : srch_log.ct_used_end_fail++;
    // }

    if ( it == used_map->end() ) 
        start ? srch_log.ct_used_start_fail++ : srch_log.ct_used_end_fail++;
    else {
        size_t used = it->second;
        count += used;
        if ( Param::verbose >1 ) 
            std::cout << "Pre-existing used count:" << used << "\n";
        start ? srch_log.ct_used_start_succ++ : srch_log.ct_used_end_succ++;
    }
    
    if ( start )
        srch_log.et_used_start += ( mytime()-t0 ) ;
    else
        srch_log.et_used_end += ( mytime()-t0 );

    return count;
}

// BoundType Extracter::refineRightBoundary(int p, 
//                                          BoundType &srch, 
//                                          std::string &srch_str,
//                                          size_t &srch_len, 
//                                          int &Nums, 
//                                          int &Used,
//                                          BoundArray &Bounds,
//                                          std::vector<int> &Ends)
// {
//     BoundType b;
//     int num, used;

//     int s = p*int(nreads/Param::nparts);
    
//     if ( srch.second >= srch.first ) {
//         assert( srch.first < gsa[p].getSize() );
//         GsaType item1 = gsa[p].getAt(srch.first);
//         SfaType rid1 = item1.doc;
//         LcpType pos1 = item1.pos;
        
//         if ( Param::verbose > 1 ) {
//             printf("rid:%d, pos:%d, len:%zu, read:%s, suffix:%s\n", rid1, pos1, srch_len, seqs[rid1], gsa[p].getSuffix(srch.first));
//         }
//         assert( s+rid1 < nreads );
//         if( strlen( seqs[s+rid1]+pos1 ) < srch_len ) {
//             std::cout << "wrong?:\tstrlen:" << strlen( seqs[s+rid1]+pos1 ) << "\tsearch length:" << srch_len << "\n";
//     }
//         assert( strlen( seqs[s+rid1]+pos1 ) >= srch_len );
        
//         GsaType item2 = gsa[p].getAt(srch.second);
//         SfaType rid2 = item2.doc;
//         LcpType pos2 = item2.pos;
        
//         size_t rbound = 0;
//         if ( seqs[s+rid1][pos1+srch_len] == '\0' && seqs[s+rid2][pos2+srch_len] == '\0' ) 
//             rbound = srch.second;
        
//         else if ( seqs[s+rid1][pos1+srch_len] == '\0' ) 
//             rbound = getReadEndFromGsa( srch.first, srch.second, srch_len, p, s );
        
//         if ( rbound >= (size_t)srch.first ) {
//             num = rbound-srch.first+1;
//             if ( Param::verbose > 1 ) {
//                 printf("l:%d, r:%zu, count:%d\n", srch.first, rbound, num);
//             }
//             Nums += num;
            
//             double ttic = mytime();
//             b = BoundType(srch.first, rbound);
//             UsedBoundMap::iterator rt = used_ends.find(b);
//             if ( rt != used_ends.end() ) {
//                 used = rt->second;
                
//                 if ( used > num ) used = num;
                
//                 Used += used;
//                 if ( Param::verbose >1 ) 
//                     std::cout << "Pre-existing read end count:" << used << "\n";
//                 srch_log.ct_read_end_succ++;
//             }
//             else srch_log.ct_read_end_fail++;
//             srch_log.et_read_end += ( mytime()-ttic );
            
//             //Check_Ends.insert( std::pair<BoundType, int>( srch, num ) );
//             Bounds.push_back( b );
//             Ends.push_back( num );
            
//         } else {
//             if ( Param::verbose > 1 ) std::cout << "NOT valid read ends\n";
//         }
//     }
    
//     if ( Param::verbose ) {
//         //#pragma omp critical
//         std::cout << "Test str:" << srch_str << "\t"
//                   << "(" << b.first << "\t" << b.second << ")\t"
//                   << "Read end:" << num << "\t"
//                   << "Used end:" << used << "\t"
//                   << "Read support:" << num-used << "\n";
//     }
    
//     return b;
// }

// Precondition: left most boundary ends with '\0'
// Find right most boundary ends with '\0'
size_t Extracter::getReadEndFromGsa( size_t left,
                                  size_t right,
                                  size_t length,
                                  size_t sfa_no,
                                  size_t rid_no )
{
    size_t middle = ceil((left+right)/2.0);
    //if ( left == middle ) return right;
    assert( middle>=0 );
    assert( middle < (size_t)gsa[sfa_no].getSize() );
    GsaType item = gsa[sfa_no].getAt(middle);
    int rid = rid_no + (int)item.doc;
    int pos = (int)item.pos;
    if ( Param::verbose > 1 ) {
        char *suffix = gsa[sfa_no].getSuffix(middle);
        printf("l:%zu, r:%zu, m:%zu, rid:%d, pos:%d, len:%zu, read:%s, suffix:%s\n", left, right, middle, rid, pos, length, seqs[rid], suffix );
    }
    assert( rid < nreads );
    assert( strlen( seqs[rid]+pos ) >= length );

    char ch = seqs[rid][pos+length];
    //bool found = seqs[rid] == '\0' ? true: false;
    if ( Param::verbose > 1 ) std::cout << "char:" << ch << "\n";

    if ( ch > '\0' ) {
        //left = middle;
        if ( middle == right ) return left;
        return getReadEndFromGsa(left, middle, length, sfa_no, rid_no);
    } else {
        if ( middle == right ) return right;
        return getReadEndFromGsa(middle, right, length, sfa_no, rid_no);
    }
    // if ( middle >= right ) {
    //     if ( !found ) return -1;
    //     return middle;
    // } 

    // if (found ) left = middle;
    // else right = middle;
 
    

}

//====================================================================
// Get the count of previoulsy used reads of given trace string.
//====================================================================
int Extracter::getIntervalCounts( std::string &trace_str, 
                                  PathIdPosMap &poss,
                                  PathIdSet &pids,
                                  int direction )
{
    if ( trace_str.size() == (size_t)Param::kmer_size+2 )
        findTraceStringFromPaths( trace_str, poss, pids );

    // //-----------------------------------------------------------------
    // // In case that the length of the current trace string is k+2, 
    // // then, find all previous paths that contains the current k+2 mer,
    // // by looking up k+2 mer map table
    // //-----------------------------------------------------------------
    // if ( trace_str.size() == (size_t)Param::kmer_size+2 ) {
    //     WordToPathMap::iterator wt = word_map.find(trace_str);
    //     if ( wt == word_map.end() ) return 0;
    //     if ( Param::verbose ) 
    //         std::cout << "Trace str:" << trace_str 
    //                   << "\t#paths:" << wt->second.size() << "\n";
    //     pids = wt->second;
    //     //for ( PathIdSet::iterator it = pids.begin(); it != pids.end(); ++it ) { 
    //     for ( auto pid : pids ) {
    //         //PathId pid = *it;
    //         assert( path_map.find(pid) != path_map.end() );
    //         size_t pos = path_map[pid].getReference().find(trace_str);
    //         assert( pos != std::string::npos );
    //         poss.insert( std::pair<PathId, int>( pid, (int)pos) );
    //     }
    // }

    if ( Param::verbose ) std::cout << "#interval paths:" << pids.size() << "\n";
    if ( pids.size() == 0 ) return 0;

    if ( Param::verbose ) {
        std::cout << "check pids:";
        for ( auto pid : pids ) 
            std::cout << pid << ":" << (poss.find(pid)!=poss.end()) << ",";
        std::cout << "\n";
    }
    PathIdSet bad;
    int count = 0;
    for ( auto pid : pids ) {
        assert( path_map.find(pid) != path_map.end() );
        std::string str = path_map[pid].getReference();

        // if ( poss.find(pid) == poss.end() ) {
        //     std::cerr << "pid:" << pid << "\tpos not found\n";
        //     exit(-1);
        // } else {
        //     if ( Param::verbose ) std::cout << "pid:" << pid << "\tpos:" << poss[pid] << "\n";
        // }
        assert( poss.find(pid) != poss.end() );

        int spos = poss[pid];
        if ( Param::verbose ) 
            std::cout << "pos:" << spos << "\tstr-size:" << str.size() << "\n";
        
        //------------------------------------
        // Invalid range. Consumed completely.
        //------------------------------------
        if ( spos < 0 || spos+(int)trace_str.size() > (int)str.size() ) {
            bad.insert(pid); continue;
        }

        std::string sub = str.substr(spos, trace_str.size());
        if ( Param::verbose ) std::cout << "sub:" << sub << "\n";
        
        //--------------
        // Partial match
        //--------------
        if ( sub != trace_str ) {
            bad.insert(pid); continue;
        }
        
        int lpos = spos;
        int rpos = spos+trace_str.size()-1;
        if ( Param::verbose ) 
            std::cout << "lpos:" << lpos << "\trpos:" << rpos << "\n";
      

        count += path_map[pid].getReadCount(lpos,rpos);
        
        if ( Param::verbose ) std::cout << "count:" << count << "\n";
    }

    //--------------------------------------------------
    // Upath paths and positions for next node extension
    //--------------------------------------------------
    // 1. Drop the partial matches or fully consumed paths
    //for ( PathIdSet::iterator it = bad.begin(); it != bad.end(); ++it ) {
    for ( auto pid : bad ) {
        pids.erase(pid);
        poss.erase(pid);
    }    
    // 2. Update start position for next node extension
    if ( direction == LEFT ) {
        for ( auto &item : poss ) 
            item.second--;
    }

    return count;

    // for ( PathIdSet::iterator it = pids.begin(); it != pids.end(); ) {
    //     PathId pid = *it;
    //     std::string str = path_map[pid].getReference();
    //     if ( trace_str.size() > (size_t)Param::kmer_size+2 ) {
    //         if ( direction == LEFT ) poss[pid]--;
    //     }
    //     int pos = poss[pid];
    //     if ( Param::verbose ) std::cout << "pos:" << pos << "\tstr-size:" << str.size() << "\n";
    //     if ( pos < 0 || pos+(int)trace_str.size() > (int)str.size() ) 
    //         pids.erase(it++);
    //     else {
    //         int spos = pos;
    //         if ( Param::verbose ) std::cout << "sub start:" << spos << "\n";

    //         std::string sub = str.substr(spos, trace_str.size());
    //         if ( Param::verbose ) std::cout << "sub:" << sub << "\n";

    //         if ( sub != trace_str ) pids.erase(it++);
    //         else {
    //             int lpos = pos;
    //             int rpos = pos+trace_str.size()-1;
    //             if ( Param::verbose ) std::cout << "lpos:" << lpos << "\trpos:" << rpos << "\n";
    //             count += path_map[pid].getReadCount(lpos,rpos);
    //             if ( Param::verbose ) std::cout << "count:" << count << "\n";
    //             ++it;
    //         }
    //     }
    // }
        
    // return count;
}

void Extracter::findTraceStringFromPaths( std::string &trace_str, 
                                          PathIdPosMap &poss,
                                          PathIdSet &pids )
{
    if ( Param::verbose ) std::cout << "Locating initial position of trace str:" << trace_str << "\n";
    //-----------------------------------------------------------------
    // In case that the length of the current trace string is k+2, 
    // then, find all previous paths that contains the current k+2 mer,
    // by looking up k+2 mer map table
    //-----------------------------------------------------------------
    //if ( trace_str.size() == (size_t)Param::kmer_size+2 ) {
    assert ( trace_str.size() == (size_t)Param::kmer_size+2 );
    WordToPathMap::iterator wt = word_map.find(trace_str);
    if ( wt == word_map.end() ) return;
    if ( Param::verbose ) 
        std::cout << "Trace str:" << trace_str 
                  << "\t#paths:" << wt->second.size() << "\n";
    pids = wt->second;
    for ( PathIdSet::iterator it = pids.begin(); it != pids.end(); ++it ) {
        PathId pid = *it;
        size_t pos = path_map[pid].getReference().find(trace_str);
        assert( pos != std::string::npos );
            poss.insert( std::pair<PathId, int>( pid, (int)pos) );
    }

    // Nop, we require entire trace-string match
    // //--------------------------------------------------------------------
    // // If the trace string length is greater than k+2, then the new paths
    // // from k+2 mer from the end of the path
    // //--------------------------------------------------------------------
    // else {
    //     std::string end_str = ( direction == LEFT ) ? 
    //         trace_str.substr(0,Param::kmer_size+2) :
    //         trace_str.substr(trace_str.size()-(Param::kmer_size+2), Param::kmer_size+2) ;
    //     WordToPathMap::iterator wt = word_map.find(end_str);
    //     if ( wt == word_map.end() ) return;

    //     if ( Param::verbose ) 
    //         std::cout << "Trace str:" << trace_str 
    //                   << "\t#paths:" << wt->second.size() << "\n";
    //     for ( auto pid : wt->second ) {
    //         size_t pos = path_map[pid].getReference().find(end_str);
    //         assert( pos != std::string::npos );
            
    //         //-------------------------------------
    //         // add only when the entry is not found
    //         //-------------------------------------
    //         if ( poss.find( pid ) == poss.end() ) 
    //             poss.insert( std::pair<PathId, int>( pid, (int)pos) );
    //     }
    // }
}

// int Extracter::getIntervalCounts(std::string &trace_str, int direction )
// {
//     if ( trace_str.size() == (size_t)Param::kmer_size+2 ) {
//         WordToPathMap::iterator wt = word_map.find(trace_str);
//         if ( wt == word_map.end() ) return 0;
//         if ( Param::verbose ) std::cout << "Trace str:" << trace_str << "\t#paths:" << wt->second.size() << "\n";
//         //pids = word_map.find(trace_str);
//         pids = wt->second;
//         for ( PathIdSet::iterator it = pids.begin(); it != pids.end(); ++it ) {
//             PathId pid = *it;
//             size_t pos = path_map[pid].getReference().find(trace_str);
//             assert( pos != std::string::npos );
// #pragma omp critical
//             check_pos.insert( std::pair<PathId, int>( pid, (int)pos) );
//         }
//     }

//     if ( pids.size() == 0 ) return 0;


//     int count = 0;
//     for ( PathIdSet::iterator it = pids.begin(); it != pids.end(); ) {
//         PathId pid = *it;
//         std::string str = path_map[pid].getReference();
//         if ( trace_str.size() > (size_t)Param::kmer_size+2 ) {
//             // if ( prev_dir == direction ) {
//             //     direction == LEFT ? check_pos[pid]-- : check_pos[pid]++;        
//             // }
//             // prev_dir = direction;
//             if ( direction == LEFT ) {
// #pragma omp atomic 
//                 check_pos[pid]--;
//             }
//         }
//         int pos = check_pos[pid];
//         if ( Param::verbose ) std::cout << "pos:" << pos << "\tstr-size:" << str.size() << "\n";
//         //if ( pos < 0 || pos >= (int)str.size() ) pids.erase(it++);
//         if ( pos < 0 || pos+(int)trace_str.size() > (int)str.size() ) {
// #pragma omp critical
//             pids.erase(it++);
//         }
//         else {
//             //int spos = (direction==LEFT)?  pos : pos-trace_str.size()+1;
//             int spos = pos;
//             if ( Param::verbose ) std::cout << "sub start:" << spos << "\n";
                
//             // std::string sub = ( direction == LEFT ) ?
//             //     str.substr(pos, trace_str.size() ) :
//             //     str.substr(pos-trace_str.size()+1, trace_str.size() );
//             //assert( spos+trace_str.size() <= str.size() );
//             std::string sub = str.substr(spos, trace_str.size());
//             if ( Param::verbose ) std::cout << "sub:" << sub << "\n";

//             if ( sub != trace_str ) {
// #pragma omp critical
//                 pids.erase(it++);
//             }
//             else {
//                 //int lpos = (direction == LEFT) ? pos : pos-trace_str.size()+1;
//                 //int rpos = (direction == LEFT) ? pos+trace_str.size()-1 : pos;
//                 int lpos = pos;
//                 int rpos = pos+trace_str.size()-1;
//                 if ( Param::verbose ) std::cout << "lpos:" << lpos << "\trpos:" << rpos << "\n";
//                 count += path_map[pid].getReadCount(lpos,rpos);
//                 if ( Param::verbose ) std::cout << "count:" << count << "\n";
//                 ++it;
//             }
//         }
//     }
        
//     return count;
// }

// bool Extracter::isCyclic( GraphPath  &gpath,
//                           std::string &nback_str,
//                           BoundArray &max_bound,
//                           int direction)
// {
//     if ( nback_str.size() < (size_t)Param::back_trace ) return false;

//     int k;
//     SuffixBounds::iterator it;
//     for (  it = gpath.bounds.begin(), k=0; it != gpath.bounds.end(); ++it, ++k ) {
//         BoundArray bound = *it;
//         bool match = true;
//         for ( size_t i = 0; i < bound.size(); i++ ) {
//             if ( bound[i].first  != max_bound[i].first ||
//                  bound[i].second != max_bound[i].second ) {
//                 match = false;
//                 break;
//             }
//         }
//         if ( match ) {
//             if ( getTraceStringLength( k, gpath.spos ) != Param::back_trace+1 ) {
//                 if ( Param::verbose ) std::cout << "Short sub cycle\n";
//                 continue;
//             }
                
//             if ( Param::verbose ) {
//                 std::cout << "Cyclic path:" << biostr::getSequenceString( gpath.path, vertex_map, Param::kmer_size) << "\n";
//                 std::cout << "Back string:" << nback_str << "\n";
//             }
//             return true;
//         }
//     }
//     return false;
// }

// int Extracter::getTraceStringLength( int curr, int spos )
// {
//     //int s = 0,  
//     int l = 0;
    
//     int diff = abs(spos-curr) + Param::kmer_size;

//     // left extension
//     if ( spos >= curr ) {
//         if ( diff > Param::back_trace ) 
//             l = Param::back_trace + 1;
//         else 
//             l = diff;
//     }
//     // right extension
//     else {
//         if ( curr > Param::back_trace ) 
//             l = -1 * ( Param::back_trace + 1 );
//         else 
//             l = -1 * diff;
//     }
//     return l;
// }

void Extracter::removeSubPath( GraphPath &gpath,
                               int direction )
{
    if ( (int)gpath.path.size() <= Param::back_trace ) {
        // gpath.path.clear();
        // gpath.bounds.clear();
        // gpath.traces.clear();
        gpath.clear();
        return;
    }

    typedef std::vector<BoundArray> TBounds;
    typedef std::vector<size_t> TSizes;
    NodeArray all_nodes = NodeArray(gpath.path.begin(), gpath.path.end());
    NodeArray sub_nodes;
    TBounds all_bounds = TBounds(gpath.bounds.begin(), gpath.bounds.end());
    TBounds sub_bounds;
    TSizes  all_sizes = TSizes(gpath.traces.begin(), gpath.traces.end());
    TSizes  sub_sizes;
    if ( direction == LEFT ) {
        sub_nodes = NodeArray(all_nodes.begin()+Param::back_trace,  all_nodes.end()) ;
        sub_bounds = TBounds(all_bounds.begin()+Param::back_trace, all_bounds.end());
        sub_sizes  = TSizes (all_sizes.begin()+Param::back_trace, all_sizes.end());

        if ( Param::verbose ) std::cout << "Old seed pos:" << gpath.spos << "\t";
        gpath.spos -= Param::back_trace;
        if ( Param::verbose ) std::cout << "New seed pos:" << gpath.spos << "\n";

        if ( gpath.spos < 0 ) {
            gpath.clear();
            return;
        }
        // if ( gpath.spos < 0 ) {
        //     std::cerr << "invalid spos:" << gpath.spos << "\n";
        //     exit(-1);
        // }
        // assert(gpath.spos >= 0 );
        
    } else {
        sub_nodes = NodeArray(all_nodes.begin(), all_nodes.begin() + all_nodes.size()-Param::back_trace);
        sub_bounds = TBounds(all_bounds.begin(), all_bounds.begin() + all_bounds.size()-Param::back_trace);
        sub_sizes  = TSizes (all_sizes.begin(), all_sizes.begin() + all_sizes.size()-Param::back_trace);
    }

    gpath.path = PathType(sub_nodes.begin(), sub_nodes.end());
    gpath.bounds = SuffixBounds(sub_bounds.begin(), sub_bounds.end());
    gpath.traces = TraceSizeList(sub_sizes.begin(), sub_sizes.end() );

    if ( Param::verbose ) std::cout << "New path:" << biostr::getSequenceString( gpath.path, vertex_map, Param::kmer_size ) << "\n";
}        

// bool Extracter::handleRepeat( //PathType &nback_path,
//                               std::string &nback_str,
//                               GraphPath  &gpath,
//                               CyclicPath &cpath,
//                               Vertex &curr,
//                               Vertex &max_node,
//                               int direction)
// {
//     typedef std::vector<BoundArray> TBounds;

//     if ( graph::formCycle(gpath.path, max_node) ) {
//         if ( graph::formCycle(cpath.path, max_node) ) {
//             if ( Param::verbose ) 
//                 std::cout << "\n** Cyclic path:" << biostr::getSequenceString(cpath.path, vertex_map, Param::kmer_size) << "\n";
//             NodeArray all_nodes = NodeArray(gpath.path.begin(), gpath.path.end());
//             NodeArray sub_nodes;
//             TBounds all_bounds = TBounds(gpath.bounds.begin(), gpath.bounds.end());
//             TBounds sub_bounds;
//             if ( direction == LEFT ) {
//                 sub_nodes = NodeArray(all_nodes.begin()+cpath.path.size(),  all_nodes.end()) ;
//                 sub_bounds = TBounds(all_bounds.begin()+cpath.path.size(), all_bounds.end());

//                 gpath.spos -= cpath.path.size();
//                 assert(gpath.spos >= 0 );
//             } else {
//                 sub_nodes = NodeArray(all_nodes.begin(), all_nodes.begin() + all_nodes.size()-cpath.path.size());
//                 sub_bounds = TBounds(all_bounds.begin(), all_bounds.begin() + all_bounds.size()-cpath.path.size() );
//             }
//             gpath.path = PathType(sub_nodes.begin(), sub_nodes.end());
//             gpath.bounds = SuffixBounds(sub_bounds.begin(), sub_bounds.end());

//             double tic = mytime();
//             PathType nback_path = getSubPath( gpath.path, direction, Param::back_trace);
//             nback_str  = biostr::getSequenceString(nback_path, vertex_map, Param::kmer_size);

//             srch_log.et_node_nback += (mytime()-tic);
//             srch_log.et_node_trace += (mytime()-tic);
            
//             return false;
//         }
//         if ( cpath.path.size() == 0 ) { 
//             if ( Param::verbose ) std::cout << "\n** Cycle start\n";
//             cpath.trigger = curr;
//             cpath.start = max_node;
//             cpath.path.push_back(max_node);
//         } else if ( cpath.path.size() > 0 ) {
//             if ( Param::verbose ) std::cout << "\n** Cycle extend\n";
//             if ( direction == LEFT ) 
//                 cpath.path.push_front(max_node);
//             else
//                 cpath.path.push_back(max_node);
//         }
//     } else {
//         if ( cpath.path.size() ) {
//             // do nothing
//         }
//         if ( cpath.trigger == curr ) {
//             if ( Param::verbose ) {
//                 std::cout << "\nRedundant path\n";
//                 std::cout << "Old:" << biostr::getSequenceString(gpath.path, vertex_map, Param::kmer_size) << "\n";
//                 std::cout << "Cyc:" << biostr::getSequenceString(cpath.path, vertex_map, Param::kmer_size) << "\n";
//             }
//             NodeArray nodes = NodeArray(gpath.path.begin(), gpath.path.end());
//             TBounds  bounds = TBounds(gpath.bounds.begin(), gpath.bounds.end());
//             if ( direction == LEFT ) {
//                 nodes.erase( nodes.begin(), nodes.begin()+cpath.path.size() );
//                 bounds.erase( bounds.begin(), bounds.begin()+cpath.path.size() );

//                 gpath.spos -= cpath.path.size();
//                 assert(gpath.spos >= 0 );
//             }
//             else {
//                 nodes = NodeArray(nodes.begin(), nodes.end()-cpath.path.size() );
//                 bounds = TBounds(bounds.begin(), bounds.end()-cpath.path.size() );
//             }
//             gpath.path = PathType(nodes.begin(), nodes.end());
//             gpath.bounds = SuffixBounds(bounds.begin(), bounds.end());
//             if ( Param::verbose ) std::cout << "New:" << biostr::getSequenceString(gpath.path, vertex_map, Param::kmer_size) << "\n";

//             double tic = mytime();
//             PathType nback_path = getSubPath( gpath.path, direction, Param::back_trace);
//             nback_str  = biostr::getSequenceString(nback_path, vertex_map, Param::kmer_size);            
//             srch_log.et_node_nback += (mytime()-tic);
//             srch_log.et_node_trace += (mytime()-tic);

//         }
//         cpath.path.clear();
//         cpath.trigger = cpath.start = NULL;
//     }
//     return true;
// }

// bool Extracter::redundantPath( PathType &full, PathType &trace, int direction )
// {
//     double t0 = mytime();
//     std::string trace_str = biostr::getSequenceString(trace, vertex_map, Param::kmer_size);

//     std::tr1::unordered_map<PathId, size_t> counts;
//     for ( PathType::iterator it = trace.begin(); it != trace.end(); ++it ) {
//         VertexToPathMap::iterator vt = discovered.find(*it);
//         if ( vt == discovered.end() ) continue;
//         for ( PathIdList::iterator pt = vt->second.begin(); pt != vt->second.end(); ++pt ) {
//             if ( counts.find(*pt) == counts.end() ) 
//                 counts.insert(std::pair<PathId,size_t>(*pt, 0));
//             counts[*pt]++;
//         }
//     }
    
//     PathIdList chklist;
//     std::tr1::unordered_map<PathId, size_t>::iterator it;
//     for ( it = counts.begin(); it != counts.end(); ++it ) {
//         if ( it->second >= trace.size() ) 
//             chklist.push_back(it->first);
//     }

//     if ( Param::verbose ) std::cout << "# path to compare:" << chklist.size() << "\n";
    
//     bool succ = false;
//     for ( PathIdList::iterator pt = chklist.begin(); pt != chklist.end(); ++pt ) {
//         //std::string path_str = biostr::getSequenceString(graph_paths[*pt].path, vertex_map, Param::kmer_size);
//         // if ( graph_paths[*pt].lstop == CYCLICPATH || 
//         //      graph_paths[*pt].rstop == CYCLICPATH ) continue;
        
//         if ( graph_paths[*pt].path.size() < full.size() ) continue;

//         std::string path_str = graph_paths[*pt].toString(Param::kmer_size);
//         size_t found = path_str.find(trace_str);
//         if ( found != std::string::npos ) {
//             if ( (int)found < graph_paths[*pt].spos && direction == LEFT ) 
//                 succ = true;
//             if ( (int)found > graph_paths[*pt].spos && direction == RIGHT ) 
//                 succ = true;
//             if ( succ ) { 
//                 if ( Param::verbose ) 
//                     std::cout << "Redundant path:" << *pt 
//                               << "\tmatch at:" << found 
//                               << "\tseed pos:" << graph_paths[*pt].spos
//                               << "\tdirection:" << direction
//                               << "\tlength:" << path_str.size() 
//                               << "\tstop:" << graph_paths[*pt].lstop 
//                               << graph_paths[*pt].rstop
//                               << "\t" << path_str << "\n";
//                 break;
//             }
//         }
//     }

//     if ( Param::verbose ) std::cout << "# elapsed:" << mytime()-t0 << "\n";
//     srch_log.et_red += (mytime()-t0);
//     srch_log.ct_red++;
//     return succ;
// }

// bool Extracter::redundantPath( PathType &trace )
// {
//     double t0 = mytime();
//     std::string trace_str = biostr::getSequenceString(trace, vertex_map, Param::kmer_size);

//     std::tr1::unordered_map<PathId, size_t> counts;
//     for ( PathType::iterator it = trace.begin(); it != trace.end(); ++it ) {
//         VertexToPathMap::iterator vt = discovered.find(*it);
//         if ( vt == discovered.end() ) continue;
//         for ( PathIdList::iterator pt = vt->second.begin(); pt != vt->second.end(); ++pt ) {
//             if ( counts.find(*pt) == counts.end() ) 
//                 counts.insert(std::pair<PathId,size_t>(*pt, 0));
//             counts[*pt]++;
//         }
//     }
    
//     PathIdList chklist;
//     std::tr1::unordered_map<PathId, size_t>::iterator it;
//     for ( it = counts.begin(); it != counts.end(); ++it ) {
//         if ( it->second >= trace.size() ) 
//             chklist.push_back(it->first);
//     }

//     if ( Param::verbose ) std::cout << "# path to compare:" << chklist.size() << "\n";
    
//     bool succ = false;
//     for ( PathIdList::iterator pt = chklist.begin(); pt != chklist.end(); ++pt ) {
//         //std::string path_str = biostr::getSequenceString(graph_paths[*pt].path, vertex_map, Param::kmer_size);
//         std::string path_str = graph_paths[*pt].toString(Param::kmer_size);
//         if ( path_str.find(trace_str) != std::string::npos ) {
//             if ( Param::verbose ) 
//                 std::cout << "Redundant path:" << *pt 
//                           << "\tlength:" << path_str.size() 
//                           << "\tstop:" << graph_paths[*pt].lstop << graph_paths[*pt].rstop
//                           << "\t" << path_str << "\n";
//             succ = true; break; 
//         }
//     }

//     if ( Param::verbose ) std::cout << "# elapsed:" << mytime()-t0 << "\n";
//     srch_log.et_red += (mytime()-t0);
//     srch_log.ct_red++;
//     return succ;
// }



void Extracter::displaySummary()
{
    srch_log.ct_assem = countAssembledReads();
    srch_log.printSummary();
}

void Extracter::writePath( )
{
    if ( ! Param::output_all ) return;

    std::string seq_file = Param::out_dir + "/path.fasta";

    std::fstream out;
    fio::openFile( out, seq_file.c_str(), std::ios::out );
    
    for ( ReadAlignerMap::iterator it = path_map.begin(); it != path_map.end(); ++it ) {
        PathId pid = it->first;
        ReadAligner gpath = it->second;
        std::string seq = gpath.getReference();
        out << ">p" << pid 
            << " len:" << seq.size() 
            << " stop:" << gpath.getLStop() << "/" << gpath.getRStop()
            << " trim:" << gpath.getLTrim() << "/" << gpath.getRTrim()
            << "\n";
        //<< " seed:" << str << "\n";
        out << seq << "\n";
    }
}

void Extracter::dump( std::string filename )
{
    std::fstream out;
    fio::openFile( out, filename.c_str(), std::ios::out | std::ios::binary );
    
    size_t size = path_map.size();
    out.write((char*)&size, sizeof(size_t));
    for ( ReadAlignerMap::iterator it = path_map.begin(); it != path_map.end(); ++it ) {
        PathId    pid  = it->first;
        ReadAligner r = it->second;
        out.write((char*)&(pid), sizeof(PathId));
        r.dump(out);
    }

    out.write((char*)&nreads, sizeof(int));
    out.write((char*)preads, sizeof(PathId)*nreads);
    out.close();
}

void Extracter::load( std::string filename )
{
    std::fstream in;
    fio::openFile( in, filename.c_str(), std::ios::in | std::ios::binary );
    
    path_map.clear();

    size_t npath;
    in.read((char*)&npath, sizeof(size_t));
    if ( Param::verbose ) std::cout << "npath:" << npath << "\n";

    for ( size_t n = 0; n < npath; n++ ) {
        ReadAligner r;
        PathId pid;
        in.read((char*)&(pid), sizeof(PathId));
        r.load(in);
        path_map.insert( std::pair<PathId, ReadAligner>( pid, r ) );
    }

    in.read((char*)&nreads, sizeof(int));
    //if ( preads == NULL ) {
    preads = new PathId[nreads];
    //}
    in.read((char*)preads, sizeof(PathId)*nreads);

    printElapsed( INIT_TIME, mytime(), "Graph path loaded" );

    in.close();

    status = true;
}



bool Extracter::validLengthPath( ReadAligner &raln )
{
    std::string seq = raln.getReference();
    if ( (int)seq.size() >= Param::min_length ) 
        return true;
    
    else {
        if ( Param::strict_length ) return false;

        int  lstop = raln.getLStop();
        int  rstop = raln.getRStop();
        int  ltrim = raln.getLTrim();
        int  rtrim = raln.getRTrim();
        
        if ( lstop == PATHEND && rstop == PATHEND ) {
            if ( ltrim==0 && rtrim==0 ) return true;
        }
    }
    return false;
}

void Extracter::trimPath( GraphPath &p, ReadAligner &r )
{
    int ltrim = r.getLTrim();
    int rtrim = r.getRTrim();

    if ( Param::verbose ) std::cout << "ltrim:" << ltrim << "\trtrim:" << rtrim << "\n";

    if ( ltrim == 0 && rtrim == 0 ) return;

    // if ( ltrim ) p.trim(ltrim,  true);
    // if ( rtrim ) p.trim(rtrim, false);

    NodeArray nodes = NodeArray( p.path.begin(), p.path.end() );
    if ( rtrim ) { 
        assert(rtrim < (int)nodes.size());
        nodes = NodeArray( nodes.begin(), nodes.end()-rtrim );
    }
    if ( ltrim ) {
        assert(ltrim < (int)nodes.size());
        nodes = NodeArray( nodes.begin()+ltrim, nodes.end() );
    }

    p.path = PathType(nodes.begin(), nodes.end());

}

bool Extracter::update( GraphPath &p, ReadAligner &r )
{
    //double tup = mytime();
    PathId pid = npaths; // = path_list.size();    
    //srch_log.et_update_pid += (mytime()-tup);

    // Multiple CPUs can be used
    bool valid = updateCoverage(p,r);
    if ( !valid ) return false;
            


    // //std::thread t2(&Extracter::updateCoverage, this, std::ref(p), std::ref(r));
    // std::thread t1(&Extracter::updateGraph, this, std::ref(p), std::ref(r));
    // std::thread t2(&Extracter::updateNback, this, std::ref(r));
    // std::thread t3(&Extracter::updateReads, this, std::ref(r), std::ref(pid));
    // std::thread t4(&Extracter::savePath, this, std::ref(p), std::ref(r), std::ref(pid));
    // std::thread t5(&Extracter::updateSeed, this, std::ref(p));

    // t1.join();
    // t2.join();
    // t3.join();
    // t4.join();
    // t5.join();


    //--------------------------
    // Update graph edge weights
    //--------------------------
    // Multiple CPUs can be used
    updateGraph(p,r);


    updateUsedBounds(r, 1);

    updateUsedBounds(r, 0);

    //--------------------------------------------------
    // Release temp objects from read placement and path
    //--------------------------------------------------
    r.release();
    p.release();


    //-----------------------
    // Update recruited reads
    //-----------------------
    updateReads(r, pid);

    //----------
    // Save path
    //----------
    savePath(p,r,pid);

    //--------------------------------------
    // Update seed kmers of the current path
    //--------------------------------------
    updateSeed(p);

    
    return true;
}

void Extracter::updateSeed( GraphPath &p )
{
    double t0 = mytime();
    for ( PathType::iterator it = p.path.begin(); it != p.path.end(); ++it ) {
        Vertex v = *it;

        if ( Param::trim_flag && deleted_nodes.find(v) != deleted_nodes.end() ) continue;

        assert( ancestors.find(v) != ancestors.end() );
        assert( vertex_map.find(v) != vertex_map.end());

        size_t sd = ancestors[v].size() +  graph::successors(graph, v).size();
        assert( sd > 0 );
        //if ( sd == 0 ) continue;
        
        BinType b;
        bool exist = seed_nodes.has(v,b);
        CoverageType count = kmer_coverage[ vertex_map[v] ];
        if ( (int)count < Param::min_seed ) {
            if ( exist ) seed_nodes.erase(v);
        }
        else {
            double weight = count / exp((double)sd);
            //weight = round(weight);
            weight = round(weight*10)/10;
            if ( exist ) seed_nodes.update( v, (BinType)weight );
            else seed_nodes.insert(v, (BinType)weight);
        }
    }
    srch_log.et_update_seed += (mytime()-t0);
    if ( Param::verbose ) printf("Seed update:%.4f\n", mytime()-t0);
}

void Extracter::savePath( GraphPath &p, ReadAligner &r, PathId pid )
{
    double t0 = mytime();

    //r.release(); // release temporary spaces
    //p.release();
    path_map.insert( std::pair<PathId, ReadAligner>( pid, r ) );

    //std::string str = biostr::getSequenceString(p.path, vertex_map, Param::kmer_size);
    std::string str = r.getReference();
    if ( Param::verbose ) 
        std::cout << "Path:" << pid << "\tstop:" << p.lstop << p.rstop 
                  << "\tlength:" << str.size() << "\t" << str << "\n";

    for ( PathType::iterator it = p.path.begin(); it != p.path.end(); ++it ) 
        discovered[*it].push_back(pid);


    int nword = str.size() - (Param::kmer_size+2) + 1;
    for ( int i = 0; i < nword; i++ ) {
        std::string key = str.substr(i, Param::kmer_size+2);
        //word_map[key].insert(pid);
        WordToPathMap::iterator wt = word_map.find(key);
        if ( wt == word_map.end() ) 
            word_map.insert( std::pair<std::string, std::set<PathId> >( key, std::set<PathId>() ) );
        else {
            if ( Param::verbose ) std::cout << "Existing k+2 mer:" << key << "\tsize:" << wt->second.size() << "\n";
        }
        word_map[key].insert(pid);
    }

    npaths++;

    srch_log.et_update_path += (mytime()-t0);
    if ( Param::verbose ) printf("Path saved:%.4f\n", mytime()-t0);

    //srch_log.ct_rec++;
    //srch_log.et_rec += (mytime()-t0);

    //return pid;
}
 
// void* Extracter::savePath_t( void *targ )
// {
//     double t0 = mytime();

//     //PathId pid = graph_paths.size();
//     //setPathKmers(p);
//     //graph_paths.insert( std::pair<PathId, GraphPath>( pid, p ) );

//     struct thread_data *data = (struct thread_data *)targ;
//     ReadAligner *r = data->raln;
//     GraphPath   *p = data->path;
//     //PathId     pid = data->pid;

//     path_list.push_back(*r);

//     std::string str = biostr::getSequenceString(p->path, vertex_map, Param::kmer_size);
//     if ( Param::verbose ) 
//         std::cout << "Path:" << pid << "\tstop:" << p->lstop << p->rstop 
//                   << "\tlength:" << str.size() << "\t" << str << "\n";

//     for ( PathType::iterator it = p->path.begin(); it != p->path.end(); ++it ) 
//         discovered[*it].push_back(*pid);

//     //srch_log.ct_rec++;
//     srch_log.et_rec += (mytime()-t0);

//     //return pid;
// }

bool Extracter::updateCoverage( GraphPath &gp, ReadAligner &raln )
{
    double t0 = mytime();

    WordFreqMap word_freqs;
    raln.getWordFrequency( word_freqs, Param::kmer_size, gsa, seqs );
    

    if ( Param::path_node_check ) {
        for ( PathType::iterator it = gp.path.begin(); it != gp.path.end(); ++it ) {
            std::string sstr = biostr::getSequenceString(*it, vertex_map, Param::kmer_size);
            if ( word_freqs.find(sstr) == word_freqs.end() ) {
                if ( Param::verbose ) std::cout << "Zero covered node:" << sstr << "\n";
                return false;
            }
        }
    }
        
    // if ( Param::verbose ) {
    //     std::cout << "Word frequencies\n";
    //     for ( WordFreqMap::iterator it = word_freqs.begin(); it != word_freqs.end(); ++it )
    //         std::cout << it->first << "\t" << it->second << "\n";
    // }

    for ( WordFreqMap::iterator it = word_freqs.begin(); it != word_freqs.end(); ++it ) {
        std::string str = it->first;
        KmerId kid = alpha::AminoAcidToInteger<KmerId>(str);
        if( kmer_coverage.find(kid) == kmer_coverage.end() ) {
            if ( Param::trim_flag && trimmed_kmers.find(kid) != trimmed_kmers.end() )
                continue;
            else {
                std::cout << "kmer-str:" << str << "\tkmer-id:" << kid << "\tcoverage mapping missing\n";
            exit(EXIT_FAILURE);
            }
        }
        
        assert( kmer_coverage.find(kid) != kmer_coverage.end() );
        if ( Param::verbose > 1 ) std::cout << "old coverage " << str << ":" << kmer_coverage[kid] << "\t";
        //kmer_coverage[kid] -= it->second;
        kmer_coverage[kid] -= it->second;
        if ( Param::verbose > 1 ) std::cout << "new coverage " << str << ":" << kmer_coverage[kid] << "\n";
        if ( kmer_coverage[kid] < 0 ) {
            if ( Param::verbose > 1 ) std::cout << "coverage warning " << "\n";
            kmer_coverage[kid] = 0;
        }
    }
        
        

    // for ( PathType::iterator it = gp.path.begin(); it != gp.path.end(); ++it ) {
    //     std::string sstr = biostr::getSequenceString(*it, vertex_map, Param::kmer_size);
    //     if ( word_freqs.find(sstr) == word_freqs.end() ) 
    //         std::cout << "[Warning] " << sstr << " - no supprt\n";
    //     //assert( word_freqs.find(sstr) != word_freqs.end() );
    //     else {
    //         assert( vertex_map.find(*it) != vertex_map.end() );
    //         KmerId kid = vertex_map[*it];
    //         kmer_coverage[kid] -= word_freqs[sstr];
    //     }
    // }
    //srch_log.et_upc += (mytime()-t0);
    srch_log.et_update_depth += (mytime()-t0);
    if ( Param::verbose )  printf("Coverage updated:%.4f\n", mytime()-t0);

    return true;
}

// void* Extracter::updateCoverage_t( void *targ ) //GraphPath *gp, ReadAligner *raln )
// {
//     double t0 = mytime();

//     struct thread_data *data = (struct thread_data *)targ;
//     ReadAligner *raln = data->raln;
//     GraphPath   *gp = data->path;
//     //PathId     pid = data->pid;

//     WordFreqMap word_freqs;
//     raln->getWordFrequency( word_freqs, Param::kmer_size, gsa, seqs );
    
//     if ( Param::verbose ) {
//         std::cout << "Word frequencies\n";
//         for ( WordFreqMap::iterator it = word_freqs.begin(); it != word_freqs.end(); ++it )
//             std::cout << it->first << "\t" << it->second << "\n";
//     }
//     for ( PathType::iterator it = gp->path.begin(); it != gp->path.end(); ++it ) {
//         std::string sstr = biostr::getSequenceString(*it, vertex_map, Param::kmer_size);
//         if ( word_freqs.find(sstr) == word_freqs.end() ) 
//             std::cout << "[Warning] " << sstr << " - no supprt\n";
//         //assert( word_freqs.find(sstr) != word_freqs.end() );
//         else {
//             assert( vertex_map.find(*it) != vertex_map.end() );
//             KmerId kid = vertex_map[*it];
//             kmer_coverage[kid] -= word_freqs[sstr];
//         }
//     }
//     srch_log.et_upc += (mytime()-t0);
// }


void Extracter::updateGraph( GraphPath &gp, ReadAligner &raln )
{
    double t0 = mytime();

    WordFreqMap word_freqs;
    raln.getWordFrequency( word_freqs, Param::kmer_size+1, gsa, seqs );

    for ( PathType::iterator it = gp.path.begin(); it != gp.path.end(); ++it ) {
        PathType::iterator jt = it; 
        jt++;
        if ( jt == gp.path.end() ) break;

        NodeArray nodes;
        nodes.push_back(*it); nodes.push_back(*jt);
        std::string str = biostr::getSequenceString(nodes, vertex_map, Param::kmer_size);
        
        //if ( Param::verbose ) std::cout << "substr:" << str << "\n";
        if ( word_freqs.find(str) == word_freqs.end() ) {
            if ( Param::verbose ) std::cout << "[Warning] no edge supporting reads found for " << str << "\n";
            continue;
        }
        //assert( word_freqs.find(str) != word_freqs.end() );
        //if ( word_freqs.find(str)  != word_freqs.end() ) {
        size_t cover = word_freqs[str];
        Edge e = edge(*it, *jt, graph).first;
        if ( Param::verbose > 1 ) std::cout << "old edge weight between " << str << ":" << graph[e].weight << "\t";
        graph[e].weight -= cover;
        if ( Param::verbose > 1 ) std::cout << "new edge weight between " << str << ":" << graph[e].weight << "\n";
        if ( graph[e].weight < 0 ) {
            if ( Param::verbose > 1 ) std::cout << "edge weight warning\n";
        }

        // std::string curr = biostr::getSequenceString(*it, vertex_map, Param::kmer_size);
        // std::string next = biostr::getSequenceString(*jt, vertex_map, Param::kmer_size);

        // //assert( word_freqs.find(curr) != word_freqs.end() );
        // //assert( word_freqs.find(next) != word_freqs.end() );

        // if ( word_freqs.find(curr) == word_freqs.end() ||
        //      word_freqs.find(next) == word_freqs.end() ) continue;

        // size_t cnum = word_freqs[curr];
        // size_t nnum = word_freqs[next];
    }

    if ( Param::trim_flag ) {
        dropWeakNodes( gp );
        dropIslands( gp );
    }
    srch_log.et_update_graph += (mytime()-t0);
    if ( Param::verbose ) printf("Graph updated:%.4f\n", mytime()-t0);
}

void Extracter::dropWeakNodes( GraphPath &gpath )
{
    for ( PathType::iterator it = gpath.path.begin(); it != gpath.path.end(); ++it ) {
        KmerId akid = vertex_map[*it];
        if ( kmer_coverage[akid] > (size_t)Param::min_share ) continue;

        // Possible same node in same path
        if ( deleted_nodes.find(*it) != deleted_nodes.end() ) {
            if ( Param::verbose > 1 ) std::cout << "Previously deleted nodes:" << alpha::IntegerToAminoAcid(akid, Param::kmer_size) << "\n";
            continue;
        }
        //assert( deleted_nodes.find(*it) == deleted_nodes.end());

        NodeArray succs = graph::successors(graph, *it);
        for ( NodeArray::iterator jt = succs.begin(); jt != succs.end(); ++jt ) {
            if ( edge( *it, *jt, graph ).second ) {
                double tic = mytime();
                remove_edge( *it, *jt, graph );
                srch_log.et_remove_edge += ( mytime()-tic );
                srch_log.ct_remove_edge++;
            }
            else {
                if ( Param::verbose > 1 ) 
                    std::cout << "Edge from " 
                              << alpha::IntegerToAminoAcid(vertex_map[*jt], Param::kmer_size) << " to "
                              << alpha::IntegerToAminoAcid(vertex_map[*it], Param::kmer_size) << " not exist\n";
            }
            graph::dropAncestor(ancestors, *jt, *it);         
        }
        NodeArray preds = ancestors[*it];
        for ( NodeArray::iterator jt = preds.begin(); jt != preds.end(); ++jt ) 
            if ( edge(*jt, *it, graph).second ) {
                double tic = mytime();
                remove_edge( *jt, *it, graph ); 
                srch_log.et_remove_edge += ( mytime()-tic );
                srch_log.ct_remove_edge++;
            }
            else {
                if ( Param::verbose ) 
                    std::cout << "Edge from " 
                              << alpha::IntegerToAminoAcid(vertex_map[*jt], Param::kmer_size) << " to "
                              << alpha::IntegerToAminoAcid(vertex_map[*it], Param::kmer_size) << " not exist\n";
            }
        ancestors.erase(ancestors.find(*it));

        double tic = mytime();
        remove_vertex(*it, graph);                    
        srch_log.et_remove_node += ( mytime()-tic );
        srch_log.ct_remove_node++;

        deleted_nodes.insert(*it);
        if ( Param::verbose > 1) std::cout << "Newly deleted nodes:" << alpha::IntegerToAminoAcid(akid, Param::kmer_size) << "\n";
        //         if ( param.verbose ) 
        //             std::cout << "\tDeleted:" 
        //                       << alpha::IntegerToAminoAcid(vertex_map[*it], param.kmer_size) 
        //                       << "\t"
        //                       << iindex.getValue(vertex_map[*it])->size//kmer2read_map[vertex_map[*it]]->size
        //                       << "\n";
    }
}

void Extracter::dropIslands( GraphPath &gpath )
{
    AdjacencyMap::iterator at;
    
    for ( PathType::iterator it = gpath.path.begin(); it != gpath.path.end(); ++it ) {
        //KmerId akid = vertex_map[*it];

        if ( deleted_nodes.find(*it) != deleted_nodes.end() ) 
            continue;


        NodeArray succs = graph::successors(graph, *it);
        NodeArray preds;
        at = ancestors.find(*it);
        if ( at != ancestors.end() )
            preds = at->second;

        if ( succs.size() == 0 && preds.size() == 0 ) {
            double tic = mytime();
            remove_vertex(*it, graph);            
            srch_log.et_remove_node += ( mytime()-tic );
            srch_log.ct_remove_node++;
            
            ancestors.erase(ancestors.find(*it));
            deleted_nodes.insert(*it);
            if ( Param::verbose ) std::cout << "Newly deleted nodes:" << alpha::IntegerToAminoAcid(vertex_map[*it], Param::kmer_size) << "\n";
        }
    }
}

// void* Extracter::updateGraph_t( void *targ ) // GraphPath *gp, ReadAligner *raln )
// {
//     double t0 = mytime();

//     struct thread_data *data = (struct thread_data *)targ;
//     ReadAligner *raln = data->raln;
//     GraphPath   *gp = data->path;
//     //PathId     pid = data->pid;

//     WordFreqMap word_freqs;
//     raln->getWordFrequency( word_freqs, Param::kmer_size+1, gsa, seqs );

//     for ( PathType::iterator it = gp->path.begin(); it != gp->path.end(); ++it ) {
//         PathType::iterator jt = it; 
//         jt++;
//         if ( jt == gp->path.end() ) break;

//         NodeArray nodes;
//         nodes.push_back(*it); nodes.push_back(*jt);
//         std::string str = biostr::getSequenceString(nodes, vertex_map, Param::kmer_size);
        
//         //if ( Param::verbose ) std::cout << "substr:" << str << "\n";
//         if ( word_freqs.find(str) == word_freqs.end() ) 
//             continue;
//         //assert( word_freqs.find(str) != word_freqs.end() );
//         //if ( word_freqs.find(str)  != word_freqs.end() ) {
//         size_t cover = word_freqs[str];
//         Edge e = edge(*it, *jt, graph).first;
//         graph[e].weight -= cover;
//         // if ( graph[e].weight < Parma::min_share ) 
//         //     remove_edge(*it,*jt,graph);

//         // std::string curr = biostr::getSequenceString(*it, vertex_map, Param::kmer_size);
//         // std::string next = biostr::getSequenceString(*jt, vertex_map, Param::kmer_size);

//         // //assert( word_freqs.find(curr) != word_freqs.end() );
//         // //assert( word_freqs.find(next) != word_freqs.end() );

//         // if ( word_freqs.find(curr) == word_freqs.end() ||
//         //      word_freqs.find(next) == word_freqs.end() ) continue;

//         // size_t cnum = word_freqs[curr];
//         // size_t nnum = word_freqs[next];
//     }
//     //srch_log.ct_edg++;
//     srch_log.et_upg += (mytime()-t0);
// }

// void Extracter::updateNback( ReadAligner &raln )
// {
//     //double t0 = mytime();
//     WordFreqMap word_freqs;
//     raln.getWordFrequency( word_freqs, Param::back_trace+1, gsa, seqs );

//     //if ( Param::verbose ) std::cout << "Updating n-back:\n";
//     WordFreqMap::iterator it, jt;
//     for ( it = word_freqs.begin(); it != word_freqs.end(); ++it ) {

//         jt = nback_freqs.find(it->first);
//         if ( jt == nback_freqs.end() ) {
//             nback_freqs.insert( std::pair<std::string, size_t>(it->first,it->second));
//             if ( Param::verbose > 1 ) std::cout << "N-back create:" << it->first << "\t" << nback_freqs[it->first] << "\n";
//         }
//         else {
//             jt->second += it->second;
//             if ( Param::verbose > 1 ) std::cout << "N-back update:" << it->first << "\t" << nback_freqs[it->first] << "\n";
//         }

//         // if ( Param::verbose ) 
//         //     std::cout << "Nback:" << it->first << "\t" << it->second << "\t" << nback_freqs[it->first] << "\n";
//         // nback_freqs[it->first] += it->second;
//     }

//     //srch_log.et_upb += (mytime()-t0);
// }

// void* Extracter::updateNback_t( void *targ )// ReadAligner *raln )
// {
//     double t0 = mytime();

//     struct thread_data *data = (struct thread_data *)targ;
//     ReadAligner *raln = data->raln;
//     //GraphPath   *p = data->path;
//     //PathId     pid = data->pid;


//     WordFreqMap word_freqs;
//     raln->getWordFrequency( word_freqs, Param::back_trace, gsa, seqs );

//     if ( Param::verbose ) std::cout << "Updating n-back:\n";
//     WordFreqMap::iterator it, jt;
//     for ( it = word_freqs.begin(); it != word_freqs.end(); ++it ) {

//         jt = nback_freqs.find(it->first);
//         if ( jt == nback_freqs.end() ) 
//             nback_freqs.insert( std::pair<std::string, size_t>(it->first,it->second));
//         else
//             jt->second += it->second;

//         if ( Param::verbose ) 
//             std::cout << "Nback:" << it->first << "\t" << it->second << "\t" << nback_freqs[it->first] << "\n";
//         // nback_freqs[it->first] += it->second;
//     }

//     srch_log.et_upb += (mytime()-t0);
// }

void Extracter::updateUsedBounds( ReadAligner &raln, bool start )
{
    double t0 = mytime();
    
    UsedBoundMap *used_cur = start ? raln.getReadBegBoundMap() : raln.getReadEndBoundMap();

    UsedBoundMap *used_all = start ? &used_starts : &used_ends;

    if ( Param::verbose > 1 ) {
        if ( start ) 
            std::cout << "New start bound size:" << used_cur->size() << "\n";
        else
            std::cout << "New end bound size:" << used_cur->size() << "\n";
    }

    for ( auto it = used_cur->begin(); it != used_cur->end(); ++it ) {
        if ( it->second == 0 ) continue;

        auto jt = used_all->find(it->first);
        if ( jt == used_all->end() ) used_all->insert(*it);
        else jt->second += it->second;
    }
    
    if ( Param::verbose ) {
        start ? 
            printf("Read starts updated:%.4f\n", mytime()-t0) :
            printf("Read ends updated:%.4f\n", mytime()-t0);
    }
    start ? srch_log.et_update_rstarts += (mytime()-t0) : srch_log.et_update_rends += (mytime()-t0) ;
}

// void Extracter::updateReadStarts(GraphPath &gp)
// {
//     double t0 = mytime();
    
//     for ( auto info_list : gp.used_begs ) {
//         for ( auto info : info_list ) {
//             auto range = used_starts.equal_range( info.first );
//             if ( range.first == range.second ) {
//                 used_starts.insert( info );
//                 continue;
//             }
            
//             for ( auto rt = range.first; rt != range.second;  ) {
//                 if ( rt->second.gsa != info.second.gsa ) continue;
//                 if ( rt->second.len != info.second.len ) continue;
                
//                 rt->second.num += info.second.num;
//             }
//         }
//     }
    
//     if ( Param::verbose ) printf("Read starts updated:%.4f\n", mytime()-t0);
//     srch_log.et_update_rstarts += (mytime()-t0);
// }

// void Extracter::updateReadEnds()
// {
//     double t0 = mytime();
    
//     for ( auto item : check_ends ) {
//         auto ret = used_ends.equal_range(item.first);
//         if ( ret.first == ret.second ) {
//             used_ends.insert( item );
//             continue;
//         }
        
//         for ( auto rt = ret.first; rt != ret.second; ) {
//             if ( rt->second.gsa != item.second.gsa ) continue;
//             if ( rt->second.len != item.second.len ) continue;
            
//             rt->second.num += item.second.num;
//         }
//     }
//     if ( Param::verbose ) printf("Read ends updated:%.4f\n", mytime()-t0);
//     srch_log.et_update_rends += (mytime()-t0);

// }


void Extracter::updateReads( ReadAligner &raln, PathId pid )
{
    double t0 = mytime();
    ReadPlacementList *places = raln.getPlacements();
    for ( ReadPlacementList::iterator it = places->begin(); it != places->end(); ++it ) {
        ReadId rid = it->rid;
        assert( preads[rid] == NOT_PATH );
        assert( preads[rid] != BAD_READ );
        preads[rid] = pid;
    }
    if ( Param::verbose ) printf("Reads updated:%.4f\n", mytime()-t0);
    srch_log.et_update_reads += (mytime()-t0);
}

// void* Extracter::updateReads_t( void *targ ) //ReadAligner *raln, PathId *pid )
// {
//     double t0 = mytime();

//     struct thread_data *data = (struct thread_data *)targ;
//     ReadAligner *raln = data->raln;
//     //GraphPath   *p = data->path;
//     PathId     pid = data->pid;

//     ReadPlacementList *places = raln->getPlacements();
//     for ( ReadPlacementList::iterator it = places->begin(); it != places->end(); ++it ) {
//         ReadId rid = it->rid;
//         assert( preads[rid] == NOT_PATH );
//         preads[rid] = pid;
//     }
//     srch_log.et_upr += (mytime()-t0);
// }

// void Extracter::convert()
// {
//     path_array = ReadAlignerArray( path_list.begin(), path_list.end() );
//     path_list.clear();
// }

void Extracter::release()
{
    //nback_freqs.clear();
    word_map.clear();
    check_pid.clear();
    check_pos.clear();
    //check_starts.clear();
    //check_ends.clear();
    seed_nodes.clear();
    //rhist_map.clear();
    //lhist_map.clear();

    traces_map.clear();    
    srch_failed.clear();
    read_failed.clear();
    save_failed.clear();
    used_seeds.clear();

    trimmed_nodes.clear();
    trimmed_kmers.clear();
    deleted_nodes.clear();

    purgeGraph();

    //------------------------------------------------------------------------
    // Now we can safely release substing match counts of k+2-mer to nback-mer
    //------------------------------------------------------------------------
    for ( ReadAlignerMap::iterator it = path_map.begin(); it != path_map.end(); ++it )
        it->second.dropCount();
}

void Extracter::purgeGraph()
{
    graph.clear();
    ancestors.clear();
    vertex_map.clear();
    kmer_coverage.clear();
    discovered.clear();
    //priority_map.clear();
}

int Extracter::countAssembledReads()
{
    int count = 0;
    for ( int i = 0; i < nreads; i++ )
        if ( preads[i] != NOT_PATH ) count++;
    return count;
}


void Extracter::addFlags( PathType &path, NodeFlags &flags )
{
    for ( PathType::iterator it = path.begin(); it != path.end(); ++it )
        flags[*it] = true;
}

// void Extracter::setAssembledReads( PathIdArray &assembled )
// {
//     assembled = PathIdArray( preads, preads + nreads );
// }

// void Extracter::setPathEntries( PathEntryMap &paths )
// {
//     paths.clear();
//     for ( auto entry : path_map ) {
//         PathId pid = entry.first;
//         std::string seq = entry.second.getReference();
//         int lstop = entry.second.getLStop();
//         int rstop = entry.second.getRStop();
//         int ltrim = entry.second.getLTrim();
//         int rtrim = entry.second.getRTrim();
        
//         PathEntry e( seq, lstop, rstop, ltrim, rtrim );
//         paths.insert( std::pair<PathId, PathEntry>( pid, e ) );
//     }
// }

void Extracter::makePathEntries( PathEntryMap &paths )
{
    for ( auto entry : path_map ) {
        PathId pid = entry.first;
        std::string seq = entry.second.getReference();
        int lstop = entry.second.getLStop();
        int rstop = entry.second.getRStop();
        int ltrim = entry.second.getLTrim();
        int rtrim = entry.second.getRTrim();
        
        PathEntry e( seq, lstop, rstop, ltrim, rtrim );
        paths.insert( std::pair<PathId, PathEntry>( pid, e ) );
    }
}
