#include "assembler.h"

//====================================================================
// Constructer
//====================================================================
Assembler::Assembler( Param &p )
{
    //---------------
    // Create objects
    //---------------
    loader = new Loader();
    finder = new Extracter();
    merger = new Merger();
    metacc = new Merger();
    ljoins = new LongJoiner();
    sjoins = new TinyJoiner();
    rjoins = new ReadJoinerSA();
    placer = new Placer();
    radder = new Recruiter();
    writer = new Reporter();

    assem_reads = NULL;  // assembled reads

    //-----------------------------------------
    // Display summary of given program options
    //-----------------------------------------
    p.print(std::cout);
}

//====================================================================
// Destructor
//====================================================================
Assembler::~Assembler()
{
    clear();
}

//====================================================================
// Release all created objects
//====================================================================
void Assembler::clear()
{
    if ( loader != NULL ) delete loader;
    if ( finder != NULL ) delete finder;
    if ( merger != NULL ) delete merger;
    if ( metacc != NULL ) delete metacc;
    if ( ljoins != NULL ) delete ljoins;
    if ( sjoins != NULL ) delete sjoins;
    if ( rjoins != NULL ) delete rjoins;
    if ( placer != NULL ) delete placer;    
    if ( radder != NULL ) delete radder;    
    if ( writer != NULL ) delete writer;    

    assem_paths.clear();

    reset();
}

//====================================================================
// Reset pointers as NULL
//====================================================================
void Assembler::reset()
{
    loader = NULL;
    finder = NULL;
    merger = NULL;
    ljoins = NULL;
    sjoins = NULL;
    rjoins = NULL;
    placer = NULL;
    radder = NULL;
    writer = NULL;

    assem_reads = NULL;
}

//====================================================================
// Assemble sequences.
//====================================================================
void Assembler::run()
{
    //-------------------------
    // Load all necessary files
    //-------------------------
    load();

    //----------------------
    // Extract initial paths
    //----------------------
    extract();

    //-----------------------
    // Extend and merge paths
    //-----------------------
    if ( Param::extend_first ) {
        connect(); // Extend paths first.
        cluster(); // Then, cluster paths.
    } else {
        cluster(); // Cluster paths, first.
        connect(); // Then, extend paths.
    }

    // recluster();

    //---------------------
    // Place reads to paths
    //---------------------
    place();

    //-------------------------
    // Recruit unassigned reads
    //-------------------------
    recruit();

    //--------------
    // Write outputs
    //--------------
    report();
}

//====================================================================
// Load objects (sequences, suffix arrays, and graph input)
//====================================================================
void Assembler::load()
{
    std::cout << "\nLoading files ...\n";
    assert( loader != NULL );
    loader->load();

    printElapsed( INIT_TIME, mytime(), "Files loaded" ); 
    printMemoryUsage();
    std::cout << std::string(Param::line_length, '=') << "\n";
}

//====================================================================
// Extract paths from a graph
//====================================================================
void Assembler::extract()
{
    if ( !Param::path_flag ) return;
    std::cout << "\nExtracting paths ...\n";

    //------------------------------------
    // Load dependent objects if necessary
    //------------------------------------
    checkDependency( EXTRACT );

    //---------------------------
    // Inititalize path extractor
    //---------------------------
    finder->init(loader);

    //--------------------------
    // Extract all initial paths
    //--------------------------
    finder->run();
    if ( Param::dump_flag || Param::release_path ) finder->dump(getFileName(EXTRACT));

    //-------------------------------------------------------
    // Save extracted paths and reads to additional storages,
    // so that subsequent steps use them for initialization,
    // without refering to extracter object
    //-------------------------------------------------------
    assem_paths.clear();
    finder->makePathEntries( assem_paths );
    assem_reads = finder->getReadMembership();

    printElapsed( INIT_TIME, mytime(), "Path extracted" ); 
    std::cout << std::string(Param::line_length, '=') << "\n";

    //-----------------------------------------------------------------
    // Purge extracter if requested, in order to use memory usage less.
    // This can be useful with large amount of inputs.
    // But, extracter should be reloaded before read placements.
    //-----------------------------------------------------------------
    if ( Param::release_path ) {
        delete finder;
        finder = NULL;
        assem_reads = NULL;
    }
}

//====================================================================
// Do clustering of similar sequences
//====================================================================
void Assembler::cluster()
{
    if ( !Param::merge_flag ) return;

    std::cout << "\nClustering paths ...\n";

    //----------------------------------------
    // Load any dependent objects if necessary
    //----------------------------------------
    checkDependency( CLUSTER );

    //---------------------------------------
    // Initialize object with assembled paths
    //---------------------------------------
    merger->init( &assem_paths );


    //--------------
    // Do clustering
    //--------------
    merger->run();

    if ( Param::verbose >= 3 ) merger->print();
    if ( Param::dump_flag ) merger->dump( getFileName(CLUSTER) );

    
    //-----------------------------------------------------------------------
    // Update assembled sequences, so that so that subsequent steps use them.
    //-----------------------------------------------------------------------
    assem_paths.clear();
    merger->makePathEntries(assem_paths);
    

    printElapsed( INIT_TIME, mytime(), "Path clustered" ); 
    std::cout << std::string(Param::line_length, '=') << "\n";
}

// void Assembler::recluster()
// {
//     if ( ! Param::recluster_flag ) return;
//     if ( Param::extend_first ) return;

//     std::cout << "\nRe-clustering paths ...\n";

//     //----------------------------------------
//     // Load any dependent objects if necessary
//     //----------------------------------------
//     checkDependency( RECLUSTER );

//     //---------------------------------------
//     // Initialize object with assembled paths
//     //---------------------------------------
//     metacc->init( assem_paths );

//     metacc->setRecall(true);

//     //--------------
//     // Do clustering
//     //--------------
//     metacc->run();

//     if ( Param::verbose ) metacc->print();
//     if ( Param::dump_flag ) metacc->dump( getFileName(RECLUSTER) );

    
//     //-----------------------------------------------------------------------
//     // Update assembled sequences, so that so that subsequent steps use them.
//     //-----------------------------------------------------------------------
//     assem_paths.clear();
//     metacc->makePathEntries(assem_paths);

//     printElapsed( INIT_TIME, mytime(), "Path re-clustered" ); 
//     std::cout << std::string(Param::line_length, '=') << "\n";

    
// }

//====================================================================
// Latch overlapping sequences
//====================================================================
void Assembler::connect()
{
    if ( !Param::extend_flag ) return;

    std::cout << "\nExtending paths ...\n";

    //----------------------
    // Extend paths by types
    //----------------------
    long_overlap_extend();  // Long overlapping path extension
    tiny_overlap_extend();  // Short overlapping path extension
    read_straddle_extend(); // Path extension with bridging reads

    printElapsed( INIT_TIME, mytime(), "Path extended" ); 
    std::cout << std::string(Param::line_length, '=') << "\n";
}

//====================================================================
// Long overlapping path extension
//====================================================================
void Assembler::long_overlap_extend()
{
    if ( Param::short_overlap_only || Param::read_overlap_only ) return;

    //---------------------------------------------------
    // Check the dependency and load objects if necessary
    //---------------------------------------------------
    checkDependency( LONGJOIN );
    
    std::cout << "Extending long overlapping paths ...\n";
    //assert( assem_reads != NULL );

    //--------------------------------------------------
    // Initialize object with loader and assembled paths
    //--------------------------------------------------
    ljoins->init( loader, &assem_paths, NULL, LONG_OVERLAP_EXTEND );

    //-------------
    // Do extension
    //-------------
    ljoins->run();

    //-----------------------------------------------------------------------
    // Update assembled sequences, so that so that subsequent steps use them.
    //-----------------------------------------------------------------------
    assem_paths.clear();
    ljoins->makePathEntries(assem_paths);

    printElapsed( INIT_TIME, mytime(), "Long overlapping paths extended" ); 
    //if ( Param::verbose ) ljoins->print(std::cout);
    if ( Param::dump_flag ) ljoins->dump( getFileName(LONGJOIN) );
}

//====================================================================
// Short overlapping path extension
//====================================================================
void Assembler::tiny_overlap_extend()
{
    if ( Param::long_overlap_only || Param::read_overlap_only ) return;

    //---------------------------------------------------
    // Check the dependency and load objects if necessary
    //---------------------------------------------------
    checkDependency( SHORTJOIN );
    
    std::cout << "\nExtending short overlapping paths ...\n";

    //--------------------------------------------------
    // Initialize object with loader and assembled paths
    //--------------------------------------------------
    sjoins->init( loader, &assem_paths, NULL, SHORT_OVERLAP_EXTEND );

    //-------------
    // Do extension
    //-------------
    sjoins->run();

    //-----------------------------------------------------------------------
    // Update assembled sequences, so that so that subsequent steps use them.
    //-----------------------------------------------------------------------
    assem_paths.clear();
    sjoins->makePathEntries(assem_paths);

    printElapsed( INIT_TIME, mytime(), "Short overlapping paths extended" ); 
    //if ( Param::verbose ) sjoins->print(std::cout);
    if ( Param::dump_flag ) sjoins->dump( getFileName(SHORTJOIN) );
}

//====================================================================
// Path extension with reads straddling both paths
//====================================================================
void Assembler::read_straddle_extend()
{
    if ( Param::long_overlap_only || Param::short_overlap_only ) return;
    if ( ! Param::read_overlap_only && ! Param::read_bridge_extend ) return;

    //---------------------------------------------------
    // Check the dependency and load objects if necessary
    //---------------------------------------------------
    checkDependency( READJOIN );
    
    std::cout << "\nConnecting read bridging paths ...\n";

    //----------------------------------------
    // Make sure assembled reads are not NULL.
    // It requries to find bridging reads.
    //----------------------------------------
    assert( assem_reads != NULL );

    //---------------------------------------------------------
    // Initialize object with loader, assembled paths and reads
    //---------------------------------------------------------
    rjoins->init( loader, &assem_paths, assem_reads, READ_BRIDGE_EXTEND );

    //-------------
    // Do extension
    //-------------
    rjoins->run();

    //-----------------------------------------------------------------------
    // Update assembled sequences, so that so that subsequent steps use them.
    //-----------------------------------------------------------------------
    assem_paths.clear();
    rjoins->makePathEntries(assem_paths);

    printElapsed( INIT_TIME, mytime(), "Read bridging paths extended" ); 
    //if ( Param::verbose ) rjoins->print(std::cout);
    if ( Param::dump_flag ) rjoins->dump( getFileName(READJOIN) );
}



//====================================================================
// Place reads to meta-paths.
//====================================================================
void Assembler::place()
{
    if ( !Param::place_flag ) return;

    std::cout << "\nPlacing reads ...\n";
    
    //-----------------------------------------
    // Check depency and load objects if needed
    //-----------------------------------------
    checkDependency( PLACE );
    
    //--------------------------------------
    // Suffix arrays are no longer necessary
    //--------------------------------------
    loader->purgeGSA();

    //------------------------------------------------------
    // Initialize objects with all assembled objects so far.
    //------------------------------------------------------
    // ! Param::recluster_flag ?
    //     placer->init( loader, finder, merger, ljoins, sjoins, rjoins, &assem_paths ) :
    //     placer->init( loader, finder, merger, metacc, ljoins, sjoins, rjoins, &assem_paths ) ;
    placer->init( loader, finder, merger, ljoins, sjoins, rjoins, &assem_paths );
    
    //------------------
    // Do read placement
    //------------------
    placer->run();

    if ( Param::dump_flag ) placer->dump( getFileName(PLACE) );
    printElapsed( INIT_TIME, mytime(), "Reads placed" ); 
    std::cout << std::string(Param::line_length, '=') << "\n";
}

//====================================================================
// Recruite unassembled reads
//====================================================================
void Assembler::recruit()
{
    if ( !Param::recruit_flag ) return;

    std::cout << "\nRecruiting reads ...\n";
    

    //---------------------------------------------------
    // Check the dependency and load objects if necessary
    //---------------------------------------------------
    checkDependency( RECRUIT );

    //-------------------------------------------
    // Release suffix arrays if not released yet.
    //-------------------------------------------
    loader->purgeGSA();

    //---------------------------------------------------------
    // Initialize object with loader, assembled paths and reads
    //---------------------------------------------------------
    radder->init( loader, finder, merger, ljoins, sjoins, rjoins, &assem_paths );

    //---------------------------------
    // Copy objects from read placement
    //---------------------------------
    radder->copy(*placer);
    
    //-------------------------------
    // Necessary objects were copied.
    // Release read placement
    //-------------------------------
    placer->clear(); // no need further
    
    //--------------
    // Recruit reads
    //--------------
    radder->run();

    /* necessary for post-processing */
    if ( Param::dump_flag ) radder->dump( getFileName(RECRUIT) );

    printElapsed( INIT_TIME, mytime(), "Reads recruited" ); 
    std::cout << std::string(Param::line_length, '=') << "\n";
}

//====================================================================
// Generate SPA outputs
//====================================================================
void Assembler::report()
{
    if ( ! Param::report_flag ) return;

    //---------------------------------------------------
    // Check the dependency and load objects if necessary
    //---------------------------------------------------
    checkDependency( REPORT );

    //-------------------------------------------
    // Release suffix arrays if not released yet.
    //-------------------------------------------
    loader->purgeGSA();

    MetaPath *mpath = radder;
    std::cout << "\nGenerating assembly outputs ...\n";
        
    //---------------------------------------------------
    // Initialize object with loader and read recruitment
    //---------------------------------------------------
    writer->init(loader, mpath);

    //--------------
    // Write outputs
    //--------------
    writer->run();

    printElapsed( INIT_TIME, mytime(), "Outputs generated" ); 
    std::cout << std::string(Param::line_length, '=') << "\n";
}

//====================================================================
// Get the binary file name of given stage.
// This is called when assembly is done stage by stage.
//====================================================================
std::string Assembler::getFileName(int stage)
{
    std::string file = Param::out_dir + "/";
    switch (stage) {
    case EXTRACT:
        file += BinFiles[EXTRACT];
        break;
    case CLUSTER:
        file += BinFiles[CLUSTER];
        break;
    // case RECLUSTER:
    //     file += BinFiles[RECLUSTER];
    //     break;
    case LONGJOIN:
        file += BinFiles[LONGJOIN];
        break;
    case SHORTJOIN:
        file += BinFiles[SHORTJOIN];
        break;
    case READJOIN:
        file += BinFiles[READJOIN];
        break;
    case PLACE:
        file += BinFiles[PLACE];
        break;
    case RECRUIT:
        file += BinFiles[RECRUIT];
        break;
    default:
        std::cerr << "Invalid stage\n";
        exit(-1);
    }
    return file;
}

//====================================================================
// Check dependency of a given stage. 
// Load objects if not resides in memory
//====================================================================
void Assembler::checkDependency( int stage )
{
    //-----------------------------------------
    // Loader object is required for all stages
    //-----------------------------------------
    assert( loader != NULL );

    //------------------------------------------------
    // Make all loader objects are loaded successfully
    //------------------------------------------------
    assert( loader->getStatus() );

    switch (stage) {
        
    case EXTRACT :
        break;

    case LONGJOIN :
        //-------------------------------------------------------------
        // In case of extension first, load clusters if not loaded  yet
        //-------------------------------------------------------------
        if ( !Param::extend_first ) initiateMerger();
        
        //----------------------------
        // In case of clustering first
        //----------------------------
        else if ( assem_paths.size() == 0 ) {
            //-------------------------------------------------
            // If assembled paths are empty, load initial paths
            //-------------------------------------------------
            initiateExtracter(); 

            //---------------------------
            // Release paths if requested
            //---------------------------
            if ( Param::release_path ) {
                delete finder; finder = NULL;
                assem_reads = NULL;
            }
        }
        break;

    case SHORTJOIN :
        //if ( !Param::extend_first ) initiateMerger();
        //-------------------------------------------
        // Load long ovlerpapping paths if not loaded
        //-------------------------------------------
        initiateLongJoiner();
        break;
        
    case READJOIN :
        //if ( finder == NULL || ) std::cout << "null\n";
        // else std::cout << "not null\n";
        //assembled reads. 

        //if ( finder == NULL || ! finder->getStatus() ) loadAssembledReads();
        
        //--------------------------------------------------------------------
        // If assembled reads is NULL, load assembled reads from initial paths
        //--------------------------------------------------------------------
        if ( assem_reads == NULL ) loadAssembledReads();
        //if ( !Param::extend_first ) initiateMerger();
        //else { 
            initiateLongJoiner();
            
            //--------------------------------
            // It uses short overlapping paths
            //--------------------------------
            initiateTinyJoiner();
            //}

        // if ( finder != NULL ) std::cout << "finder is not null\n";
        // else std::cout << "Finder is still null \n";
        // if ( assem_reads == NULL ) 
        //     std::cout << "assembled reads is null\n";
        // else             std::cout << "assembled reads is NOT null\n";
        break;

    case CLUSTER :
        if ( Param::extend_first ) {
            //initiateLongJoiner();
            //initiateTinyJoiner();
            if ( !Param::read_bridge_extend ) initiateTinyJoiner();
            else initiateReadJoiner();
        } else if ( assem_paths.size() == 0 ) initiateExtracter();
        break;

    // case RECLUSTER :
    //     if ( !Param::read_bridge_extend ) initiateTinyJoiner();
    //     else initiateReadJoiner();
    //     break;
    case PLACE : case RECRUIT :
        if ( finder == NULL || ! finder->getStatus() ) reloadExtracter();
        if ( Param::extend_first ) {
            initiateLongJoiner();
            initiateTinyJoiner();
            if ( Param::read_bridge_extend ) initiateReadJoiner();
            initiateMerger();
        } else {
            initiateMerger();
            initiateLongJoiner();
            initiateTinyJoiner();
            if ( Param::read_bridge_extend ) initiateReadJoiner();
        }
        if ( stage == RECRUIT ) initiatePlacer();
        break;

    case REPORT :
        initiateRecruiter();
        break;

    default:
        std::cerr << "Invalid stage\n";
        exit(-1);
    }
}

//====================================================================
// Reload initial paths.
// NOTE: Do not update path entries
//====================================================================
void Assembler::reloadExtracter()
{
    if ( finder != NULL ) delete finder;
    finder = new Extracter();
    finder->load( getFileName(EXTRACT) );
}

//====================================================================
// Load initial paths.
// Save initial paths to assembled paths.
// Copy assembled reads to another storages.
//====================================================================
void Assembler::initiateExtracter()
{
    assert( finder != NULL );
    //if ( finder == NULL ) finder = new Extracter();
    if ( !finder->getStatus() ) {
        finder->load( getFileName(EXTRACT) );

        assem_paths.clear();
        finder->makePathEntries( assem_paths );
        assem_reads = finder->getReadMembership();
    }
}

//====================================================================
// Load cluster objects.
// Copy clustered paths to another storage (assem_paths)
//====================================================================
void Assembler::initiateMerger()
{
    //if ( merger == NULL ) merger = new Merger();
    assert( merger != NULL );
    if ( !merger->getStatus() ) {
        merger->load( getFileName(CLUSTER) );
        
        assem_paths.clear();
        merger->makePathEntries( assem_paths );
    }
}

//====================================================================
// Load long overlapping path extension objects.
// Copy extended paths to another storage (assem_paths)
//====================================================================
void Assembler::initiateLongJoiner()
{
    //if ( ljoins == NULL ) ljoins = LongJoiner();
    assert( ljoins != NULL );
    if ( ! ljoins->getStatus() ) {
        ljoins->load( getFileName(LONGJOIN) );

        assem_paths.clear();
        ljoins->makePathEntries( assem_paths );
    }
}

//====================================================================
// Load shorte overlapping path extension objects.
// Copy extended paths to another storage (assem_paths)
//====================================================================
void Assembler::initiateTinyJoiner()
{
    assert( sjoins != NULL );
    if ( ! sjoins->getStatus() ) {
        sjoins->load( getFileName(SHORTJOIN) );
        
        assem_paths.clear();
        sjoins->makePathEntries( assem_paths );
    }
}

//====================================================================
// Load read bridging path extension objects.
// Copy extended paths to another storage (assem_paths)
//====================================================================
void Assembler::initiateReadJoiner()
{
    assert( rjoins != NULL );
    if ( ! rjoins->getStatus() ) {
        rjoins->load( getFileName(READJOIN) );
        
        assem_paths.clear();
        rjoins->makePathEntries( assem_paths );
    }
}

//====================================================================
// Load read placement objects.
//====================================================================
void Assembler::initiatePlacer()
{
    assert( placer != NULL );
    if ( ! placer->getStatus() ) 
        placer->load( getFileName( PLACE ) );
}

//====================================================================
// Load read recruitment objects.
//====================================================================
void Assembler::initiateRecruiter()
{
    assert( radder != NULL );
    if ( ! radder->getStatus() ) 
        radder->load( getFileName( RECRUIT ) );
}

//====================================================================
// Load extracter and update assembled reads
//====================================================================
void Assembler::loadAssembledReads()
{
    if ( finder == NULL ) finder = new Extracter();
    finder->load( getFileName(EXTRACT) );
    
    assem_reads = finder->getReadMembership();
    assert( assem_reads != NULL );
}
