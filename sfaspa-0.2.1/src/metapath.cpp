#include "metapath.h"

MetaPath::MetaPath()
{
    preads = NULL;
    status = false;
    reset();
}

MetaPath::~MetaPath()
{
    clear();
}

/**
 * Copy constructor 
 */
MetaPath::MetaPath(const MetaPath &source)
{
    copy(source);
}

/**
 * Operator overloading
 */
MetaPath& MetaPath::operator= (const MetaPath &source)
{
    copy(source);
    return *this;
}

void MetaPath::copy( const MetaPath &other )
{
    Alns      = other.Alns;
    Id2Paths  = other.Id2Paths;
    Path2Ids  = other.Path2Ids;
    bad_paths = other.bad_paths;

    assert( preads != NULL );
    assert( other.preads != NULL );
    memcpy( preads, other.preads, nreads*sizeof(PathId));

    if ( Param::verbose >= 1 ) 
        std::cout << "Other metapaths copied\n";
}

void MetaPath::reset()
{
    clear();
}

void MetaPath::clear()
{
    if ( Alns.size() ) Alns.clear();
    if ( Id2Paths.size() ) Id2Paths.clear();
    if ( Path2Ids.size() ) Path2Ids.clear();
    if ( bad_paths.size() ) bad_paths.clear();
    if ( preads != NULL ) delete[] preads;
    preads = NULL;
    created   = false;
    status    = false;
    if ( Param::verbose ) std::cout << "MetaPath cleaned\n";
}

void MetaPath::init( Loader *loader,
                     Extracter *finder,
                     Merger *merger,
                     Connecter *joiner,
                     Connecter *jshort,
                     Connecter *rjoiner,
                     PathEntryMap *pentries )
                     // Summary *assems)
{
    /* Read and suffix array */
    seqs          = loader->getReads();
    nreads        = loader->getCount();

    read_aligns = finder->getReadAlignerMap();

    clusters = merger->getClusters();
    // path_entries = assems->getPathEntryMap();
    path_entries = pentries;

    tiny_exts = jshort->getExtenders();
    long_exts = joiner->getExtenders();
    
    tiny_idmap = jshort->getPathIdMap();
    long_idmap = joiner->getPathIdMap();
    cluster_idmap = merger->getPathIdMap();


    if ( Param::read_bridge_extend ) {
        read_exts = rjoiner->getExtenders();
        read_idmap = rjoiner->getPathIdMap();
        added_reads = rjoiner->getAddedReads();
    }

    size_t i;
    if ( !Param::extend_first ) {
        ExtenderMap::iterator it;
        if ( Param::read_bridge_extend ) {                              
            assert( rjoiner != NULL );
            for ( it = read_exts->begin(), i=0; it != read_exts->end(); ++it, ++i ) {
                Id2Paths.insert( std::pair<size_t, PathId>(i, (PathId)it->first) );
                Path2Ids.insert( std::pair<PathId, size_t>((PathId)it->first, i) );
            }
            npaths = read_exts->size();
        } else {
            for ( it = tiny_exts->begin(), i=0; it != tiny_exts->end(); ++it, ++i ) {
                Id2Paths.insert( std::pair<size_t, PathId>(i, (PathId)it->first) );
                Path2Ids.insert( std::pair<PathId, size_t>((PathId)it->first, i) );
            }
            npaths = tiny_exts->size();
        }
    } else {
        Clusters::iterator it;
        for ( it = clusters->begin(), i=0; it != clusters->end(); ++it, ++i ) {
            Id2Paths.insert( std::pair<size_t, PathId>(i, (PathId)it->first) );
            Path2Ids.insert( std::pair<PathId, size_t>((PathId)it->first, i) );
        }
        npaths = clusters->size();
    }

    //std::cout << "#paths:" << npaths << "\n";

    Alns = PathAligners( npaths, PathAligner() );
    bad_paths = FlagArray(npaths, false);

    preads = new PathId[nreads];

    created = true;

    if ( Param::verbose >= 1 ) 
        std::cout << "Metapaths inititalized\n";
}


void MetaPath::init( Loader *loader,
                     Extracter *finder,
                     Merger *merger,
                     Merger *metacc,
                     Connecter *joiner,
                     Connecter *jshort,
                     Connecter *rjoiner,
                     PathEntryMap *pentries )
                     // Summary *assems)
{
    /* Read and suffix array */
    seqs          = loader->getReads();
    nreads        = loader->getCount();

    read_aligns = finder->getReadAlignerMap();

    clusters = merger->getClusters();
    meta_ccs = metacc->getClusters();

    // path_entries = assems->getPathEntryMap();
    path_entries = pentries;

    tiny_exts = jshort->getExtenders();
    long_exts = joiner->getExtenders();
    
    tiny_idmap = jshort->getPathIdMap();
    long_idmap = joiner->getPathIdMap();

    if ( Param::read_bridge_extend ) {
        read_exts = rjoiner->getExtenders();
        read_idmap = rjoiner->getPathIdMap();
        added_reads = rjoiner->getAddedReads();
    }

    size_t i;
    Clusters::iterator it;
    for ( it = meta_ccs->begin(), i=0; it != meta_ccs->end(); ++it, ++i ) {
        Id2Paths.insert( std::pair<size_t, PathId>(i, (PathId)it->first) );
        Path2Ids.insert( std::pair<PathId, size_t>((PathId)it->first, i) );
    }
    npaths = clusters->size();

    //std::cout << "#paths:" << npaths << "\n";

    Alns = PathAligners( npaths, PathAligner() );
    bad_paths = FlagArray(npaths, false);

    preads = new PathId[nreads];

    created = true;

    if ( Param::verbose >= 1 ) 
        std::cout << "Metapaths inititalized\n";
}

void MetaPath::run()
{
    /* Do nothing */
}


void MetaPath::write()
{
    if ( ! Param::output_all ) return;

    writeAlignment();
    writePlacement();
    writeSequences();
}

void MetaPath::writePlacement()
{
    std::fstream *outs = new std::fstream[Param::ncpus];

    std::string place_type = ( derived == DERIVED_PLACE ) ?
        PlacerFiles[PLACE_TXT] : RecruiterFiles[PLACE_TXT];
    
    int i;
    for ( i = 0; i < Param::ncpus; i++ ) {
        std::string filename = Param::out_dir + "/"
            + place_type + "."
            //+ boost::lexical_cast<std::string>(i);
            + std::to_string(i);
        fio::openFile( outs[i], filename.c_str(), std::ios::out );
    }
    
#pragma omp parallel for schedule(dynamic, 1) if(Param::ncpus>1) private(i) num_threads(Param::ncpus)    
    for ( i = 0; i < (int)npaths; i++ ) {

        if ( Param::debug_id != -1 && Param::debug_id != Id2Paths[i] ) continue;
        if ( bad_paths[i] ) continue;
        
        int thread = omp_get_thread_num();
        
        std::string raw = Alns[i].getSequence();
        std::string con = Alns[i].getConsensus();
        outs[thread] << "ID:" << i << "\n";
        outs[thread] << "Path:" << Id2Paths[i] << "\n";
        outs[thread] << "Assembled:" << raw << "\n";
        if ( con != "" ) outs[thread] << "Consensus:" << con << "\n";
        Alns[i].printPlacement(outs[thread]);
        outs[thread] << "//\n";
    }
    
    for ( i = 0; i < Param::ncpus; i++ ) outs[i].close();
    delete[] outs;
    
    combine( place_type );
}

void MetaPath::combine( std::string name )
{
    std::string nfile = Param::out_dir + "/" + name;
    if ( Param::ncpus == 1 ) {
        std::string ofile = Param::out_dir + "/" + name + ".0";
        rename( ofile.c_str(), nfile.c_str() );
        return;
    }
    
    std::fstream out;
    fio::openFile( out, nfile.c_str(), std::ios::out );

    std::string line;
    for ( int i = 0; i < Param::ncpus; i++ ) {
        std::stringstream ofile;
        ofile << Param::out_dir 
              << "/" <<  name << "." << i;
        std::ifstream sin( ofile.str().c_str() );
        while ( getline(sin, line ) ) 
            out << line << "\n";

        if ( remove( ofile.str().c_str() ) != 0 ) 
            perror( "Error deleting file" );
    }

    out.close();
}


void MetaPath::writeAlignment()
{
    if ( ! Param::output_all ) return;
    if ( ! Param::align_flag  && ! Param::profile_flag ) return;
    
    std::fstream *aouts = new std::fstream[Param::ncpus];
    std::fstream *pouts = new std::fstream[Param::ncpus];

    std::string aln_type = ( derived == DERIVED_PLACE ) ?
        PlacerFiles[ALIGN] : RecruiterFiles[ALIGN];
    std::string pro_type = ( derived == DERIVED_PLACE ) ?
        PlacerFiles[PROFILE] : RecruiterFiles[PROFILE];

    int i;
    for ( i = 0; i < Param::ncpus; i++ ) {
        std::string filename = Param::out_dir + "/"
            + aln_type + "."
            //            + boost::lexical_cast<std::string>(i);
            + std::to_string(i);
        fio::openFile( aouts[i], filename.c_str(), std::ios::out );
        filename = Param::out_dir + "/"
            + pro_type + "."
            //+ boost::lexical_cast<std::string>(i);
            + std::to_string(i);
        fio::openFile( pouts[i], filename.c_str(), std::ios::out );
    }
    
#pragma omp parallel for schedule(dynamic, 1) if(Param::ncpus>1) private(i) num_threads(Param::ncpus)    
    for ( i = 0; i < (int)npaths; i++ ) {
        if ( bad_paths[i] ) continue;

        if ( Param::debug_id != -1 && Param::debug_id != Id2Paths[i] ) continue;

        int thread = omp_get_thread_num();

        aouts[thread] << "ID:" << i << "\n";
        aouts[thread] << "Path:" << Id2Paths[i] << "\n";
        MSA msa( &Alns[i], seqs) ;
        msa.print(aouts[thread],  Param::line_length);
        aouts[thread] << "//\n";
        
        if ( Param::verbose > 1 ) {
            std::cout << "MSA created\n";
            std::cout << "Pivot-ref:" << Alns[i].getSequence() << "\n";
            std::cout << "Pivot-con:" << Alns[i].getConsensus() << "\n";
            std::cout << "Pivot-msa:" << msa.getPivot() << "\n";
        }

        if ( !Param::profile_flag )  continue;
        Profile pro(&msa, AA_TYPE);
        std::string con = pro.makeConsensus();
        pouts[thread] << "ID:" << i << "\n";
        pouts[thread] << "Path:" << Id2Paths[i] << "\n";
        pouts[thread] << "Consensus:" << con << "\n";
        pro.print(pouts[thread]);
        Alns[i].setConsensus(con);
        pouts[thread] << "//\n";
        if ( pro.replacedStopCodon() ) {
            if ( Param::verbose > 1 ) 
                std::cout << "Id:" << i 
                          << "\tpid:" << Id2Paths[i] 
                          << "\tstop codon replaced\n";
        }
        if ( Param::verbose > 1 ) {
            std::cout << "Profile created\n";
            std::cout << "Pivot-ref:" << Alns[i].getSequence() << "\n";
            std::cout << "Pivot-con:" << Alns[i].getConsensus() << "\n";
        }

        Alns[i].setDirty(false);
    }

    for ( i = 0; i < Param::ncpus; i++ ) aouts[i].close();
    for ( i = 0; i < Param::ncpus; i++ ) pouts[i].close();
    delete[] aouts;
    delete[] pouts;
    
    combine( aln_type );
    combine( pro_type );

}

void MetaPath::dump(std::string filename)
{
    double t0 = mytime();

    std::fstream out;
    fio::openFile( out, filename.c_str(), std::ios::out | std::ios::binary );

    out.write( (char*)&nreads, sizeof(int) );


    size_t count = Alns.size();
    if ( Param::verbose ) std::cout << "Path-count:" << count << "\n";
    out.write( (char*)&count, sizeof(size_t) );
    for ( size_t i = 0; i < count; i++ ) {
        PathId pid = Id2Paths[i];
        out.write( (char*)&pid, sizeof(PathId) );
        Alns[i].dump(out);
        if ( Param::verbose > 1 ) std::cout << i << "\tpath:" << pid << "\tbad:" << bad_paths[i] << "\t" << Alns[i].getSequence() << "\n";
    }

    for ( size_t i = 0; i < count; i++ ) {
        bool flag = bad_paths[i];
        //char flag = (char)bad_paths[i];
        if ( Param::verbose > 1 )
            printf( "No:%zu\tPid:%d\tFlag:%d\n", i, Id2Paths[i], flag );
        out.write( (char*)&flag, sizeof(bool) );
        //out.write( (char*)&flag, sizeof(char) );
    }

    for ( int i = 0; i < nreads; i++ )
        out.write( (char*)&preads[i], sizeof(PathId) );

    if ( Param::verbose ) 
        std::cout << "Read alignment dumped:" << mytime()-t0 << " sec\n";
}

void MetaPath::load(std::string filename)
{
    std::fstream in;
    fio::openFile( in, filename.c_str(), std::ios::in | std::ios::binary );

    in.read( (char*)&nreads, sizeof(int) );

    // create preads
    if ( !created ) preads = new PathId[nreads];

    size_t count;
    in.read( (char*)&count, sizeof(size_t ) );
    if ( Param::verbose ) std::cout << "Path-count:" << count << "\n";
    
    // create 
    if ( !created ) Alns = PathAligners(count, PathAligner());
    for ( size_t i = 0; i < count; i++ ) {
        PathId pid;
        in.read( (char*)&pid, sizeof(PathId) );
        Id2Paths[i] = pid;
        Path2Ids[pid] = i;
        PathAligner aln;
        aln.load(in);
        Alns[i] = aln;

        if ( Param::verbose > 1 ) std::cout << i << "\tpath:" << pid << "\t" << aln.getSequence() << "\n";
    }

    bad_paths = FlagArray(count, false);
    bool flag;
    //char flag;
    for ( size_t i = 0; i < count; i++ ) {
        in.read( (char*)&flag, sizeof(bool) );
        //in.read( (char*)&flag, sizeof(char) );
        if ( Param::verbose > 1 )
            printf( "No:%zu\tPid:%d\tFlag:%d\n", i, Id2Paths[i], flag );
        bad_paths[i] = (bool)flag;
    }
    
    assert( preads != NULL );
    PathId pid;
    for ( int i = 0; i < nreads; i++ ) {
        in.read( (char*)&pid, sizeof(PathId) );
        preads[i] = pid;
    }

    in.close();
    printElapsed( INIT_TIME, mytime(), "Read placement loaded" );

    status = true;
}

void MetaPath::writeSequences()
{
    std::string seq_type = ( derived == DERIVED_PLACE ) ?
        PlacerFiles[SEQUENCE] : RecruiterFiles[SEQUENCE];

    std::string file = Param::out_dir + "/" + seq_type;
    std::fstream out;
    fio::openFile( out, file.c_str(), std::ios::out );

    //int count = 0;
    PathLengthMap plen_map = getPathLengths();
    PathLengthMap::reverse_iterator it;
    for ( it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        PathId pid = it->second;
        size_t num = Path2Ids[pid];
        if ( bad_paths[num] ) continue;
        if ( Param::debug_id != -1 && Param::debug_id != pid ) continue;

        std::string ref = Alns[num].getSequence();
        std::string con = Alns[num].getConsensus();
        std::string seq = (con=="") ? ref : con;

        // remove gaps
        for ( int i = (int)seq.size()-1; i>=0; i-- )
            if ( seq[i] == '-' ) seq.erase(i,1);
        
        int rsize = Alns[num].getSize();
        int lstop = Alns[num].getLStop();
        int rstop = Alns[num].getRStop();
        
        // bool ltrim = Alns[num].getLTrim();
        // bool rtrim = Alns[num].getRTrim();

        // char lchar = ltrim ? 't' : 'f';
        // char rchar = rtrim ? 't' : 'f';
        //out << ">spa_" << num
        
        
        derived == DERIVED_PLACE ? out << ">m" << num : out << ">r" << num;
        if ( Param::extend_first ) out << " pid:c" << pid;
        else out << " pid:s" << pid;
        out << " len:" << seq.size() 
            << " reads:" << rsize
            << " stop:" << lstop << "/" << rstop << "\n";
            // << " trim:" << ltrim << "/" << rtrim << "\n";
        //<< " trim:" << lchar << rchar << "\n";
        out << seq << "\n";
    }
}

PathLengthMap MetaPath::getPathLengths()
{
    PathLengthMap plen_map;
    for ( size_t i = 0; i < npaths; i++ ) {
        if ( bad_paths[i] ) continue;
        PathId pid = Id2Paths[i];
        size_t len = Alns[i].getSequence().size();
        plen_map.insert( std::pair<size_t, PathId>( len, pid ) );
    }    
    return plen_map;
}


void MetaPath::updatePathMembership()
{
    // reset
    memset( preads, NOT_PATH, nreads*sizeof(PathId));

    int dup = 0;
    for ( size_t i = 0; i < npaths; i++ ) {
        if ( bad_paths[i] ) continue;
        
        ReadPlacementList *places = Alns[i].getPlacements();        
        ReadPlacementList::iterator it;
        for ( it = places->begin(); it != places->end(); ++it ) {
            if ( preads[it->rid] != NOT_PATH ) {
                std::cout << "Redundant placement\tread:" << it->rid << "\told path:" << preads[it->rid] << "\tnew path:" << Id2Paths[i] << "\tsequence:" << seqs[it->rid] << "\n";
                dup++;
            }
            preads[it->rid] = Id2Paths[i];
        }
    }
    if ( dup > 0 ) std::cout << "Redundant placement:" << dup << "\n";
    
    if ( Param::verbose >= 1 ) 
        std::cout << "Path membership updated\n";
}

int MetaPath::countAssembledReads()
{
    updatePathMembership();

    int count = 0;
    for ( int i = 0; i < nreads; i++ ) 
        if ( preads[i] != NOT_PATH ) count++;
    return count;
}

