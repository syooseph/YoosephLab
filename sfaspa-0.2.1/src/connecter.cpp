#include "connecter.h"

//====================================================================
// Default constructor
//====================================================================
Connecter::Connecter()
{
    status = false;
}

//====================================================================
// Destructor
//====================================================================
Connecter::~Connecter()
{
    clear();
    if ( Param::verbose ) std::cout << "Connecter cleaned\n";
}

//====================================================================
// What requires:
// Loader: reads and suffix arrays
// Assembled paths so far
// Assembled reads so far - read bridging extension only
// Extension type
//====================================================================
void Connecter::init(Loader *l, PathEntryMap *e, PathId *r, int t)
{
    //-------------------
    // Set extension type
    //-------------------
    extension_type = t;

    if ( Param::verbose >=1 ) std::cout << "Extension type:" << extension_type << "\n";

    //-----------------------
    // Set pointers to loader
    //-----------------------
    loader        = l;
    seqs          = loader->getReads();
    nreads        = loader->getCount();
    pairs         = loader->getPairs();
    //strands       = loader->getStrands();
    gsa           = loader->getGSA();
    for ( int i = 0; i < Param::nparts; i++ ) {
        int s = i*int(nreads/Param::nparts);
        gsa[i].setSequences(seqs+s);
    }

    //-----------------------
    // Assembled paths so far
    //-----------------------
    ipaths = e;

    //-----------------------
    // Assembled reads so far
    //-----------------------
    preads = r;

    //---------------------------------------
    // Assign new path IDs to assembled paths
    //---------------------------------------
    initMembership();

    //----------------------------------
    // Is it long overlapping extension?
    //----------------------------------
    log.setLJoin( extension_type == LONG_OVERLAP_EXTEND );

    //-------------------------------
    // Set pair-end read support flag
    //-------------------------------
    if ( extension_type == LONG_OVERLAP_EXTEND ) 
        log.setPESupport( Param::ljoin_pe_support );
    else 
        log.setPESupport( Param::sjoin_pe_support );

    //------------------------
    // Estimate library length
    //------------------------
    setLibraryRange();

    //-----------------------------
    // left/right extension success
    //-----------------------------
    lcount = rcount = 0;

    //-------------------------
    // Initital success status
    //-------------------------
    status = false;
}

//====================================================================
// There are no dynamically allocated objects.
//====================================================================
void Connecter::clear()
{
    ipaths = NULL;
    preads = NULL;
    status = false;
}

//====================================================================
// Set up connecter object by assigning new Path ID to existing path
//====================================================================
void Connecter::initMembership()
{
    double t0 = mytime();

    //---------------------------------------------------
    // Get the maximum numerical ID of the assembled path
    //---------------------------------------------------
    PathId pivot = getMaxPathId();
    if ( Param::verbose ) std::cout << "Max pid:" << pivot << "\n";
    if ( Param::verbose ) std::cout << "Size: " << ipaths->size() << "\n";
    
    PathEntryMap::const_iterator it;
    for ( it = ipaths->begin(); it != ipaths->end(); ++it ) {
        //-----------------------
        // Check numeric ID limit
        //-----------------------
        assert( pivot < std::numeric_limits<PathId>::max() );
        ++pivot;
        
        //-------------------------------------------------
        // Create an connect path object from a path entry
        //-------------------------------------------------
        LatchId pid = it->first;          // existing Path ID
        std::string str = it->second.seq; // path sequence
        Extender e( str );                // initiate new object
        e.members.insert(pivot);          // set itself as a member
        e.positions[pivot] = 0;           // set alignment position of itself.
        e.lstop = it->second.lstop;       // left stop condition
        e.rstop = it->second.rstop;       // right stop condition
        e.ltrim = it->second.ltrim;       // # of trimmed bases from left end
        e.rtrim = it->second.rtrim;       // # of trimmed bases from right end
        
        //-----------------------------
        // Save the just created object
        //-----------------------------
        extenders.insert( std::pair<LatchId, Extender>(pivot, e) );

        //--------------------------------------------
        // Path ID mapping between of old/new path ID
        //--------------------------------------------
        id_map.insert( std::pair<PathId, PathId>( pivot, pid ) );
    }
    
    if ( Param::verbose >=1 ) 
        std::cout << "Extension entries initialized:" << mytime()-t0 << " sec\n";
}

//====================================================================
// Get the maximum numerical path ID from entries
//====================================================================
PathId Connecter::getMaxPathId()
{
    PathId max = NOT_PATH;
    PathEntryMap::const_iterator it;
    for ( it = ipaths->begin(); it != ipaths->end(); ++it ) {
        PathId pid = it->first;
        assert( pid >= 0 );
        if ( pid > max ) max = pid;
    }
    return max;
}

//====================================================================
// Estimate range of sequence library in amino acid space
//====================================================================
void Connecter::setLibraryRange()
{
    //----------------------------
    // Maximum peptide read length
    //----------------------------
    int maxl = loader->getMaxLength();
    
    //---------------------------------
    // Library size in amino acid space
    //---------------------------------
    int size = Param::insert_size/3 + 2*maxl;
    int sd   = Param::insert_sd/3;
   
    //---------------------
    // 2 standard deviation
    //---------------------
    max_libsize = size + 2*sd;
    min_libsize = size - 2*sd;
    if ( min_libsize < 0 ) min_libsize = 0;
}

//====================================================================
// Runner
// Do path extension
//====================================================================
void Connecter::run()
{
    if ( Param::verbose ) { 
        std::cout << "\nConnecting overlapping paths ...\n";
        std::cout << "#Total paths:" << extenders.size() << "\n";
    }
    
    //----------------------------------------------------------
    // Virtual method. Right one be called from the child object
    //----------------------------------------------------------
    initMaps();

    //--------------------------------
    // Do all possible path extensions
    //--------------------------------
    joinPaths();

    //------------------
    // Drop merged paths
    //------------------
    dropMergedPaths();

    //--------------------------
    // Write assembled sequences
    //--------------------------
    writeSequences();

    //----------------------------------------
    // Release any temporarily created objects
    //----------------------------------------
    purgeTemp();

    status = true;
}

//====================================================================
// Join all latchable paths
// Latching is proceeded in the order of sequences lengths 
// from the longest to shortest paths.
//====================================================================
void Connecter::joinPaths()
{
    //---------------------------
    // Initialize progress status
    //---------------------------
    progress = Progress( 1.0, 0, extenders.size(), mytime() );

    //--------------------------------------------
    // Extend path from the longest to the shortet
    //--------------------------------------------
    PathLengthMap plen_map = getPathLengths();
    for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        //---------------------------------------------------------
        // Print progress status at every percent of path proceeded
        //---------------------------------------------------------
        double prev = progress.ratio;
        progress.count++;
        progress.showProgress();
        if ( ( Param::verbose >= 1 ) && progress.ratio > prev ) 
            log.printSummary();


        PathId sbjct_pid  = it->second;
        if ( Param::debug_id != -1 && Param::debug_id != sbjct_pid ) continue; 

        //----------------------
        // Check sequence length
        //----------------------
        std::string sbjct = extenders[sbjct_pid].sequence;
        if ( (int)sbjct.size() < Param::extend_length ) continue;
        
        //-----------------------------------------------
        // Check whether already merged to the other path
        //-----------------------------------------------
        if ( merged_paths.find( sbjct_pid ) != merged_paths.end() ) continue;

        if ( Param::verbose ) {
            std::cout << "\nIteration:" << progress.count << "\n";
            std::cout << "Pivot:" << sbjct_pid << "\tlength:" << it->first << "\n";
            std::cout << sbjct << "\n";
        }

        //-----------------------------------------
        // Reset left/right extension succes counts
        //-----------------------------------------
        lcount=rcount=0;

        //-----------
        // Latch left
        //-----------
        double tic = mytime();
        double lt0 = mytime();
        latchPath(sbjct_pid, LEFT);
        log.et_left += (mytime()-lt0);

        //------------
        // Latch right
        //------------
        lt0 = mytime();
        latchPath(sbjct_pid, RIGHT);
        log.et_right += (mytime()-lt0);

        log.et_total += (mytime()-tic);

        //================================================
        // Do not delete k-mers or end strings from maps.
        // Shorter paths are allowed to join longer paths.
        //================================================
        /*
        sbjct = extenders[sbjct_pid].sequence;

        if ( !short_flag ) {
            dropPathIds( path_map, sbjct_pid, sbjct );
        } else {
            deleteEndString(sbjct, sbjct_pid);
        }
        */
        if ( Param::verbose ) std::cout << "elapsed:" << mytime()-tic << "\n";
    }

    std::cout << "# Extended paths:" << merged_paths.size() << "\n";
}

//====================================================================
// Make a path sequence length table
// K: length, V: Path ID
//====================================================================
PathLengthMap Connecter::getPathLengths()
{
    PathLengthMap plen_map;
    for ( ExtenderMap::iterator it = extenders.begin(); it != extenders.end(); ++it ) {
        PathId pid = it->first;
        size_t len = extenders[pid].sequence.size();
        plen_map.insert( std::pair<size_t, PathId>( len, pid ) );
    }
    return plen_map;
}

//====================================================================
// Find and latch path, recursively.
//==================================================================== 
void Connecter::latchPath( PathId &spid, int direction )
{
    std::string pivot = extenders[spid].sequence;
    assert( pivot.size() > 0 );

    if ( Param::verbose ) std::cout << "Pivot Pid:" << spid << "\n"
                                    << "Sequence:" << pivot << "\n"
                                    << "Direction:"  << direction << "\n";
    //---------------------------------------------------------
    // Check existance of stop codon in case of right extension
    //---------------------------------------------------------
    if ( direction == RIGHT && pivot[pivot.size()-1] == alpha::STOP_CODON ) {
        if ( Param::verbose ) std::cout << "Path ends\n";
        if ( !Param::ignore_stop ) return;
        if ( extension_type == SHORT_OVERLAP_EXTEND ) return;
    }

    //----------------
    // Find and extend
    //----------------
    bool success = latchOverlapPaths( spid, direction );
    
    //------------------
    // Latch recursively
    //------------------
    if ( success ) latchPath( spid, direction );
    else return;
}

//====================================================================
// Drop paths that was merged during path extension
//====================================================================
void Connecter::dropMergedPaths()
{
    for ( ExtenderMap::iterator it = extenders.begin(); it != extenders.end(); ) 
        merged_paths.find(it->first) == merged_paths.end() ?
            ++it : extenders.erase(it++);
    std::cout << "# Paths:" << extenders.size() << "\n";
}

//====================================================================
// Write extended paths
//====================================================================
void Connecter::writeSequences()
{
    if ( ! Param::output_all ) return;

    //-------------------------------------
    // Make a file name wrt extenstion type
    //-------------------------------------
    std::string file = Param::out_dir;
    switch (extension_type) {
    case LONG_OVERLAP_EXTEND:
        file += "/ljoin.fasta" ;
        break;
    case  SHORT_OVERLAP_EXTEND:
        file += "/sjoin.fasta" ;
        break;
    case READ_BRIDGE_EXTEND:
        file += "/rjoin.fasta";
        break;
    }

    //----------
    // Open file
    //----------
    std::fstream out;
    fio::openFile( out, file.c_str(), std::ios::out );
    
    print(out);
}

//====================================================================
// Print out all paths
//====================================================================
void Connecter::print( std::ostream &out ) 
{
    int count = 0;
    PathLengthMap plen_map = getPathLengths();
    for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        PathId pid = it->second;
        if ( merged_paths.find(pid) != merged_paths.end() ) continue;

        std::string seq = extenders[pid].sequence;
        int lstop = extenders[pid].lstop;
        int rstop = extenders[pid].rstop;

        //----------------------------------------------
        // Make path header type based on extension type
        //----------------------------------------------
        char seq_type;
        switch (extension_type) {
        case LONG_OVERLAP_EXTEND:
            seq_type = 'l';
            break;
        case  SHORT_OVERLAP_EXTEND:
            seq_type = 's';
            break;
        case READ_BRIDGE_EXTEND:
            seq_type = 'b';
            break;
        }

        //------------------------------------------------
        // Make member header type based on extension type
        //------------------------------------------------
        char mem_type;
        switch (extension_type) {
        case LONG_OVERLAP_EXTEND:
            mem_type = Param::extend_first ? 'p' : 'c';
            break;
        case  SHORT_OVERLAP_EXTEND:
            mem_type = 'l';
            break;
        case READ_BRIDGE_EXTEND:
            mem_type = 's';
            break;
        }

        //---------------------
        // Write a FASTA header
        //--------------------- 
        out << ">" << seq_type << pid
            << " len:" << seq.size() 
            << " stop:" << lstop << "/" << rstop
            << " mem:" << extenders[pid].members.size() << " ";
        //-------------------
        // Append all members
        //-------------------
        for ( PathIdSet::iterator jt = extenders[pid].members.begin(); jt != extenders[pid].members.end(); ++jt ) {
            PathId mem = id_map[*jt];
            out << mem_type << mem << ";";
        }
        out << "\n";

        //---------------
        // Write sequence
        //---------------
        out << seq << "\n";

        count++;
    }
}

//====================================================================
// Check whether stop conditions of two paths conflicts each other.
// In latchable region, if both paths were terminated because of 
// terminal nodes, then it is a conflicting case.
// *: terminal node (source, sink, stop codon [right end])
// RIGHT
// pivot .................ooooo*
//                       *ooooo.............. match
// LEFT
//                       *ooooo.............. pivot
// match .................ooooo*
//====================================================================
bool Connecter::stopConflict( PathId spid, PathId mpid, int direction )
{
    int lend,rend;
    if ( direction == LEFT ) {
        lend = extenders[spid].lstop;
        rend = extenders[mpid].rstop;
    } else {
        lend = extenders[mpid].lstop;
        rend = extenders[spid].rstop;
    }
    if ( lend == PATHEND || rend == PATHEND ) return true;
    return false;
}

//====================================================================
// Connect two path.
// 1. Generate an extended path
// 2. Update index tables for subsequence searches
//====================================================================
void Connecter::connect(PathId spid,
                        PathId mpid,
                        std::string &pivot,
                        std::string &match,
                        AlignSummary &summary,
                        int direction)
{
    double t0 = mytime();

    //--------------------------
    // Generate an extended path
    //--------------------------
    std::string old_str = pivot;
    extend(spid, mpid, pivot, match, summary, direction );

    //------------------------------------
    // Update index table (virtual method)
    // Details are in child class
    //------------------------------------
    std::string new_str = extenders[spid].sequence;
    updateIndex(spid, old_str, new_str, direction);

    //-------------------------------------------
    // now mpid is merged and won't be used again
    //-------------------------------------------
    merged_paths.insert(mpid);
    
    direction == LEFT ? lcount++ : rcount++;
    if ( Param::verbose ) 
        printf("Success\tpid:%d\tmid:%d\tscore:%f\tlength:%d\tdirection:%d\tlcount:%d\trcount:%d\n", spid, mpid, summary.posrate, summary.length, direction, lcount, rcount);

    log.et_connect += (mytime()-t0);
}

//====================================================================
// Simple comparition of two paths in same lengths.
// Similarity is calcluated by comparing base by base.
//====================================================================
void Connecter::compareBases( AlignSummary &summary, 
                              std::string &pivot, 
                              std::string &match)
{
    //-------------------------
    // Should have same lengths
    //-------------------------
    assert( pivot.size() == match.size() );

    summary.length = pivot.length();

    //------------------------------------------------------
    // Count bases with positive score using BLOSUM62 matrix
    //------------------------------------------------------
    summary.positive = scoring::countPositive(pivot, match, BLOSUM62);

    //--------------------
    // Positive base ratio
    //--------------------
    summary.posrate = (double)summary.positive/summary.length;

    //-----------------------------
    // Zero leading & trailing gaps
    //-----------------------------
    summary.lgap.first = summary.lgap.second = 0;
    summary.egap.first = summary.egap.second = 0;
}

void Connecter::alignPair( AlignSummary &summary,
                           std::string &pivot,
                           std::string &match,
                           int direction )
{
    double tic = mytime();
   
    int anchor  = direction == RIGHT ? ANCHOR_RIGHT : ANCHOR_LEFT;

    int band_size = int ( Param::band_ratio * pivot.size() );
    SemiGlobalAlign aln = ! Param::banded_align ?
        SemiGlobalAlign(pivot, match, anchor, Param::gap_ext, Param::gap_open) :
        SemiGlobalAlign(pivot, match, anchor, Param::gap_ext, Param::gap_open, -1*band_size, band_size) ;

    summary = aln.getSummary();
    
    if ( Param::verbose ) {
        std::cout << "alignment time:" << mytime()-tic << "\n";
        aln.printAlignment(std::cout);  
        summary.print(std::cout);                             
    }
}

//====================================================================
// Make sure the correctness of the followings.
// 1. Prefix or suffix
// 2. Alignment start and gap positions of existing members
// 3. Alignment start and gap positions of new members
// 4. Stop/trim update
//==================================================================== 
void Connecter::extend( PathId &spid, 
                        PathId &mpid, 
                        std::string &pivot, 
                        std::string &match,
                        AlignSummary &summary, 
                        int direction )
{
    //--------------------------
    // Generate an extended path
    //--------------------------
    bool success = updateSequence(spid, mpid, pivot, match, summary, direction);
    if ( !success ) return;

    updateMembers(spid, mpid, summary, direction);
    updateStop(spid, mpid, direction);
    updateTrim(spid, mpid, direction);
}

//====================================================================
// Update a pivot sequence with the extended one
//====================================================================
bool Connecter::updateSequence( PathId &spid, 
                                PathId &mpid, 
                                std::string &pivot, 
                                std::string &match,
                                AlignSummary &summary, 
                                int direction )
{
    if ( direction == LEFT ) {
        if ( summary.lgap.second > 0 ) return false;
        
        //----------------------------------------------------
        // pivot ---------oooooxxxxxxxxx lgap.first = l
        //       |   l   |
        // match xxxxxxxxxooooo          lgap.second = 0
        //----------------------------------------------------
        int lgap = summary.lgap.first;            
        std::string prefix = match.substr(0, lgap);
        if ( Param::verbose ) std::cout << "Prefix:" << prefix << "\n";
        
        //-----------------------------------------
        // Prepend a portion of sequence from match
        //-----------------------------------------
        extenders[spid].sequence = prefix + pivot;
        if ( Param::verbose ) std::cout << "New Pivot:" << extenders[spid].sequence << "\n";
    }
    else {
        if ( summary.lgap.first > 0 ) return false;

        //----------------------------------------------------
        // pivot xxxxxxxxxxxxxxxooooo--------- egap.first = l
        //                           |   l   |
        // match                oooooxxxxxxxxx egap.second = 0
        //----------------------------------------------------
        int tgap = summary.egap.first; 
        std::string suffix = match.substr(match.size()-tgap);
        if ( Param::verbose ) std::cout << "Suffix:" << suffix << "\n";
        
        //----------------------------------------
        // Append a portion of sequence from match
        //----------------------------------------
        extenders[spid].sequence += suffix;
        if ( Param::verbose ) std::cout << "New Pivot:" << extenders[spid].sequence << "\n";
    }
    return true;
}

//====================================================================
// Update members:
// 1. Update existing members of current path only if direction is 
//    left and lgap is greater than 0.
// 2. Insert new alignment between sbjct & match.
//    [WARNING] Do not reset any existing position of match. 
//    It leads to inconsistent read placement.
// 3. Update existing members of matching path of if direction is
//    right and lgap is greater than 0.
//====================================================================
void Connecter::updateMembers(  PathId &spid, 
                                PathId &mpid, 
                                AlignSummary &summary, 
                                int direction )
{    
    int lgap = direction == LEFT ? summary.lgap.first : summary.lgap.second;
    if ( Param::verbose ) std::cout << "lgap:" << lgap << "\n";

    //----------------------------------------------
    // Update alignment position of existing members
    //----------------------------------------------
    updateExistingMembers( spid, mpid, lgap, direction );

    //=======================================================
    /* Shouldn't do this
     * extenders[spid].members.insert(mpid);
     * Nop. Position should not be reset.
     * extenders[spid].positions[mpid] = summary.lgap.second;
     * mpid will be inserted from the following codes of 
     * new member addition
     */
    //=======================================================
    //-----------------------------------------
    // Take members from match and add to pivot
    //-----------------------------------------
    addNewMembers( spid, mpid, lgap, direction );

    //-----------------------------------------
    // Include alignment position of match path
    //-----------------------------------------
    extenders[spid].aligns[mpid] = summary;

    if ( Param::verbose ) {
        std::cout << "#Members:" << extenders[spid].members.size() << "\n";
        for ( PathIdSet::iterator it = extenders[spid].members.begin();
              it != extenders[spid].members.end(); ++it ) {
            std::cout << *it << "\t" << extenders[spid].positions[*it] << "\n";
        }
    }
}

//====================================================================
// Update existing members of pivot path
//====================================================================
void Connecter::updateExistingMembers( const PathId &spid,
                                       const PathId &mpid,
                                       const int    &lgap,
                                       const int    &direction )
{
    if ( direction == LEFT && lgap > 0 ) {
        for ( PathIdSet::iterator mt = extenders[spid].members.begin();
              mt != extenders[spid].members.end(); ++mt ) {
            int opos = extenders[spid].positions[*mt];

            extenders[spid].positions[*mt] += lgap;                
            if ( Param::verbose ) 
                std::cout << "pid:" << *mt << "\topos:" << opos << "\t"
                          << "\tnpos:" << extenders[spid].positions[*mt] << "\n";
            if ( spid != *mt ) 
                extenders[spid].aligns[*mt].shift(lgap,true);
        }
    }
}

//====================================================================
// Add new members that are members of match path
//====================================================================
void Connecter::addNewMembers( const PathId &spid,
                               const PathId &mpid,
                               const int    &lgap,
                               const int    &direction )
{
    AlignSummary osum;
    int          opos;
    PathAlignMap::iterator pt;
    StartPosMap::iterator st;
    for ( auto opid : extenders[mpid].members ) {
        if ( Param::verbose ) std::cout << "opid:" << opid << "\n";
        //-------------------------------------------------
        // Include a member path (including match) to pivot
        //-------------------------------------------------
        extenders[spid].members.insert(opid);

        //---------------------------------------------
        // Add alignment summary of new member to pivot
        //---------------------------------------------
        if ( mpid != opid ) {
            pt = extenders[mpid].aligns.find(opid);
            assert( pt != extenders[mpid].aligns.end() );
            osum = pt->second;
            extenders[spid].aligns[opid] = osum;
        }

        //--------------------------------------
        // Append aligned position of new member
        //--------------------------------------
        st = extenders[mpid].positions.find(opid);
        assert( st != extenders[mpid].positions.end() );
        opos = st->second;
        extenders[spid].positions[opid] = opos;
        
        //------------------------------------------------
        // If direction is right with non-zero gaps,
        // new members are need to update aligned position
        //------------------------------------------------
        if ( direction == RIGHT && lgap > 0 ) {
            extenders[spid].positions[opid] += lgap;                
            extenders[spid].aligns[opid].shift(lgap,true);
        }
    }
}

//====================================================================
// Update stop codition of pivot path
//====================================================================
void Connecter::updateStop( PathId spid,
                            PathId mpid,
                            int direction )
{
    if ( direction == LEFT ) 
        extenders[spid].lstop = extenders[mpid].lstop;
    else 
        extenders[spid].rstop = extenders[mpid].rstop;
}

//====================================================================
// Update trim bases of pivot path
//====================================================================
void Connecter::updateTrim( PathId spid,
                            PathId mpid,
                            int direction )
{
    if ( direction == LEFT ) 
        extenders[spid].ltrim = extenders[mpid].ltrim;
    else 
        extenders[spid].rtrim = extenders[mpid].rtrim;
}

//====================================================================
// Write assembled objects to a binary file
//====================================================================
void Connecter::dump( std::string filename )
{
    double t0 = mytime();

    std::fstream out;
    fio::openFile( out, filename.c_str(), std::ios::out | std::ios::binary );

    out.write((char*)&extension_type, sizeof(int));

    size_t count = extenders.size();
    out.write( (char*)&count, sizeof(size_t) );
    for ( ExtenderMap::iterator it = extenders.begin(); it != extenders.end(); ++it ) {
        LatchId lid = it->first;
        out.write( (char*)&lid, sizeof(LatchId) );
        Extender ex = it->second;
        ex.dump(out);
    }

    count = merged_paths.size();
    out.write((char*)&count, sizeof(size_t));
    for ( PathIdSet::iterator it = merged_paths.begin(); it != merged_paths.end(); ++it )
        out.write((char*)&(*it), sizeof(PathId));

    count = id_map.size();
    out.write((char*)&count, sizeof(size_t));
    for ( PathIdMap::iterator it = id_map.begin(); it != id_map.end(); ++it ) {
        out.write((char*)&(it->first),  sizeof(PathId));
        out.write((char*)&(it->second), sizeof(PathId));
    }
    
    count = added_reads.size();
    out.write((char*)&count, sizeof(size_t));
    for ( PathReadsMap::iterator it = added_reads.begin(); it != added_reads.end();++it ) {
        out.write((char*)&(it->first), sizeof(PathId));
        ReadEntryList::iterator jt;
        count = it->second.size();
        out.write((char*)&count, sizeof(size_t));
        for ( jt = it->second.begin(); jt != it->second.end(); ++jt ) 
            jt->dump(out);
    }

    // out.write((char*)&nreads, sizeof(int));
    // out.write((char*)&preads, nreads*sizeof(PathId));

    if ( Param::verbose ) std::cout << "Overlap object dumped:" << mytime()-t0 << " sec\n";
    out.close();
}

//====================================================================
// Read assembled objects from a binary file
//====================================================================
void Connecter::load( std::string filename )
{
    std::fstream in;
    fio::openFile( in, filename.c_str(), std::ios::in | std::ios::binary );

    in.read((char*)&extension_type, sizeof(int));

    size_t count;
    in.read( (char*)&count, sizeof(size_t ) );
    for ( size_t i = 0; i < count; i++ ) {
        LatchId id;
        in.read( (char*)&id, sizeof(LatchId) );
        Extender ex;
        ex.load(in);
        extenders.insert( std::pair<LatchId, Extender>( id, ex) );

        //if ( Param::verbose ) std::cout << "pid:" << id << "\tseq:" << ex.sequence << "\n";
    }

    PathId pid;
    in.read((char*)&count, sizeof(size_t));
    for ( size_t i = 0; i < count; i++ ) {
        in.read((char*)&pid, sizeof(PathId));
        merged_paths.insert(pid);
    }

    PathId nid, oid;
    in.read((char*)&count, sizeof(size_t));
    for ( size_t i = 0; i < count; i++ ) {
        in.read((char*)&nid, sizeof(PathId));
        in.read((char*)&oid, sizeof(PathId));
        id_map.insert( std::pair<PathId, PathId>(nid, oid) );
    }

    in.read((char*)&count, sizeof(size_t));
    for ( size_t i = 0; i < count; i++ ) {
        in.read((char*)&(pid), sizeof(PathId));
        
        size_t cread;
        in.read((char*)&cread, sizeof(size_t));
        for ( size_t j = 0; j < cread; j++ ) {
            ReadEntry entry;
            entry.load(in);
            added_reads[pid].push_back(entry);
        }
    }

    in.close();
    if ( Param::verbose >= 1 ) std::cout << "Extension type:" << extension_type << "\n";
    printElapsed( INIT_TIME, mytime(), "Extended paths loaded " );

    status = true;
}


//====================================================================
// Make assembled path entries, so that the next assembly step uses 
// them during initialization stage.
//====================================================================
void Connecter::makePathEntries( PathEntryMap &paths )
{
    setbuf(stdout, NULL);

    for ( auto item : extenders ) {
        PathId pid = item.first;

        if ( merged_paths.find(pid) != merged_paths.end() ) continue;

        std::string seq = item.second.sequence;
        int lstop = item.second.lstop;
        int rstop = item.second.rstop;
        int ltrim = item.second.ltrim;
        int rtrim = item.second.rtrim;

        // if ( Param::verbose ) 
        //     std::cout << "path:" << pid 
        //               << "\ttype:" << extension_type 
        //               << "\tseq:" << seq << "\n";;        

        PathEntry e(seq, lstop, rstop, ltrim, rtrim );
        paths.insert( std::pair<PathId, PathEntry>( pid, e ) );
    }
}

//====================================================================
// Find no. of read pairs from pair-end library.
//====================================================================
size_t Connecter::findPairendReads( ReadFlagMap &pivot_reads,
                                    std::string &pivot, 
                                    std::string &match,
                                    int direction )
{
    std::string sub_match;

    //--------------------------------------------------
    // Extract match substring in sequence library range
    //--------------------------------------------------
    ReadFlagMap match_reads;
    sub_match = getStringInRange(match, direction, false );
    if ( Param::verbose ) {
        //std::cout << "pivot sub:" << sub_pivot << "\n";
        std::cout << "Sub-match in range:" << sub_match << "\tlength" << sub_match.size() << "\n";
    }

    //--------------------------------------------------
    // Find supporting reads in range using suffix array
    //--------------------------------------------------
    findReads( sub_match, match_reads );
    double tic = mytime();

    //------------------------------
    // Find the pair-end reads pairs
    //------------------------------
    size_t npairs = getPairedReads( pivot_reads, match_reads );
    if ( Param::verbose ) std::cout << "paired read search:" << mytime()-tic << "\n";

    return npairs;
}

//====================================================================
// Extract sequence in library range 
//====================================================================
std::string Connecter::getStringInRange(std::string &seq, 
                                        int direction, 
                                        bool pivot)
{
    /*
    |  5'  |     insert      |  3'  |
      maxl                     maxl
           ..........................................
           | dist (string in ragne) |
    */

    // Maximum read length
    int maxl = loader->getMaxLength();

    int dist = Param::insert_size/3 + maxl; 
    int sd   = Param::insert_sd/3;
    
    int max = dist + 2*sd;

    if ( direction == LEFT ) {
        if ( pivot ) {
            if ( max >= (int)seq.size() ) return seq;
            else return seq.substr(0, max);
        } else {
            if ( max >= (int)seq.size() ) return seq;
            else return seq.substr(seq.size()-max);
        }
    } else {
        if ( pivot ) {
            if ( max >=(int) seq.size() ) return seq;
            else return seq.substr(seq.size()-max);
        }
        else {
            if ( max >= (int)seq.size() ) return seq;
            else return seq.substr(0, max);
        }
    }
}


//====================================================================
// First, make overlapping sequence in chunk size.
// Then, find each chunk in suffix arrays.
//====================================================================
void Connecter::findReads( std::string &query,
                            ReadFlagMap &reads )
{
    typedef std::vector<BoundType>  BoundArray;
    typedef std::vector<BoundArray> BoundArray2D;
    
    if ( (int)query.size() < Param::chunk_size ) return;
    
    std::vector<std::string> strs;
    strs.reserve(query.size());

    //--------------------------------------------------
    // Make overlapping chunks of sequences from a query
    //--------------------------------------------------
    for ( size_t i = 0; i <= query.size()-Param::chunk_size; i++) {
        std::string str = query.substr(i, Param::chunk_size);
        strs.push_back(str);
    }

    double t0 = mytime();

    //------------------------------------------------------
    // For each sequence chunk, search it over suffix arrays
    //------------------------------------------------------
    size_t k, n = strs.size();
    BoundArray2D bound_matrix( n, BoundArray(Param::nparts, BoundType()) );
#pragma omp parallel for schedule(dynamic, 1) if (Param::ncpus>1) private(k) num_threads(Param::ncpus)
    for ( k = 0; k < n; k++ ) {
        std::string str = strs[k];
        for ( int i = 0; i < Param::nparts; i++ ) {
            BoundType srch = gsa[i].search((SfaChar*)str.c_str(), str.size());
            bound_matrix[k][i] = srch;
//             if ( Param::verbose ) 
// #pragma omp critical
//                 std::cout << "search:" << str << "\t" 
//                           << srch.first << "\t" << srch.second << "\n";
        }
    }

    if ( Param::verbose ) std::cout << "SFA search:" << mytime()-t0 << " sec\n";

    //----------------------------------------------------------------
    // For each read in suffix array match, update the matching length
    //----------------------------------------------------------------
    t0 = mytime();
    std::tr1::unordered_map<ReadId, size_t> match_lengths;
    for ( k = 0; k < n; k++ ) {
        for ( int i = 0; i < Param::nparts; i++ ) {
            BoundType srch = bound_matrix[k][i];
            if ( srch.first > srch.second ) continue;
            for ( int j = srch.first; j <= srch.second; j++ ) {
                GsaType g = gsa[i].getAt(j);
                if ( match_lengths.find(g.doc) == match_lengths.end() )
                    match_lengths.insert(std::pair<ReadId, size_t>(g.doc, Param::chunk_size));
                else
                    match_lengths[g.doc]++;
            }
        }
    }

    //----------------------------------------
    // Only count reads with sufficient length
    //----------------------------------------
    size_t count = 0;
    std::tr1::unordered_map<ReadId, size_t>::iterator it;
    for ( it = match_lengths.begin(); it != match_lengths.end(); ++it ) {
        std::string rstr = seqs[it->first];
        // if ( Param::verbose ) std::cout << "rid:" << it->first << "\tmat-len:" << it->second << "\tstr-len:" << rstr.size() << "\n";
        if ( (double)it->second/rstr.size() < Param::read_align_ratio )  continue;
        reads[it->first] = true;
        count++;
    }

    if ( Param::verbose ) std::cout << "#good reads:" << count << "\tread verification:" << mytime()-t0 << " sec\n";
}

//====================================================================
// Count the number of paired reads
//====================================================================
size_t Connecter::getPairedReads( ReadFlagMap & pivot_reads, 
                                  ReadFlagMap & match_reads )
{
    double t0 = mytime();
    if ( Param::verbose ) {
        std::cout << "#pivot_reads:" << pivot_reads.size() << "\n";
        std::cout << "#match_reads:" << match_reads.size() << "\n";
    }
    if ( pivot_reads.size() == 0 || match_reads.size() == 0 ) return 0;

    if ( Param::verbose ) {
        std::cout << "#reads (pivot):" << pivot_reads.size() << "\t";
        std::cout << "#reads (match):" << match_reads.size() << "\n";
    }
    size_t count = 0;
    for ( ReadFlagMap::iterator it = match_reads.begin();
          it != match_reads.end(); ++it ) {
        ReadId rid = it->first;
        ReadId pid = pairs[rid];
        if ( pivot_reads.find(pid) != pivot_reads.end() ) count++;
    }
    if ( Param::verbose ) {
        std::cout << "#paired reads:" << count << "\n";
        std::cout << "elapsed:" << mytime()-t0 << " sec\n";
    }
    return count;
}


