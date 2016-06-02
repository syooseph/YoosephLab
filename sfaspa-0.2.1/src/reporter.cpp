#include "reporter.h"

Reporter::Reporter()
{
}

Reporter::~Reporter()
{
    //seqs      = NULL;
    // Alns      = NULL;
    // Id2Paths  = NULL;
    // Path2Ids  = NULL;
    // bad_paths = NULL;
    if ( Param::verbose ) std::cout << "Reporter cleaned\n";
}

void Reporter::init(Loader *loader, MetaPath *mpath)
{
    seqs      = loader->getReads();

    meta_paths = mpath;
    Alns       = meta_paths->getAlignments();
    Id2Paths   = meta_paths->getPathIds();
    Path2Ids   = meta_paths->getPathNums();
    bad_paths  = meta_paths->getBadPaths();

    npaths     = Alns->size();

    if ( Param::verbose ) std::cout << "#paths:" << npaths << "\n";
}

void Reporter::run()
{
    write();
}

void Reporter::write()
{
    writeAlignment();
    writePlacement();
    writeSequences();

    //----------------------
    // binary read placement
    //---------------------- 
    std::string filename = Param::out_dir + "/" + SpaFiles[PLACE_BIN];
    meta_paths->dump(filename);
}

void Reporter::writeAlignment()
{
    //if ( ! Param::align_flag  && ! Param::profile_flag ) return;
    
    std::fstream *aouts = new std::fstream[Param::ncpus];
    std::fstream *pouts = new std::fstream[Param::ncpus];

    int i;
    for ( i = 0; i < Param::ncpus; i++ ) {
        if ( Param::align_flag ) {
            std::string filename = Param::out_dir + "/"
                + SpaFiles[ALIGN] + "."
                + std::to_string(i);
            fio::openFile( aouts[i], filename.c_str(), std::ios::out );
        }
        if ( Param::profile_flag ) {
            std::string filename = Param::out_dir + "/"
                + SpaFiles[PROFILE] + "."
                + std::to_string(i);
            fio::openFile( pouts[i], filename.c_str(), std::ios::out );
        }
    }
    
#pragma omp parallel for schedule(dynamic, 1) if(Param::ncpus>1) private(i) num_threads(Param::ncpus)    
    for ( i = 0; i < (int)npaths; i++ ) {
        size_t num = i;
        PathId pid = (*Id2Paths)[i];

        if ( (*bad_paths)[i] ) continue;

        int thread = omp_get_thread_num();

        if ( Param::align_flag ) {
            aouts[thread] << "ID:" << i << "\n";
            MSA msa( &(*Alns)[i], seqs) ;
            msa.print(aouts[thread],  Param::line_length);
            aouts[thread] << "//\n";
            
            if ( !Param::profile_flag ) continue;
            Profile pro(&msa, AA_TYPE);
            std::string con = pro.makeConsensus();
            pouts[thread] << "ID:" << i << "\n";
            pouts[thread] << "Consensus:" << con << "\n";
            pro.print(pouts[thread]);
            (*Alns)[i].setConsensus(con);
            (*Alns)[i].setDirty(false);
            pouts[thread] << "//\n";
            
            if ( pro.replacedStopCodon() ) {
                if ( Param::verbose > 1 ) 
                    std::cout << "Id:" << num 
                              << "\tpid:" << pid 
                              << "\tstop codon replaced\n";
            }
        } 
        // If there is any dirty path alignment, it needs to regenerate MSA and profile
        else {
            // already updated
            if ( (*Alns)[i].getDirty() == false ) continue;
            MSA msa( &(*Alns)[i], seqs) ;
            Profile pro(&msa, AA_TYPE);
            std::string con = pro.makeConsensus();
            (*Alns)[i].setDirty(false);
        }
    }

    if ( Param::align_flag ) 
        for ( i = 0; i < Param::ncpus; i++ ) aouts[i].close();
    if ( Param::profile_flag ) 
        for ( i = 0; i < Param::ncpus; i++ ) pouts[i].close();
    delete[] aouts;
    delete[] pouts;
    
    if ( Param::align_flag )   combine( SpaFiles[ALIGN] );
    if ( Param::profile_flag ) combine( SpaFiles[PROFILE] );

}

void Reporter::writePlacement()
{
    std::fstream *outs = new std::fstream[Param::ncpus];
    
    int i;
    for ( i = 0; i < Param::ncpus; i++ ) {
        std::string filename = Param::out_dir + "/"
            + SpaFiles[PLACE_TXT] + "."
            + std::to_string(i);
        fio::openFile( outs[i], filename.c_str(), std::ios::out );
    }
    
#pragma omp parallel for schedule(dynamic, 1) if(Param::ncpus>1) private(i) num_threads(Param::ncpus)    
    for ( i = 0; i < (int)npaths; i++ ) {
        // debugging
        if ( Param::debug_id != -1 && Param::debug_id != (*Id2Paths)[i] ) continue;
        if ( (*bad_paths)[i] ) continue;
        
        int thread = omp_get_thread_num();
        
        std::string raw = (*Alns)[i].getSequence();
        std::string con = (*Alns)[i].getConsensus();
        outs[thread] << "ID:" << i << "\n";
        outs[thread] << "Assembled:" << raw << "\n";
        if ( con != "" ) outs[thread] << "Consensus:" << con << "\n";
        (*Alns)[i].printPlacement(outs[thread]);
        outs[thread] << "//\n";
    }
    
    for ( i = 0; i < Param::ncpus; i++ ) outs[i].close();
    delete[] outs;
    
    combine( SpaFiles[PLACE_TXT] );
}

void Reporter::writeSequences()
{
    std::string file = Param::out_dir + "/" + SpaFiles[SEQUENCE];
    std::fstream out;
    fio::openFile( out, file.c_str(), std::ios::out );

    PathLengthMap plen_map = getPathLengths();
    PathLengthMap::reverse_iterator it;
    for ( it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        PathId pid = it->second;
        size_t num = (*Path2Ids)[pid];
        if ( (*bad_paths)[num] ) continue;

        std::string ref = (*Alns)[num].getSequence();
        std::string con = (*Alns)[num].getConsensus();
        std::string seq = (con=="") ? ref : con;

        // remove gaps
        for ( int i = (int)seq.size()-1; i>=0; i-- )
            if ( seq[i] == '-' ) seq.erase(i,1);
        
        if ( Param::clip_stop ) {
            size_t pos = seq.find('*');
            if ( pos != std::string::npos && pos < seq.size()-1 ) { 
                std::string old = seq;
                seq = seq.substr(0, pos);
                if ( Param::verbose >= 1 ) 
                    std::cout << "Id:" << num 
                              << "\tpid:" << pid 
                              << "\tsequence clipped\n";
                if ( Param::verbose ) {
                    std::cout << "old:" << old << "\n";
                    std::cout << "new:" << seq << "\n";
                }
            }
        }
        int rsize = (*Alns)[num].getSize();
        
        out << ">spa" << num;
        out << " len:" << seq.size() 
            << " reads:" << rsize << "\n";
        out << seq << "\n";
    }
}

PathLengthMap Reporter::getPathLengths()
{
    PathLengthMap plen_map;
    for ( size_t i = 0; i < npaths; i++ ) {
        if ( (*bad_paths)[i] ) continue;
        PathId pid = (*Id2Paths)[i];
        size_t len = (*Alns)[i].getSequence().size();
        plen_map.insert( std::pair<size_t, PathId>( len, pid ) );
    }    
    return plen_map;
}




void Reporter::combine( std::string name )
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

