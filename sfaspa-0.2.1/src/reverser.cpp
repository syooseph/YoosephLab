/** 
 * @file reverser.cpp
 * @date Thu 2013-09-26 12:06:39 PM
 * @author Youngik Yang
 * @version 0.002
 * @brief Reverse translation of protein sequences to DNA sequences
 * @details 
 * This modules takes assembled protein sequences and reconstruct
 * DNA sequences. To peform this, DNA profile is generated from
 * short peptide read placements to a protein sequences. Then, 
 * DNA consensus sequence is generated from the profile.
 */

#include "reverser.h"

Reverser::Reverser( std::string ffn,
                    std::string faa,
                    std::string aln,
                    std::string out,
                    int cpu,
                    bool a,
                    bool p,
                    bool v )
{
    ffn_file = ffn;
    faa_file = faa;
    aln_file = aln;
    out_dir  = out;
    ncpus    = cpu;

    align_flag   = a;
    profile_flag = p;
    verbose      = v;

    gap_ext = GAP_EXT;
    gap_open = GAP_OPEN;

    banded_align = true;
    band_ratio = BAND_RATIO;
    //lower_band = LOWER_BAND;
    //upper_band = UPPER_BAND;

    orfs = NULL;
    dnas = NULL;
    init_time = mytime();

    outs = new std::fstream[ncpus];
    logs = new std::fstream[ncpus];
    if ( align_flag )   alns = new std::fstream[ncpus];
    if ( profile_flag ) pros = new std::fstream[ncpus];
    initFiles();

    //run();
}

Reverser::~Reverser()
{
    if ( orfs != NULL ) seq::purge( orfs, norfs );
    if ( dnas != NULL ) seq::purge( dnas, norfs );

    orfs = NULL;
    dnas = NULL;

    delete[] outs;
    delete[] logs;
    if ( align_flag )   delete[] alns;
    if ( profile_flag ) delete[] pros;
}

void Reverser::run()
{   
    setbuf(stdout, NULL); // no buffering

    loadPlacements();

    loadSequences();

    trimDnaReads();

    generate();

    combine();

    printElapsed( init_time, mytime(), "Done" );     
}

void Reverser::loadPlacements()
{
    Param param;
    Param::verbose = verbose;
    mpaths.load( aln_file.c_str() );
}

void Reverser::loadSequences()
{
    norfs = seq::getSequenceCount( faa_file.c_str() );
    if ( verbose ) std::cout << "#orfs:" << norfs << "\n";

    orfs = new char *[norfs];
    int count = 0;
    seq::loadSequenceFile( faa_file, NULL, orfs, count, SEQONLY );
    if ( verbose ) std::cout << "ORFs loaded\n";
    assert(count == norfs);

    count = 0;
    dnas = new char *[norfs];
    seq::loadSequenceFile( ffn_file, NULL, dnas, count, SEQONLY );
    if ( verbose ) std::cout << "DNA read loaded\n";
    assert(count == norfs);
    
    printElapsed( init_time, mytime(), "Sequences loaded\n" );
}


void Reverser::trimDnaReads()
{
    std::cout << "Trimming Reads ...\n";
    setbuf(stdout, NULL); // no buffering

    size_t i;// count = 0;

    progress = Progress(1.0, 0, (size_t)norfs, mytime());

#pragma omp parallel for schedule(dynamic, 1) if(ncpus>1) private(i) num_threads(ncpus)
    for ( i = 0; i < (size_t)norfs; i++ ) {
        trimSingleRead(i);
#pragma omp critical 
{
        ++progress.count;
        progress.showProgress();

    }
}
    printElapsed( init_time, mytime(), "DNA read trimmed\n" );
}


void Reverser::trimSingleRead(int i)
{
    if ( strlen(orfs[i]) == 0 ) {
        std::cout << "[Warning] empty faa:" << i << "\n";
        return;
    }
    
    if ( strlen(dnas[i]) == 0 ) {
        std::cout << "[Warning] empty ffn:" << i << "\n";
        return;
    }
    
    std::string faa = orfs[i];
    std::string ffn = dnas[i];
    std::string orf;
    size_t found = std::string::npos;
    
    
    if ( ffn.size() < faa.size()*3 ) {
        printf("[Warning] ffn-size(%zu) < faa-size(%zu)*3\n", ffn.size(), faa.size());
        return;
    }

    
    while ( ffn.size() >= faa.size()*3 ) {
        orf = codon.translate(ffn);            
        if ( orf.size() < faa.size() ) break;
        found = orf.find(faa);
        if (found != std::string::npos) break;
        ffn.erase(0, 1);
        if ( ffn.size() < 3 ) break; // not a codon
    }
    
    if ( found != std::string::npos ) {     
        if ( ffn == dnas[i] ) return;
        delete[] dnas[i];
        dnas[i] = new char[ faa.size()*3+1];
        strncpy( dnas[i], ffn.c_str(), faa.size()*3 );
        dnas[i][faa.size()*3] = '\0';
    }
    // Just leave it. Not able to handle this.
    else {
        std::cout << "[Warning] match does not found:" << i << "\n";
        std::cout << "\tfaa:" << orfs[i] << "\n";
        std::cout << "\tffn:" << dnas[i] << "\n";
    }
}

void Reverser::initFiles()
{
    std::string file;
    for ( int i = 0; i < ncpus; i++ ) {
        file = out_dir + "/spa.nuc.fasta." + std::to_string(i);
        fio::openFile( outs[i], file.c_str(), std::ios::out );

        file = out_dir + "/spa.nuc.log." + std::to_string(i);
        fio::openFile( logs[i], file.c_str(), std::ios::out );
        
        if ( align_flag ) {
            file = out_dir + "/spa.align.nuc." + std::to_string(i);
            fio::openFile( alns[i], file.c_str(), std::ios::out );
        }
        if ( profile_flag ) {
            file = out_dir + "/spa.profile.nuc." + std::to_string(i);
            fio::openFile( pros[i], file.c_str(), std::ios::out );
        }
    }
}

void Reverser::closeFiles()
{
    for ( int i = 0; i < ncpus; i++ ) {
        outs[i].close();
        logs[i].close();
        if ( align_flag ) alns[i].close();
        if ( profile_flag ) pros[i].close();
    }
}

void Reverser::generate()
{
    std::cout << "Generating nucleotide sequences ...\n";

    size_t npaths = mpaths.getAlignments()->size();
    NucAlns = PathAligners(npaths, PathAligner());

    PathAligners *Alns = mpaths.getAlignments();
    NumToPathMap *Pids = mpaths.getPathIds();
    FlagArray    *bads = mpaths.getBadPaths();

    size_t i, n = (*Alns).size(), count = 0;

    progress = Progress( 1.0, count, n, mytime() );

#pragma omp parallel for schedule(dynamic, 1) if(ncpus>1) private(i) num_threads(ncpus)
    for ( i = 0; i < n; i++ ) {
#pragma omp atomic
        ++progress.count;

        int thread = omp_get_thread_num();

        PathId pid = (*Pids)[i];
        if ( Param::verbose ) 
            std::cout << "PathId:" << pid << "\tthread:" << thread << "\n";
        if ( !(*bads)[i] ) 
            rtranslate( thread, i, (*Alns)[i], NucAlns[i] );
#pragma omp critical 
{
    progress.showProgress();
}
    }

    closeFiles();

    printElapsed( init_time, mytime(), "Nucleotide sequences generated" );         
}

void Reverser::rtranslate( int thread,
                           int nid,
                           PathAligner &AacAln,
                           PathAligner &NucAln )
{
    std::vector<std::string> DnaReads;
    ReadPlacementList *places = AacAln.getPlacements();
    std::string spa_raw       = AacAln.getSequence();
    std::string spa_con       = AacAln.getConsensus();
    std::string spa_seq = (spa_con != "" ) ? spa_con : spa_raw;
    if ( Param::verbose ) std::cout << "seq:" << spa_seq << "\n";
    if ( spa_seq.size() == 0 ) {
        std::cout << "Zero length path:" << nid << "\n";
        exit(EXIT_FAILURE);
    }

    size_t nuc_len = spa_seq.size() * 3;
    std::string nuc_seq = std::string(nuc_len, '.');
    NucAln.setSequence(nuc_seq);

    if ( places->size() == 0 ) {
        std::cerr << "Zero read supported path:" << nid << "\n";
        exit(EXIT_FAILURE);
    }

    ReadPlacementList nuc_places;
    for ( ReadPlacementList::iterator it = places->begin(); it != places->end(); ++it ) {
        ReadPlacement aac_place = *it;
        ReadPlacement nuc_place;// = *it;
        nuc_place.rid      =   aac_place.rid;
        nuc_place.read_pos = 3*aac_place.read_pos;
        nuc_place.path_pos = 3*aac_place.path_pos;
        nuc_place.length   = 3*aac_place.length;
        
        std::list<Mismatch>::iterator mt;
        for ( mt = it->ilist.begin(); mt != it->ilist.end(); ++ mt ) {
            for ( size_t i = 0; i < 3; i++ ) {
                Mismatch mm = *mt;
                mm.qry_pos = 3*mt->qry_pos + i;
                mm.ref_pos = 3*mt->ref_pos;
                nuc_place.ilist.push_back(mm);
            }
        }
        for ( mt = it->dlist.begin(); mt != it->dlist.end(); ++ mt ) {
            for ( size_t i = 0; i < 3; i++ ) {
                Mismatch mm = *mt;
                mm.qry_pos = 3*mt->qry_pos;
                mm.ref_pos = 3*mt->ref_pos + i;
                nuc_place.dlist.push_back(mm);
            }
        }
        nuc_places.push_back(nuc_place);
    }
    
    NucAln.setPlacements(nuc_places);

    NucAln.setLStop(AacAln.getLStop());
    NucAln.setRStop(AacAln.getRStop());
    NucAln.setLTrim(AacAln.getLTrim());
    NucAln.setRTrim(AacAln.getRTrim());


    MSA msa( &NucAln, dnas) ;    

    Profile profile(&msa, DNA_TYPE);    
    std::string dna_con = profile.makeConsensus();
    assert(dna_con.size());
    msa.setPivot(dna_con);

    if ( align_flag ) {
        alns[thread] << "ID:" << nid << "\n";
        msa.print(alns[thread],  Param::line_length);
        alns[thread] << "//\n";
    }

    if ( profile_flag ) {
        pros[thread] << "ID:" << nid << "\n";
        pros[thread] << "Consensus:" << dna_con << "\n";        
        profile.print(pros[thread]);
        pros[thread] << "//\n";
    }


    for ( int i = (int)dna_con.size()-1; i>=0; i-- ) 
        if ( dna_con[i] == '-' ) dna_con.erase(i,1);
    
    assert(dna_con.size());
    NucAln.setSequence(dna_con);

    std::string transeq = codon.translate(dna_con);            

    if ( verbose ) {
        logs[thread] << "Consensus:" << dna_con << "\n";
        logs[thread] << "Profile:\n";
        profile.print(logs[thread]);
    }

    logs[thread] << "ID:" << nid << "\n";
    logs[thread] << "SPASequence:" << spa_seq << "\n";
    logs[thread] << "Nucleotides:" << dna_con << "\n";
    logs[thread] << "Translation:" << transeq << "\n";
    
    double positive = 0;
    double identity = 0;
    if ( transeq == spa_seq ) {
        positive = identity = 100;
    } else {
    int band_size = int ( band_ratio * spa_seq.size() );
        SemiGlobalAlign galn = banded_align ? 
            SemiGlobalAlign(spa_seq, transeq, ANCHOR_CENTER) :
            SemiGlobalAlign(spa_seq, transeq, ANCHOR_CENTER, gap_ext, gap_open, -1*band_size, band_size );
        AlignSummary summary = galn.getSummary();
        positive = summary.posrate * 100;
        identity = (double)summary.match/summary.length * 100;
    }

    //if ( verbose ) {
    char buff[100];
    sprintf(buff, "%.2f", identity);
    logs[thread] << "identity:" << buff << "\n";
    sprintf(buff, "%.2f", positive);
    logs[thread] << "positive:" << buff << "\n";
    logs[thread] << "//\n";
    //}
    
    outs[thread] << ">spa" << nid << "\n" << dna_con << "\n";
}


void Reverser::combine()
{
    std::cout << "Writing output ...\n";
    std::string nseq = out_dir + "/spa.nuc.fasta";
    std::string nlog = out_dir + "/spa.nuc.log";
    std::string naln = out_dir + "/spa.align.nuc";
    std::string npro = out_dir + "/spa.profile.nuc";
    if ( ncpus == 1 ) {
        std::string oseq = out_dir + "/spa.nuc.fasta.0";
        rename( oseq.c_str(), nseq.c_str() );
        std::string olog = out_dir + "/spa.nuc.log.0";
        rename( olog.c_str(), nlog.c_str() );
        if ( align_flag ) {
            std::string oaln = out_dir + "/spa.align.nuc.0";
            rename( oaln.c_str(), naln.c_str() );
        }
        if ( profile_flag ) {
            std::string opro = out_dir + "/spa.profile.nuc.0";
            rename( opro.c_str(), npro.c_str() );
        }
        return;
    }
    
    std::fstream out;
    fio::openFile( out, nseq.c_str(), std::ios::out );
    std::fstream log;
    fio::openFile( log, nlog.c_str(), std::ios::out );
    std::fstream aln;
    if ( align_flag ) fio::openFile( aln, naln.c_str(), std::ios::out );
    std::fstream pro;
    if ( profile_flag ) fio::openFile( pro, npro.c_str(), std::ios::out );

    std::string line;
    for ( int i = 0; i < ncpus; i++ ) {
        std::stringstream oseq;
        oseq << out_dir << "/spa.nuc.fasta." << i;
        std::ifstream sin( oseq.str().c_str() );
        while ( getline(sin, line ) ) 
            out << line << "\n";

        if ( remove( oseq.str().c_str() ) != 0 ) 
            perror( "Error deleting file" );

        std::stringstream olog;
        olog << out_dir << "/spa.nuc.log." << i;
        std::ifstream lin( olog.str().c_str() );
        while ( getline(lin, line ) ) 
            log << line << "\n";

        if ( remove( olog.str().c_str() ) != 0 ) 
            perror( "Error deleting file" );

        if ( align_flag ) {
            std::stringstream oaln;
            oaln << out_dir << "/spa.align.nuc." << i;
            std::ifstream ain( oaln.str().c_str() );
            while ( getline(ain, line ) ) 
                aln << line << "\n";
            
            if ( remove( oaln.str().c_str() ) != 0 ) 
                perror( "Error deleting file" );
        }

        if ( profile_flag ) {
            std::stringstream opro;
            opro << out_dir << "/spa.profile.nuc." << i;
            std::ifstream pin( opro.str().c_str() );
            while ( getline(pin, line ) ) 
                pro << line << "\n";
            
            if ( remove( opro.str().c_str() ) != 0 ) 
            perror( "Error deleting file" );
        }
    }

    out.close();
    log.close();
    if ( align_flag ) aln.close();
    if ( profile_flag ) pro.close();

    printElapsed( init_time, mytime(), "Output written" );     
}

