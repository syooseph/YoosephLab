#include "preprocess.h"

int main(int argc, char* argv[]) 
{
    double time1 = mytime();

    setbuf(stdout, NULL); // no buffering

    std::string log;

    int k, nparts;
    std::vector<std::string> input_files;
    std::string out_dir;
    bool reverse_sfa = false;
    bool sfa_only = false;
    bool verbose = false;
    args::parse_cmd_args_prepare(argc, argv, k, input_files, out_dir, nparts, reverse_sfa, sfa_only, verbose );

    setVerbose(verbose);

    printLocalTime();
    args::printCommand(argc, argv, std::cout);

    //---------------
    // Sequence count
    //---------------
    int nreads = seq::totalSequenceCount( input_files );
    //std::cout << "#seqs:" << nreads << "\n";

    //------------------
    // IDs and sequences
    //------------------
    char **seqs = new char *[nreads];
    
    //----------------------
    // load sequence file(s)
    //----------------------
    //seq::loadSequences(input_files, tags, seqs, TAGSEQ);
    seq::loadSequences(input_files, NULL, seqs, SEQONLY);

    //log = boost::lexical_cast<std::string>(nreads) + " sequences loaded"; 
    log = std::to_string(nreads) + " sequences loaded"; 
    printElapsed( time1, mytime(), log.c_str() );

    //-------------------
    // Build suffix array
    //------------------- 
    buildSuffixArray( seqs, nreads, nparts, out_dir );

    //log = boost::lexical_cast<std::string>(nparts) + " suffix array created";
    log = std::to_string(nparts) + " suffix array created";
    printElapsed( time1, mytime(), log.c_str() );
 
    char **rseqs;
    if ( reverse_sfa ) {
        rseqs = new char *[nreads];
        reverseReads( rseqs, seqs, nreads );
        buildReverseSuffixArray( rseqs, nreads, nparts, out_dir );
    }

    //-------------------------------------
    // Amino acids in front of current kmer
    //-------------------------------------
    if ( !sfa_only ) {
        //-----------------------
        // Find invalid sequences
        //-----------------------
        std::set<int> bad_index;
        trim( bad_index, seqs, nreads, k );
        std::string str = "Bad sequence identified (count:" + 
            //boost::lexical_cast<std::string>( bad_index.size() ) + ")";
            std::to_string( bad_index.size() ) + ")";
        printElapsed( time1, mytime(), str.c_str() );
        
        LeftAAsCoverageMap prev_AAs;
        CoverageMap coverage_map;
        setKmerLinks( prev_AAs, coverage_map, seqs, nreads, bad_index, k );

        //--------------------------
        // DeBruijn graph input dump
        //--------------------------
        std::string graph_file = out_dir + "/gin";
        writeGraphInput( graph_file.c_str(), prev_AAs, coverage_map );
        printElapsed( time1, mytime(), "Graph created" );
         
        if ( bad_index.size() ) {
            std::string bad_file = out_dir + "/bad";
            writeBadIndex( bad_file.c_str(), bad_index );
            printElapsed( time1, mytime(), "Bad read ID written" );
        }
        
        //------------------------
        // Now release graph input
        //------------------------ 
        trash( prev_AAs );
    }

    //---------------
    // Reclaim memory
    //---------------
    seq::purge( seqs, nreads );
    if ( reverse_sfa ) seq::purge( rseqs, nreads );

    printElapsed( time1, mytime(), "Preprocessing completed" );
    return 0;
}

void reverseReads( char **rseqs,
                   char **seqs,
                   int nreads )
{
    for ( int i = 0; i < nreads; i++ ) {
        int len = strlen(seqs[i]);
        //char *tmp = new char[len+1];
        std::string tmp(seqs[i]);
        std::string rev = std::string( tmp.rbegin(), tmp.rend() );
        //strcpy( tmp, seqs[i] );
        //strrev( tmp );
        rseqs[i] = new char[len+1];
        strcpy( rseqs[i], rev.c_str() );
        //delete[] tmp;
    }
}

void buildSuffixArray( char **seqs, 
                       int nreads,
                       int nparts,
                       //int ncpus,
                       std::string out_dir)
{
    for ( int i = 0; i < nparts; i++ ) {
        int sr = i*int(nreads/nparts);
        int nr = (i < nparts-1) ? nreads/nparts : nreads - int(nreads/nparts)*(nparts-1);
        // std::string lcp_file = out_dir + "/lcp." + boost::lexical_cast<std::string>(i);
        // std::string mcp_file = out_dir + "/mcp." + boost::lexical_cast<std::string>(i);
        // std::string gsa_file = out_dir + "/gsa." + std::to_string(i);
        std::string lcp_file = out_dir + "/lcp." + std::to_string(i);
        std::string mcp_file = out_dir + "/mcp." + std::to_string(i);
        std::string gsa_file = out_dir + "/gsa." + std::to_string(i);
        
        GSA gsa( &seqs[sr],
                 nr, 
                 true,
                 gsa_file.c_str(),
                 lcp_file.c_str(),
                 mcp_file.c_str() );
    }
}

void buildReverseSuffixArray( char **rseqs, 
                              int nreads,
                              int nparts,
                              //int ncpus,
                              std::string out_dir)
{
    for ( int i = 0; i < nparts; i++ ) {
        int sr = i*int(nreads/nparts);
        int nr = (i < nparts-1) ? nreads/nparts : nreads - int(nreads/nparts)*(nparts-1);
        std::string lcp_file = out_dir + "/rlcp." + std::to_string(i);
        std::string mcp_file = out_dir + "/rmcp." + std::to_string(i);
        std::string gsa_file = out_dir + "/rgsa." + std::to_string(i);
        
        GSA gsa( &rseqs[sr],
                 nr, 
                 true,
                 gsa_file.c_str(),
                 lcp_file.c_str(),
                 mcp_file.c_str() );
    }
}

void trash( LeftAAsCoverageMap &prev_AAs ) 
{
    LeftAAsCoverageMap::iterator it;
    for ( it = prev_AAs.begin(); it != prev_AAs.end(); ++it ) 
        delete[] prev_AAs[it->first];

    prev_AAs.clear();
}

void setKmerLinks( LeftAAsCoverageMap &prev_AAs,
                 CoverageMap &coverage_map,
                char **seqs,
                 int nreads,
                 std::set<int> &bad_index,
                 int k)
{
    for ( int n = 0; n < nreads; n++ ) {
        if ( bad_index.find(n) != bad_index.end() ) continue;

        std::string seq = seqs[n];
        assert((int)seq.size() >= k);

        for (int i = 0; i <= (int)seq.size() - k; i++) {
            KmerType kmer = seq.substr(i,k);
            KmerId   kmer_id = alpha::AminoAcidToInteger<KmerId>(kmer);
            if ( prev_AAs.find(kmer_id) == prev_AAs.end() ) {
                prev_AAs[kmer_id] = new CoverageType[ alpha::COUNT_AA ];
                for ( int j = 0; j < alpha::COUNT_AA; j++ )
                    prev_AAs[kmer_id][j] = 0;
            }            
            if ( i > 0 ) {
                unsigned aa_num = alpha::AsciiToAA[(int)seq[i-1]] - 1;
                prev_AAs[kmer_id][aa_num] += 1;
            }

            if ( coverage_map.find(kmer_id) ==  coverage_map.end() )
                coverage_map.insert(std::pair<KmerId, CoverageType>(kmer_id, 0));
            coverage_map[kmer_id]++;
            
        }
    }
}


void trim( std::set<int> &bad_index,
           char **seqs, 
           int size,
           int k )
{
    int ctshort = 0, ctbadch = 0;
    for ( int i = 0; i < size; i++ ) {
        std::string seq = seqs[i];
        if ( (int)seq.size() < k ) {
            ctshort++;
            bad_index.insert(i);
            if ( VERBOSE ) std::cout << i+1 << ": Short\n";
        }
        if ( !seq::validSequence(seq) ) {
            bad_index.insert(i);
            ctbadch++;
            if ( VERBOSE ) std::cout << i+1 << ": Invalid\n";
        }
    }
    if ( VERBOSE ) 
        std::cout << "no. of invalid sequences:" << bad_index.size() 
                  << " (short:" << ctshort << ", bad char:" << ctbadch << ")\n";
}



/**
 * \pre The length of a sequence is at least the size of k-mer length.
 */
KmerSet getKmerSet( char *cseq, int k )
{
    std::string seq = cseq;
    std::list<KmerId> kids;
    for (size_t i = 0; i <= seq.size() - k; i++) {
        KmerType kmer = seq.substr(i,k);
        kids.push_back( alpha::AminoAcidToInteger<KmerId>(kmer) );
    }
    return KmerSet( kids.begin(), kids.end() );
}

void writeGraphInput( const char* file,
                      LeftAAsCoverageMap &prev_AAs,
                      CoverageMap &coverage_map )
{
    std::fstream out;
    fio::openFile( out, file, std::ios::out | std::ios::binary );
    LeftAAsCoverageMap::iterator it;
    for ( it = prev_AAs.begin(); it != prev_AAs.end(); ++it ) {
        out.write((char*)&(it->first), sizeof(KmerId));
        out.write((char*)&coverage_map[it->first], sizeof(CoverageType));
        for ( int i = 0; i < alpha::COUNT_AA; i++ ) {
            out.write((char*)&(prev_AAs[it->first][i]), sizeof(CoverageType));
        }
    }
    out.close();
}

void writeBadIndex( const char* file,
                    std::set<int> &bad_index )
{
    std::fstream out;
    fio::openFile( out, file, std::ios::out | std::ios::binary );
    size_t nbad = bad_index.size();
    out.write((char*)&nbad, sizeof(size_t));
    std::set<int>::iterator it;
    for ( it = bad_index.begin(); it != bad_index.end(); ++it )
        out.write((char*)&(*it), sizeof(int));
    out.close();
}
