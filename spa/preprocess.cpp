#include "preprocess.h"

int main(int argc, char* argv[]) 
{
    double time1 = mytime();

    setbuf(stdout, NULL); // no buffering

    int k;
    std::vector<std::string> input_files;
    std::string kmer_file, graph_file, index_file, out_dir;
    args::parse_cmd_args_prepare(argc, argv, k, input_files, graph_file, index_file, out_dir);

    printLocalTime();
    args::printCommand(argc, argv);

    //---------------
    // Sequence count
    //---------------
    int nreads = seq::totalSequenceCount( input_files );
    std::cout << "#seqs:" << nreads << "\n";

    //------------------
    // IDs and sequences
    //------------------
    char **seqs;
    seqs = new char *[nreads];

    
    //----------------------
    // load sequence file(s)
    //----------------------
    seq::loadSequences(input_files, NULL, seqs, SEQONLY);


    //-----------------------
    // Find invalid sequences
    //-----------------------
    std::set<int> bad_index;
    trim( bad_index, seqs, nreads, k );


    //-------------------------------------
    // Amino acids in front of current kmer
    //-------------------------------------
    LeftAAsCoverageMap prev_AAs;
    CoverageMap coverage_map;
    setKmerLinks( prev_AAs, coverage_map, seqs, nreads, bad_index, k );
    
    //--------------------------
    // DeBruijn graph input dump
    //--------------------------
    writeGraphInput( graph_file.c_str(), prev_AAs, coverage_map );

    //------------------------
    // Now release graph input
    //------------------------ 
    trash( prev_AAs );


    //---------------------
    // Inverted index build
    //---------------------
    InvertedIndex InvInd;
    buildIndex( InvInd, seqs, nreads, coverage_map, bad_index, k );
    InvInd.write( index_file.c_str() );
    
    cov::coverageSummary(InvInd, std::cout);

    InvInd.clear();

    
    //---------------
    // Reclaim memory
    //---------------
    seq::purge( seqs, nreads );


    std::cout << "\nPreprocessing:" << mytime()-time1 << " sec\n";
    return 0;
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
    double time1 = mytime();

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
        if ( n > 0 && n % 100000 == 0 ) 
            std::cout << "\t" << n << " seqs\t(" << mytime()-time1 << " sec)\n";
    }
    
    std::cout << "graph input created:\t(" << mytime()-time1 << " sec)\n\n";
}


void trim( std::set<int> &bad_index,
           char **seqs, 
           int size,
           int k )
{
    double time1 = mytime();

    for ( int i = 0; i < size; i++ ) {
        std::string seq = seqs[i];
        if ( (int)seq.size() < k ) {
            bad_index.insert(i);
            std::cout << i+1 << ": Short\n";
        }
        if ( !seq::validSequence(seq) ) {
            bad_index.insert(i);
            std::cout << i+1 << ": Invalid\n";
        }
    }
    std::cout << "no. of invalid sequences:" << bad_index.size() << "\t(";
    std::cout << mytime()-time1 << " sec)\n\n";
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

/**
 * \pre Sequence lengths are greater than k-mer.
 */
void buildIndex( InvertedIndex &InvInd, 
                 char **seqs, 
                 int nreads, 
                 CoverageMap &coverage_map,
                 std::set<int> &bad_index,
                 int k)
{
    double time1 = mytime();

    int i;
    for ( i = 0; i < nreads; i++ ) {
        if ( bad_index.find(i) != bad_index.end() ) continue;
        ReadId rid = (ReadId)i;
        KmerSet kmers = getKmerSet(seqs[i], k);
        Record r( rid, kmers );
        InvInd.add(r, coverage_map);
        if ( i && i%100000 == 0 ) std::cout << "\t" << i << " seqs\t(" << mytime()-time1 << " sec)\n";
    }
    std::cout << "index built:\t(" << mytime()-time1 << " sec)\n\n";
}


void writeGraphInput( const char* file,
                      LeftAAsCoverageMap &prev_AAs,
                      CoverageMap &coverage_map )
{
    std::fstream out;
    io::openFile( out, file, std::ios::out | std::ios::binary );
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

//
//void computeReadKmerStat(char **seqs,
//                         InvertedIndex &iindex,
//                         std::set<int> &bad_index,
//                         int nreads,
//                         int k,
//                         std::fstream &out)
//{
//    out << "No\t#kmers\tMax\tMin\tMean\tMedian\n";
//    for ( int n = 0; n < nreads; n++ ) {
//        if ( bad_index.find(n) != bad_index.end() ) continue;
//        std::string seq = seqs[n];
//        size_t nkmer = seq.size() - k + 1;
//        std::vector<double> supports(nkmer);
//        for (size_t i = 0; i < nkmer; i++) {
//            KmerType kmer = seq.substr(i,k);
//            KmerId   kmer_id = alpha::AminoAcidToInteger<KmerId>(kmer);
//
//            if ( ! iindex.has(kmer_id) ) {
//                std::cout << "Index does not exist:" << kmer << "\n";
//                exit(1);
//            }
//            supports[i] = iindex.getValue(kmer_id)->size;
//        }
//        out << n+1 << "\t"
//                  << nkmer << "\t"
//                  << math::max(&supports[0], nkmer) << "\t"
//                  << math::min(&supports[0], nkmer) << "\t"
//                  << math::mean(&supports[0], nkmer) << "\t"
//                  << math::median(&supports[0], nkmer, false)
//                  << "\n";
//    }
//}

