/** 
 * @file spa.app
 * @date Thu 2011-10-20 12:39:09 PM
 * @author Youngik Yang
 * @version 0.001
 * @brief Main function of SPA assembler.
 * @details 
 * This is main handler of SPA assembler. 
 * It takes sequences and generates ORFs as output
 */

#include "main.h"

int main( int argc, char* argv[] ) 
{	
    /*===============*/
    /* Local objects */
    /*===============*/
	DeBruijnGraph graph;            // k-mer de Bruijn graph
    CoverageMap   kmer_coverage;    // kmer-id to frequency mapping
    VertexToKmerMap vertex_map;     // vertex to kmer-id mapping
    std::set<KmerType> debug_kmers; // kmers of interest
    PathToAlnMap path2aln_map;      // path to alignment mapping
    InvertedIndex iindex;           // Inverted index
    int nreads;                     // Total number of reads
    BitString *bstrs;               // Reads as bitstream
    ReadId *pairs;                  // Read pairing
    char *strands;                  // Read strands
    PathId *used_reads;             // Path Ids where each read belongs to
	Param param;                    // Parameters

    double t0 = mytime();
    setbuf(stdout, NULL); // no buffering

    /* Load data */
	load( argc, argv, bstrs, strands, pairs, nreads, used_reads, debug_kmers, param );

	/* Build graph & index */
	build( graph, vertex_map, kmer_coverage, path2aln_map, iindex, nreads, used_reads, param );

	/* Assemble reads */
	assemble( graph, iindex, path2aln_map, vertex_map, kmer_coverage, bstrs, strands, pairs, nreads, used_reads, debug_kmers, param );

	/* Write results */
	write( graph, iindex, path2aln_map, vertex_map, kmer_coverage, bstrs, strands, pairs, used_reads, nreads, param );

    /* Clean up */
    release( path2aln_map, iindex, strands, bstrs, used_reads, pairs, param );

    std::cout << "\nTotal time elapsed:" << mytime()-t0 << " sec\n";

    return 0;
}

void load( int argc, 
		   char *argv[], 
		   BitString *&bstrs, 
		   char *&strands, 
		   ReadId *&pairs, 
		   int &nreads, 
		   PathId *&used_reads, 
		   std::set<KmerType> &debug_kmers, 
		   Param &param )
{
    //--------------------
    // Get program options
    //--------------------
    readParams(argc, argv, param);

    //--------------------------
    // Load sequence information
    //--------------------------
    loadReads( bstrs, strands, pairs, nreads, param );

    //-------------------------
    // Reads to Path membership
    //-------------------------
    initReadMembership(used_reads, nreads);

    //--------------------------------
    // Kmers of interest for debugging
    //--------------------------------
   if ( param.debug_file != "" ) loadDebugKmers( debug_kmers, param );

   
}

/**
 * Load reads.
 * Read in all reads and save them as bit strings.
 * Parse sequence tags and record strands and paired-end information.
 * \param bstrs		An array of bit strings
 * \param strands	An array of sequence strands
 * \param pairs		An array of paired end sequence IDs
 * \param nreads	Integer, number of sequence reads
 * \pre		Empty arrays with size nreads
 * \post	Arrays with values
 */
void loadReads( BitString *&bstrs, 
                char      *&strands, 
                ReadId    *&pairs, 
                int        &nreads,
                Param     &param)
{
    //--------------------------
    // Get total number of reads
    //--------------------------
    nreads = seq::totalSequenceCount( param.input_files );
    if ( param.verbose ) std::cout << "Total sequences:" << nreads << "\n";

    //-------------------------------------------
    // Load sequences and tags from input file(s)
    //-------------------------------------------
    char **tags = new char *[nreads];
    char **seqs = new char *[nreads];
    seq::loadSequences(param.input_files, tags, seqs, TAGSEQ);


    strands = new char[nreads];
    bstrs   = new BitString[nreads];
    setBitstrs( bstrs, seqs, nreads );   // save sequences as bit strings
    setStrands( strands, tags, nreads ); // get strands information

    //-----------------
    // Paired-end reads
    //-----------------
    if ( param.pair_flag ) {
        pairs = new ReadId[nreads];
        setReadPairs( pairs, tags, nreads );
    }

    for ( int i = 0; i < nreads; i++ ) delete[] seqs[i];
    delete[] seqs;

    for ( int i = 0; i < nreads; i++ ) delete[] tags[i];
    delete[] tags;

}

void setStrands( char *strands, char **tags, int nreads )
{
    for ( int i = 0; i < nreads; i++ ) {
        std::string rid = tags[i];
        char s = rid[rid.size()-1];
        if ( s == '+' || s == '-' ) strands[i] = s;
        else strands[i] = '\0';
    }
}

void setBitstrs( BitString *bstrs, char **seqs, int nreads )
{
    for ( int i = 0; i < nreads; i++ ) {
        bstrs[i] = BitString(seqs[i]);
    }
}

void setReadPairs( ReadId *pairs, 
                   char   **tags,
                   int    nreads )
{
    dropPairInfo(tags, nreads);

    //-------------------------------------------
    // initialize each with total number of reads
    //-------------------------------------------
    for ( int i = 0; i < nreads; i++ ) pairs[i] = NOT_PAIR;
    
    for ( int i = 0; i < nreads-1; i++ ) {
        if ( pairs[i] != NOT_PAIR ) continue;
        
        if ( strcmp(tags[i], tags[i+1]) == 0 ) {
            pairs[i] = i+1;
            pairs[i+1] = i;
        }
    }
}

void initReadMembership( PathId *&used_reads, int nreads )
{
    used_reads = new PathId[nreads];
    initReadFlags(used_reads, nreads);
}

void initReadFlags( PathId *used_reads, int nreads )
{
    for ( int i = 0; i < nreads; i++ ) used_reads[i] = NOT_PATH;
}

void dropPairInfo(char **tags, int nreads)
{
    for ( int i = 0; i < nreads; i++ ) {
        std::string rid = tags[i];
        size_t p = rid.find("/");
        if ( p >= rid.size() ) {
            std::cerr << "[Error] Invalid paired end tag\n"; 
            exit(1);
        }
        delete[] tags[i];

        std::string nrid = rid.substr(0, p);
        int nlen = nrid.size()+1;
        tags[i] = new char[nlen];
        strncpy( tags[i], nrid.c_str(), nlen );
    }
}

void build( DeBruijnGraph &graph, 
			VertexToKmerMap &vertex_map,
			CoverageMap &kmer_coverage,
			PathToAlnMap &path2aln_map, 
			InvertedIndex &iindex, 
			int nreads, 
			PathId *used_reads, 
			Param &param )
{
    //-------------------
    // Graph construction
    //-------------------
    if ( param.file_suffix != "" )
        loadDump(graph, vertex_map, kmer_coverage, path2aln_map, used_reads, nreads, param);
    else 
        makeGraph( graph, vertex_map, kmer_coverage, path2aln_map, used_reads, nreads, param );

    //----------------------------------------------
    // Load inverted index (k-mer to a set of reads)
    //----------------------------------------------
    buildIndex(iindex, param);

}

void buildIndex(InvertedIndex &iindex, Param &param)
{
    std::cout << "\nLoading inverted index ...\n";
    double t0 = mytime();
    iindex.load(param.index_file.c_str());
    std::cout << "Inverted index loaded:" << mytime()-t0 << " sec\n";
}

void assemble( DeBruijnGraph &graph, 
			   InvertedIndex &iindex, 
			   PathToAlnMap &path2aln_map, 
			   VertexToKmerMap &vertex_map,
			   CoverageMap &kmer_coverage,
			   BitString *bstrs, 
			   char *strands, 
			   ReadId *pairs, 
			   int nreads, 
			   PathId *used_reads,
			   std::set<KmerType> &debug_kmers,
			   Param &param )
{
	assem::proceed(graph, iindex, path2aln_map, vertex_map, kmer_coverage, bstrs, strands, pairs, nreads, used_reads, debug_kmers, param );

}

void write( DeBruijnGraph &graph, 
			InvertedIndex &iindex,
			PathToAlnMap &path2aln_map, 
			VertexToKmerMap &vertex_map,
			CoverageMap &kmer_coverage,
			BitString *bstrs, 
			char *strands, 
			ReadId *pairs, 
			PathId *used_reads,
			int nreads, 
			Param &param )
{
    //-----------------------------
    // Dump objects to binary files
    //-----------------------------
    if ( param.dump_flag )
		dumpObjects(graph, vertex_map, kmer_coverage, iindex, path2aln_map, param );

    //----------------------------
    // Write results in text files
    //----------------------------
    if (param.output_flag) writeResult(path2aln_map, bstrs, strands, pairs, used_reads, nreads, param );

}

void release( PathToAlnMap &path2aln_map,
              InvertedIndex &iindex,
              char *strands,
              BitString *bstrs,
              PathId *used_reads,
              ReadId *pairs,
              Param &param )
{
    PathToAlnMap::iterator mt;
    for ( mt = path2aln_map.begin(); mt != path2aln_map.end(); ++mt ) 
        delete mt->second;

    iindex.clear();
    delete[] strands;
    delete[] bstrs;
    delete[] used_reads;
    if ( param.pair_flag ) delete[] pairs;
}


void writeResult( PathToAlnMap &path2aln_map,
                  BitString    *bstrs,
                  char         *strands,
                  ReadId       *pairs, 
                  PathId       *used_reads,
                  int           nreads,
                  Param       &param)
{

	/* set verbosity false temporarily */
	bool verbose = param.verbose;
	if ( param.verbose ) param.verbose = false;

    setbuf(stdout, NULL); // no buffering
    double time1 = mytime();
    std::cout << "\nWriting outputs ...\n";

    std::string consensus_file = param.out_dir + "/spa.fasta";
    std::string statistic_file = param.out_dir + "/spa.cover";
    std::string placement_file = param.out_dir + "/spa.place";
    std::string alignment_file = param.out_dir + "/spa.align";
    std::string profile_file   = param.out_dir + "/spa.profile";
    std::fstream cons, stat, plac, alig, prof;

    io::openFile(cons, consensus_file.c_str(), std::ios::out);
    io::openFile(stat, statistic_file.c_str(), std::ios::out);
    io::openFile(plac, placement_file.c_str(), std::ios::out);
    if ( param.profile_flag ) io::openFile(prof, profile_file.c_str(), std::ios::out);
    if ( param.align_flag ) io::openFile(alig, alignment_file.c_str(), std::ios::out);
    
    int count = 0;
    //int nproc = 0;
    //int npath = path2aln_map.size();
    //double pratio = 1;

    //double lt0 = mytime();
//     for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end();  ) {

//         nproc++;
//         if ( nproc/(double)npath*100 >= pratio ) {
//             fprintf( stdout, "%d/%d: %.2f sec\n", nproc, int(npath), mytime()-lt0);
//             pratio += 1;
//         }

//         SpaPath *spath = it->second;
//         writeConsensus( cons, spath, count );
//         writePlacement( plac, spath, count );
//         writeStatistic( stat, spath, count );
//         if ( param.align_flag )   writeAlignment( alig, spath, count, bstrs );
//         if ( param.profile_flag ) writeProfile( prof, spath, count );
//         ++count;
//         ++it;
        
// //         int error = 0;

// //         //---------------------------------------------------------------------------------
// //         // It might be redundant but do drop low covered base and bases after stop codon
// //         // to get consistent output for fasta, coverage, placement, profile, and alignment.
// //         //---------------------------------------------------------------------------------
// // 		MSA msa( *spath, it->first, used_reads, bstrs, strands, pairs, PILE|TRIM|CODON, error, param);
        
// //         if ( error ) {
// //             if ( param.verbose ) std::cerr << "\n** MSA failed\n"; 
// //             spath->resetReads(used_reads, NOT_PATH);
// //             delete spath;
// //             it = path2aln_map.erase(it);
// //         }
// //         else {
// //             writeConsensus( cons, msa, count );
// //             writePlacement( plac, spath, count );
// //             writeStatistic( stat, msa, count );
// //             if ( param.align_flag )   writeAlignment( alig, msa, spath, count );
// //             if ( param.profile_flag ) writeProfile( prof, msa, count );
// //             ++count;
// //             ++it;
// //         }
//     }
    std::cout << "Writing consensus ...\n";
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
        SpaPath *spath = it->second;
        writeConsensus( cons, spath, count );
        count++;
    }
    cons.close();
    
    count = 0;
    std::cout << "Writing read placment ...\n";
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
        SpaPath *spath = it->second;
        writePlacement( plac, spath, count );
        count++;
    }
    plac.close();

    count = 0;
    std::cout << "Writing path statistic ...\n";
    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
        SpaPath *spath = it->second;
        writeStatistic( stat, spath, count );
        count++;
    }
    stat.close();

    if ( param.align_flag ) {
        count = 0;
        std::cout << "Writing alignment ...\n";
        for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
            SpaPath *spath = it->second;
            writeAlignment( alig, spath, count, bstrs );
            count++;
        }
        alig.close();
    }

    if ( param.profile_flag ) {
        count = 0;
        std::cout << "Writing alignment ...\n";
        for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
            SpaPath *spath = it->second;
            writeProfile( prof, spath, count );
            count++;
        }
        prof.close();
    }


//     cons.close(); stat.close(); plac.close();
//     if ( param.align_flag ) alig.close();
//     if ( param.profile_flag ) prof.close();
    
    writeUnusedReads( used_reads, nreads, param );

    std::cout << "Outputs written:" << mytime()-time1 << " sec\n";

	if ( verbose ) param.verbose = true; /* reset */
}

//void writeConsensus( std::fstream &out, MSA &msa, int count )
void writeConsensus( std::fstream &out, SpaPath *spath , int count )
{
    //std::string orf = spath->getConsensus();
    //std::string orf = spath->getConsensusString(); //??? this causing problem??? 
    //assert( orf.size() > 0 );

    //orf = biostr::stripGap(orf);

//     out << ">" << count << "\n";
//     out << orf << "\n";


    char *orf = spath->getConsensus();
    assert( strlen(orf) > 0 ) ;

    out << ">" << count << "\n";
    out << orf << "\n";

    out.flush();
}


void writePlacement( std::fstream &out, SpaPath *spath, int count )
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

void writeStatistic( std::fstream &out, SpaPath *spath, int count )
{
    out << "ID:" << count << "\n";
    //std::string raw = spath->getConsensus();
    //std::string orf = biostr::stripGap(raw);
    //out << "Length:\t" << "Raw:" << raw.length() << "\tClean:" << orf.length() << "\n";
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
void writeAlignment( std::fstream &out, SpaPath *spath, int count, BitString *bstrs )
{
    out << "ID:" << count << "\n";
    spath->printAlignment( out, 100, bstrs );
    out << "//\n";
}

void writeProfile( std::fstream &out, SpaPath *spath, int count )
{
    out << "ID:" << count << "\n";
    ProfileVector *vprof = spath->getProfileVector();
    //Profile profile;
    //vprof->convert(profile);
    Profile profile = vprof->convert();
    profile.print(out);
    out << "//\n";
    out.flush();
    profile.clean();
}
// void writeProfile( std::fstream &out, MSA &msa, int count )
// {
//     out << "ID:" << count << "\n";
//     Profile profile = msa.getProfile();
//     msa.printProfile( out );
//     out << "//\n";
// }

void writeUnusedKmers( const char *file,
                       DeBruijnGraph &graph,
                       VertexToKmerMap &vertex_map,
                       InvertedIndex &iindex,
                       Param &param)
{
    std::fstream out;
    io::openFile(out, file, std::ios::out);
    Vertex_iter vc, ve;
    for ( boost::tie(vc,ve) = vertices(graph); vc != ve; ++vc ) {
        KmerId kid = vertex_map[*vc];
        int cov = iindex.getValue(kid)->size;
        out << alpha::IntegerToAminoAcid(kid, param.kmer_size) << "\t" << cov << "\n";
    }
    out.close();
}

void writeUnusedReads(PathId *used_reads,
                       int nreads,
                       Param &param
                       )
{
    std::string file = param.out_dir + "/unused_reads.txt";
    std::fstream out;
    io::openFile(out, file.c_str(), std::ios::out);

    ReadIdList unused_reads = assem::getUnusedReads( used_reads, nreads );
    for ( ReadIdList::iterator it = unused_reads.begin(); it != unused_reads.end(); ++it )
        out << *it << "\n";
    out.close();
}



//==============================================================================
// parse program options
// variables for cmd args in header
//==============================================================================
/**
 * Read command line arguments
 */
void readParams(int argc, char **argv, Param &param)
{
    args::parse_cmd_args_assembly( argc, argv, param );
    printLocalTime();
    args::printCommand(argc, argv);
}

void makeGraph( DeBruijnGraph &graph, 
                VertexToKmerMap &vertex_map,
                CoverageMap &kmer_coverage,
                PathToAlnMap &path2aln_map,
                PathId *used_reads,
                int nreads,
                Param &param)
{
    std::cout << "\nBuilding graph ...\n";
    double time1 = mytime();
//     if ( param.file_suffix != "" )
//         loadDump(graph, vertex_map, kmer_coverage, path2aln_map, used_reads, nreads, param);
//     else 
    graph::build(graph, kmer_coverage, vertex_map, param.kmer_size, param.graph_input, param.verbose);

    std::cout << "Graph built:" << mytime()-time1 << " sec\n";
	graph::graphSummary(graph, std::cout);

//     UnDeBruijn ungraph;
//     std::tr1::unordered_map<Vertex, UnVertex> v2v_map;
//     convert(graph, ungraph, v2v_map);

//     computeCC(ungraph);
}

void loadCoverages( std::string &coverage_file,
                    CoverageMap &kmer_coverage,
                    Param &param )
{
    std::fstream fstrm;
    io::openFile(fstrm, coverage_file.c_str(), std::ios::in | std::ios::binary );

    struct stat st;
    int fsize;
    stat(coverage_file.c_str(), &st);
    fsize = st.st_size;
    if ( fsize == 0 ) {
        fstrm.close();
        std::cerr << "[Error] Coverage file is empty\n";
        exit(1);
    }

    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    if ( io::getFileExtension( std::string(coverage_file) ) == "gz" ) 
        in.push(boost::iostreams::gzip_decompressor());
    in.push(fstrm);
    std::istream iin(&in);
                
    CoverageType coverage;
    KmerId kid;
    iin.read((char*)&kid, sizeof(KmerId));
    iin.read((char*)&coverage, sizeof(CoverageType));
    while ( !iin.eof()  )  {
        kmer_coverage[kid] = coverage;
        iin.read((char*)&kid, sizeof(KmerId));
        iin.read((char*)&coverage, sizeof(CoverageType));
    }
    fstrm.close();
    if ( param.verbose ) std::cout << "... k-mer loaded (" << kmer_coverage.size() << ")\n";

}

void loadVertices( std::string &vertex_file,
                   DeBruijnGraph &graph,
                   VertexToKmerMap &vertex_map,
                   std::tr1::unordered_map<void*, void*> &Old2New, 
                   Param &param )
{
    std::fstream fstrm;
    io::openFile(fstrm, vertex_file.c_str(), std::ios::in | std::ios::binary );

    struct stat st;
    int fsize;
    stat(vertex_file.c_str(), &st);
    fsize = st.st_size;
    if ( fsize == 0 ) {
        fstrm.close();
        std::cerr << "[Error] Vertex file is empty\n";
        exit(1);
    }
        
        
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    if ( io::getFileExtension( std::string(vertex_file) ) == "gz" ) 
        in.push(boost::iostreams::gzip_decompressor());
    in.push(fstrm);
    std::istream iin(&in);


    Vertex v;
    void* old;
    KmerId kid;
    iin.read((char*)&old, sizeof(void*));
    iin.read((char*)&kid, sizeof(KmerId));
    while ( !iin.eof() ) {
        v = add_vertex(graph);
        Old2New[old] = v;
        vertex_map[v] = kid;
            
        iin.read((char*)&old, sizeof(void*));
        iin.read((char*)&kid, sizeof(KmerId));
    }
    fstrm.close();
    if ( param.verbose ) std::cout << "... Vertices loaded (" << num_vertices(graph) << ")\n";

}

void loadEdges(std::string &edge_file, 
               DeBruijnGraph &graph,
               std::tr1::unordered_map<void*, void*> &Old2New, 
               Param &param)
{
    std::fstream fstrm;
    io::openFile(fstrm, edge_file.c_str(), std::ios::in | std::ios::binary );

    struct stat st;
    int fsize;
    stat(edge_file.c_str(), &st);
    fsize = st.st_size;
    if ( fsize == 0 ) {
        fstrm.close();
        std::cerr << "[Error] Edge file is empty\n";
        exit(1);
    }
        
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    if ( io::getFileExtension( std::string(edge_file) ) == "gz" ) 
        in.push(boost::iostreams::gzip_decompressor());
    in.push(fstrm);
    std::istream iin(&in);

    void* os;
    void* oe;
    EdgeWeight weight;

    iin.read((char*)&os, sizeof(void*));
    iin.read((char*)&oe, sizeof(void*));
    iin.read((char*)&weight, sizeof(EdgeWeight));
    while ( !iin.eof() ) {
        Edge e = add_edge( Old2New[os], Old2New[oe], graph ).first;
        graph[e].weight = weight;
            
        iin.read((char*)&os, sizeof(void*));
        iin.read((char*)&oe, sizeof(void*));
        iin.read((char*)&weight, sizeof(EdgeWeight));
    }
    fstrm.close();

    if (param.verbose) std::cout << "... Edges loaded (" << num_edges(graph) << ")\n";
}

void loadPaths( std::string &path_file,
                PathToAlnMap &path2aln_map,
                PathId *used_reads,
                Param &param)
{
    std::fstream fstrm;
    io::openFile(fstrm, path_file.c_str(), std::ios::in | std::ios::binary );

    struct stat st;
    int fsize;
    stat(path_file.c_str(), &st);
    fsize = st.st_size;
    if ( fsize == 0 ) {
        fstrm.close();
        if (param.verbose) std::cerr << "\tPath file is empty\n";
    } else {
        boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
        if ( io::getFileExtension( std::string(path_file) ) == "gz" ) 
            in.push(boost::iostreams::gzip_decompressor());
        in.push(fstrm);
        std::istream iin(&in);
        
        PathId pid;
        iin.read((char*)&pid, sizeof(PathId));
        SpaPath *spath = new SpaPath();
        spath->load(iin);
        while ( !iin.eof() ) {
            path2aln_map.insert( std::pair<PathId, SpaPath*>(pid, spath) );
            iin.read((char*)&pid, sizeof(PathId));
            if ( iin.eof() ) break;
            spath = new SpaPath();
            spath->load(iin);
        }
    }
    fstrm.close();
    if ( param.verbose ) std::cout << "... Paths loaded (" << path2aln_map.size() << ")\n";

    for ( PathToAlnMap::iterator it = path2aln_map.begin(); it != path2aln_map.end(); ++it ) {
        PathId pid = it->first;
        SpaPath *spath = it->second;
        ReadId *rid = spath->getReads();
        size_t nread = spath->getReadCount();
        for ( size_t i = 0; i < nread; i++ ) 
            used_reads[rid[i]] = pid;
    }
}

void loadDump( DeBruijnGraph &graph, 
               VertexToKmerMap &vertex_map,
               CoverageMap &kmer_coverage,
               PathToAlnMap &path2aln_map,
               PathId *used_reads,
               int nreads,
               Param &param )
{
    double time1 = mytime();

    param.index_file = param.out_dir + "/iindex." + param.file_suffix;
    
    std::string coverage_file = param.out_dir + "/coverage." + param.file_suffix;
    std::string vertex_file = param.out_dir + "/vertex." + param.file_suffix;
    std::string edge_file = param.out_dir + "/edge." + param.file_suffix;
    std::string path_file = param.out_dir + "/path." + param.file_suffix;

    if ( param.path_flag ) {
        std::tr1::unordered_map<void*, void*> Old2New;
        loadVertices(vertex_file, graph, vertex_map, Old2New, param);
        loadEdges(edge_file, graph, Old2New, param);
        loadCoverages(coverage_file, kmer_coverage, param);
    }
    loadPaths(path_file, path2aln_map, used_reads, param);
        
    if ( param.verbose ) std::cout << "\nDump loaded:" << mytime()-time1 << " sec\n";
}

void loadDebugKmers( std::set<KmerType> &debug_kmers, Param &param )
{
    KmerType kmer;
    std::fstream in;
    io::openFile(in, param.debug_file.c_str(), std::ios::in );
    in >> kmer;
    while ( !in.eof() ) {
        debug_kmers.insert(kmer); 
        in >> kmer;
    }
    in.close();
}


void dumpObjects(DeBruijnGraph &graph, 
               VertexToKmerMap &vertex_map,
               CoverageMap &kmer_coverage,
               InvertedIndex &iindex, 
               PathToAlnMap &path2aln_map,
               Param &param ) 
{
    if ( !param.dump_flag ) return;

    double time1 = mytime();
    if ( param.verbose ) std::cout << "\nDumping objects ...\n";
    if ( param.dump_suffix == "" )
        param.dump_suffix = timeStamp();
    std::string coverage_file = param.out_dir + "/coverage." + param.dump_suffix;
    std::string vertex_file = param.out_dir + "/vertex." + param.dump_suffix;
    std::string edge_file = param.out_dir + "/edge." + param.dump_suffix;
    std::string index_file  = param.out_dir + "/iindex." + param.dump_suffix;
    std::string path_file   = param.out_dir + "/path." + param.dump_suffix;

    // K-mer coverage
    std::fstream out;

    if ( param.path_flag ) {
        io::openFile(out, coverage_file.c_str(), std::ios::out | std::ios::binary );
        for ( CoverageMap::iterator it = kmer_coverage.begin(); it != kmer_coverage.end(); ++it ) {
            out.write((char*)&(*it).first, sizeof(KmerId));
            out.write((char*)&(*it).second, sizeof(CoverageType));
        }
        out.close();
        
        // Vertices
        io::openFile(out, vertex_file.c_str(), std::ios::out | std::ios::binary );
        Vertex_iter vc, ve;
        for ( boost::tie(vc, ve) = vertices(graph); vc != ve; ++vc ) {
            out.write((char*)&(*vc), sizeof(void*));
            out.write((char*)&vertex_map[*vc], sizeof(KmerId));
        }
        out.close();
        
        // Edges
        io::openFile(out, edge_file.c_str(), std::ios::out | std::ios::binary );
        Edge_iter ec, ee;
        for ( tie(ec,ee) = edges(graph); ec != ee; ++ec )  {
            void* s = source(*ec, graph);
            void* t = target(*ec, graph);
            out.write((char*)&s, sizeof(void*));
            out.write((char*)&t, sizeof(void*));
            out.write((char*)&graph[*ec].weight, sizeof(EdgeWeight));
        }
        out.close();
    }

    // kmer to read mapping
    InvertedIndexMap k2r_map = iindex.getInvertedIndex();
    InvertedIndexMap::iterator it;
    io::openFile(out, index_file.c_str(), std::ios::out | std::ios::binary );
    for ( it = k2r_map.begin(); it != k2r_map.end(); ++it ) {
        out.write((char*)&(*it).first, sizeof(KmerId));   // write a kmer ID

        size_t size = (*(it->second)).size;
        out.write((char*)&size, sizeof(unsigned)); // write a count of read IDs
        for ( size_t i = 0; i < size; i++ ) {
            ReadId rid = (*(*it).second).rid[i];
            out.write((char*)&(rid), sizeof(ReadId));
        }
    }
    out.close();

    io::openFile(out, path_file.c_str(), std::ios::out | std::ios::binary);
    for ( PathToAlnMap::iterator pt = path2aln_map.begin(); pt != path2aln_map.end(); ++pt ) {
        out.write((char*)&(pt->first), sizeof(PathId));
        SpaPath *spath = pt->second;
        spath->dump(out);
    }
    out.close();

    if (param.verbose) std::cout << "Objcts written:" << mytime()-time1 << " sec\n";    
}

