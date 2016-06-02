#include "loader.h"

Loader::Loader() 
{
    nreads  = 0;
    gsa     = NULL;
    seqs    = NULL;
    tags    = NULL;
    strands = NULL;
    pairs   = NULL;

    status = false;
}

Loader::~Loader()
{
    clear();
}

void Loader::purgeGSA()
{
    if ( gsa != NULL ) delete[] gsa;
    gsa = NULL;
}

void Loader::clear()
{
    if ( gsa     != NULL ) delete[] gsa;
    if ( seqs    != NULL ) seq::purge( seqs, nreads );
    if ( strands != NULL ) delete[] strands;
    if ( pairs   != NULL ) delete[] pairs;
    gsa     = NULL;
    seqs    = NULL;
    tags    = NULL;
    strands = NULL;
    pairs   = NULL;
    nreads  = 0;
    status  = false;
    if ( Param::verbose ) std::cout << "Loader cleaned\n";
}

void Loader::load()
{
    nparts = Param::nparts;

    loadReads();
    getLengths();
    setStrands(tags, nreads);
    setReadPairs(tags, nreads, Param::pair_flag);

    trimReads();

    //-----------------
    // No further needs
    //-----------------
    if (tags!=NULL)seq::purge(tags, nreads);

    loadSuffixArray();

    status = true;
}

void Loader::loadReads()
{
    nreads = seq::totalSequenceCount( Param::input_files );
    if ( Param::verbose ) std::cout << "Total sequences:" << nreads << "\n";

    //-------------------------------------------
    // Load sequences and tags from input file(s)
    //-------------------------------------------
    tags = new char *[nreads];
    seqs = new char *[nreads];
    seq::loadSequences(Param::input_files, tags, seqs, TAGSEQ);
    printElapsed( INIT_TIME, mytime(), "Sequence loaded" );
}

void Loader::getLengths()
{
    min_len = 1000000;
    max_len = 0;
    for ( int i = 0; i < nreads; i++ ) {
        int len = strlen(seqs[i]);
        if ( len > max_len ) max_len = len;
        if ( len < min_len ) min_len = len;
    }
}

void Loader::trimReads()
{
    std::string bad_file = Param::out_dir + "/bad";
    if ( ! fio::existFile(bad_file.c_str()) ) return;
    
    std::fstream in;
    fio::openFile(in, bad_file.c_str(), std::ios::in | std::ios::binary );
    size_t nbad;
    in.read((char*)&nbad, sizeof(size_t));
    bad_reads.reserve(nbad);
    int rid;
    for ( size_t i = 0; i < nbad; i++ ) {
        in.read((char*)&rid, sizeof(int));
        assert(rid<nreads);
        bad_reads.push_back(rid);
    }
    printElapsed( INIT_TIME, mytime(), "Invalid sequences loaded" );
}

void Loader::loadSuffixArray()
{
    gsa = new GSA[Param::nparts];

    std::string sfa_file, doc_file, lcp_file, mcp_file, gsa_file;
    for ( int i = 0; i < Param::nparts; i++ ) {
        // lcp_file = Param::out_dir + "/lcp." + boost::lexical_cast<std::string>(i);
        // mcp_file = Param::out_dir + "/mcp." + boost::lexical_cast<std::string>(i);
        // gsa_file = Param::out_dir + "/gsa." + boost::lexical_cast<std::string>(i);
        lcp_file = Param::out_dir + "/lcp." + std::to_string(i);
        mcp_file = Param::out_dir + "/mcp." + std::to_string(i);
        gsa_file = Param::out_dir + "/gsa." + std::to_string(i);

        gsa[i].load( lcp_file.c_str(),
                     mcp_file.c_str(),
                     gsa_file.c_str());
    }
    printElapsed( INIT_TIME, mytime(), "Suffix array loaded" );
}

void Loader::setStrands(char **tags, int nreads)
{
    strands = new char[nreads];
    for ( int i = 0; i < nreads; i++ ) {
        std::string rid = tags[i];
        char s = rid[rid.size()-1];
        if ( s == '+' || s == '-' ) strands[i] = s;
        else strands[i] = '\0';
    }
}

void Loader::setReadPairs(char **tags, int nreads, bool pair_flag)
{
    if ( !pair_flag ) {
        pairs = NULL; return;
    }

    pairs = new ReadId[nreads];

    char **ntag = new char*[nreads];
    for (int i = 0; i < nreads; i++ ) {
        int len = strlen(tags[i]);
        ntag[i] = new char[len+1];
        strncpy( ntag[i], tags[i], len+1 );
    }
        
    dropPairInfo(ntag, nreads);

    //-------------------------------------------
    // initialize each with total number of reads
    //-------------------------------------------
    for ( int i = 0; i < nreads; i++ ) pairs[i] = NOT_PAIR;
    
    for ( int i = 0; i < nreads-1; i++ ) {
        if ( pairs[i] != NOT_PAIR ) continue;
        
        if ( strcmp(ntag[i], ntag[i+1]) == 0 ) {
            pairs[i] = i+1;
            pairs[i+1] = i;
        }
    }

    for ( int i = 0; i < nreads; i++ )
        delete[] ntag[i];
    delete[] ntag;
}

void Loader::dropPairInfo(char **tags, int nreads)
{
    for ( int i = 0; i < nreads; i++ ) {
        std::string rid = tags[i];
        size_t p = rid.find("/");
        if ( p >= rid.size() ) {
            std::cerr << "[Error] Invalid paired end tag\n"; 
            exit(EXIT_FAILURE);
        }
        delete[] tags[i];

        std::string nrid = rid.substr(0, p);
        int nlen = nrid.size()+1;
        tags[i] = new char[nlen];
        strncpy( tags[i], nrid.c_str(), nlen );
    }
}
