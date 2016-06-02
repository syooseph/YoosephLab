#include "postfgs.h"
#include "dna.h"
#include "translation.h"

using namespace std;
using namespace boost;

/**
   Main function of FragGeneScan post processor.

   @return int
   @param argc The number of command line argument.
   @param argv Command line arguments.
   @pre
   -# numberOfPositions contains a positive integer.
   -# firstInitial contains a printable character at least numberOfPositions
   positions from the end of the printable character set.
   @post
   -# firstInitial contains the character numberOfPositions to the right of
   its original value.
   -# stringOfInitials contains a string of length numberOfPositions
   in which each character is the original value of firstInitial.
*/
int main( int argc, char **argv )
{
    double t0 = mytime();

    readParams(argc, argv);

    int nreads = seq::totalSequenceCount( input_files );
    if (verbose) cout << "#seqs:" << nreads << "\n";

    char **seqs = new char *[nreads];
    char **tags = new char *[nreads];
    seq::loadSequences(input_files, tags, seqs, TAGSEQ );

    int size = getSize();
    if (verbose) cout << "Size:" << size << "\n";
    FgsGene *gene = new FgsGene[size];
    loadout(gene, tags, nreads);
    if (verbose) cout << fgs_out.c_str() << " loaded\n";

    char *beg_codons = new char[size]; 
    memset(beg_codons, '0', sizeof(char) * size);
    bool *end_codons = new bool[size];
    memset(end_codons, false, sizeof(bool) * size);
    findStopCodons(gene, size, tags, seqs, end_codons);

    writeNewFaa(beg_codons, end_codons, gene, size, tags);
    
    delete[] gene;
    delete[] beg_codons;
    delete[] end_codons;
    seq::purge(seqs, nreads);
    seq::purge(tags, nreads);

    double t1 = mytime();
    if ( verbose ) cout << "\nElapsed:" << t1-t0 << " sec\n";
    return 0;
}

void parseTag(string &tag,
              int &beg,
              int &end,
              char &strand )
{
    tag.erase(0,1); // delete '>'
    strand = tag[ tag.length()-1 ];
    
    tag.erase( tag.length()-2, 2 ); // delete '_+' or '_-'

    size_t pos = tag.rfind("_");
    if (pos==string::npos) {
        if ( verbose ) cout << "[Error] Sequence tag parsing error\n"; exit(1);
    }
    end = atoi( tag.substr(pos+1, tag.length()-(pos+1)).c_str() );
    tag.erase( pos, tag.length()-pos );

    pos = tag.rfind("_");
    if (pos==string::npos) {
        if ( verbose ) cout << "[Error] Sequence tag parsing error\n"; exit(1);
    }
    beg = atoi( tag.substr(pos+1, tag.length()-(pos+1)).c_str() );
    tag.erase( pos, tag.length()-pos );
}

int correctStart( int start, int frame )
{
    int mod = start % 3;
    if ( verbose ) cout << "old:" << start << "\t";
    start += (frame-mod);
    if ( verbose ) cout << "new:" << start << "\n";
    return start;
}

string getRefSeq( char **seqs,
                  FgsGene *gene,
                  int index,
                  int beg,
                  int end,
                  char strand )
{
    int sid = gene[index].rseq;
    string ref = seqs[sid];

    ref = ref.substr(beg-1, end-beg+1);
    if ( strand == '-' ) ref = dna::reverseComplement(ref);
    return ref;
}

bool shifted( string &seq,
              string &ref,
              int direction )
{
    string query = seq;
    string sbjct = ref;
    if ( direction == RIGHT ) {
        reverse( query.begin(), query.end() );
        reverse( sbjct.begin(), sbjct.end() );
    }

    if ( sbjct.size() >= query.size()+3 ) {
        if ( sbjct.substr(3, query.size()) == query ) return true;
    }
    return false;
}

bool verify( char **tags,
             int index,
             string &tag,
             int beg,
             int end,
             char strand,
             FgsGene *gene )
{
    if ( tag != tags[gene[index].rseq] ) return false; 
    if ( beg != gene[index].spos ) return false;
    if ( end != gene[index].epos ) return false;
    
    char refstr = '+';
    if ( !gene[index].pstr ) refstr = '-';
    if ( strand != refstr ) return false;
    return true;
}

void matchRegion( int &ms,
                  int &me,
                  string &seq,
                  string &ref,
                  int min )
{
    ms = me = -1;
    string sub = seq.substr(seq.size()-min, min);
    size_t found = ref.rfind(sub);
    if (found != string::npos) {
        ms = found;
        me = found+sub.length()-1;
    }
}

void getStops( fstream &log,
               FgsGene *gene,
               char **tags,
               char **seqs,
               std::string  tag, // do not modify
               std::string &seq,
               //char *bcod,
               bool *ecod,
               int index )
{
    char strand;
    int beg, end;
    parseTag(tag, beg, end, strand);
    if ( !verify( tags, index, tag, beg, end, strand, gene ) ) {
        if ( verbose ) cout << "[Error] Tag verification error\n"; exit(1);
    }
    
    string ref = getRefSeq(seqs, gene, index, beg, end, strand);
    
    to_upper( seq );
    to_upper( ref );

    if ( verbose ) {
        std::cout << "\n" << tag << "\tspos:" << beg << "\tepos:" << end << "\tstrand:" << strand << "\tframe:" << gene[index].fram << "\tinsertions:" << gene[index].cins << "\tdeletions:" << gene[index].cdel << "\n";
        std::cout << "Query:" << seq << "\n";
        std::cout << "Sbjct:" << ref << "\n";
    }
    
    bool indel = false;
    if ( gene[index].cins + gene[index].cdel > 0 ) indel = true;

    //=======================================================
    // In case of no indel, whole sequence match
    // With existence of indel, last MIN_MATCH sequence match
    //=======================================================
    size_t min = seq.size();
    if ( indel && min > MIN_MATCH ) min = MIN_MATCH;
    int ms, me;
    matchRegion( ms, me, seq, ref, min );

    if ( verbose ) std::cout << "Match range:" << ms << "\t" << me << "\n";
    if ( ms == -1 || me == -1 ) return;

    if ( me <  (int)seq.size()-1 ) return;

    int i = me+1;
    if ( i+3 > (int)ref.size() ) return;
    
    char codon;
    string tri = ref.substr(i, 3);
    codon = tran::codon.table[tri];
    if ( codon == '*' ) {
        if ( verbose ) std::cout << "True stop codon:" << tri << "\n";
        ecod[index] = true; 
        seq += tri; // append dnas.
        log << tags[index] << "\t" << beg << "\t" << end << "\t" << strand << "\t" << indel << "\t" << ms+1 << "\t" << me+1 << "\t" << tri << "\n";
        return;
    }
    else if ( verbose ) cout << "Not a stop codon:" << tri << "\n";
    
}

// void getTips( fstream &log,
//               FgsGene *gene,
//               char **tags,
//               char **seqs,
//               std::string &tag,
//               std::string &seq,
//               char *bcod,
//               bool *ecod,
//               int index )
// {
//     char strand;
//     int beg, end;
//     parseTag(tag, beg, end, strand);
//     if ( !verify( tags, index, tag, beg, end, strand, gene ) ) {
//         if ( verbose ) cout << "[Error] Tag verification error\n"; exit(1);
//     }
    
//     string ref = getRefSeq(seqs, gene, index, beg, end, strand);
    
//     to_upper( seq );
//     to_upper( ref );

//     if ( verbose ) {
//         std::cout << "SeqId:" << tag << "\n";
//         std::cout << "Start:" << beg << "\tEnd:" << end << "\tStrand:" << strand << "\tFrame:" << gene[index].fram << "\tIns:" << gene[index].cins << "\tDel:" << gene[index].cdel << "\n";
//         std::cout << "Seq:" << seq << "\n";
//         std::cout << "Ref:" << ref << "\n";
//     }

    
//     char bc = '-';
//     char ec = '-';
//     if ( shifted(seq, ref, LEFT) ) {
//         if ( verbose ) {
//             cout << "seq:" << seq << "\n";
//             cout << "ref:" << ref << "\n";
//         }
//         bc = tran::codon.table[ ref.substr(0,3) ];
//         /* Start codons: M(ATG), V(GTG), L(TTG), I(ATT), L(CTG) */
//         //if ( verbose ) cout << "Start AA:" << c << "\n";
//         if ( bc == 'M' || bc == 'V' || bc == 'L' || bc == 'I' )
//             bcod[index] = bc;
//         else if ( verbose ) cout << "Not a start codon:" << bc << "\n";
            
//     }
//     if ( shifted(seq, ref, RIGHT) ) {
//         if ( verbose ) {
//             cout << "seq:" << seq << "\n";
//             cout << "ref:" << ref << "\n";
//         }
//         ec = tran::codon.table[ ref.substr(ref.size()-3,3) ];
//         //if ( verbose ) cout << "End AA:" << c << "\n";
//         if ( ec == '*' ) ecod[index] = true;
//         else if ( verbose ) cout << "Not a stop codon:" << ec << "\n";
//     }
    
//     log << gene[index].rseq << "\t" << beg << "\t" << end << "\t" << strand << "\t" << bc << "\t" << ec << "\n";
    
//}

string dropGzipSuffix( string gzfile )
{
    return gzfile.substr( 0, fgs_out.find_last_of(".") );
}

void findStopCodons( FgsGene *gene,
                     int size,
                     char **tags,
                     char **seqs,
                     //char *bcod,
                     bool *ecod )
{
    //string nffn = fgs_ffn + ".post";
    //fstream out;
    //io::openFile( out, nffn.c_str(), ios::out );
    
    string log_file;
    fgs_out.substr( fgs_out.find_last_of(".") + 1 ) == "gz" ? log_file = dropGzipSuffix(fgs_out) + ".post" : log_file = fgs_out + ".post";
    if ( verbose ) std::cout << "logfile:" << log_file << "\n";

    fstream log;
    io::openFile(log, log_file.c_str(), ios::out);
    log << "#Sequence-ID\tStart\tEnd\tFrame\tIndel\tMatch-start\tMatch-end\tStop-codon\n";
    
    //string ffn_new = fgs_ffn + ".post";
    string ffn_new;
    fgs_ffn.substr( fgs_ffn.find_last_of(".") + 1 ) == "gz" ?  ffn_new = dropGzipSuffix(fgs_ffn) + ".post" : ffn_new = fgs_ffn + ".post";
    if ( verbose ) std::cout << "ffnnew:" << ffn_new << "\n";

    fstream out;
    io::openFile( out, ffn_new.c_str(), ios::out );
    
    string line, tag, seq;
    int curr = 0;
    std::ifstream fstrm(fgs_ffn.c_str(), std::ios_base::in | std::ios_base::binary);
    if ( !fstrm ) {
        std::cout << "Can't open " << fgs_ffn << "\n";
        exit (1);
    }

    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    if ( fgs_ffn.substr( fgs_ffn.find_last_of(".") + 1 ) == "gz" )
        in.push(boost::iostreams::gzip_decompressor());
    in.push(fstrm);
    std::istream istrm(&in);

//     fstream in;
//     io::openFile( in, fgs_ffn.c_str(), ios::in );
//     while ( !in.eof() ) {
    while ( std::getline(istrm, line ) ) {
        //getline(in, line);
        if ( line[0] == '>' ) {
            if (  seq != "" ) {
                //getTips(log, gene, tags, seqs,  tag,  seq, bcod, ecod, curr++ );
                getStops(log, gene, tags, seqs,  tag,  seq, ecod, curr++ );
                out << tag << "\n" << seq << "\n";
            }
            tag = line;  seq = ""; 
        }
        else  seq += line;
    }
    if (  tag != "" &&  seq != "") {
        //getTips(log, gene, tags, seqs,  tag,  seq, bcod, ecod, curr++ );
        getStops(log, gene, tags, seqs,  tag,  seq, ecod, curr++ );
        out << tag << "\n" << seq << "\n";
    }

    //in.close();
    fstrm.close();
    out.close();
//     int count = 0;
//     // Start codons: M(ATG), V(GTG), L(TTG), I(ATT), L(CTG)
//     int cm, cv, cl, ci;
//     cm = cv = cl = ci = 0;
//     for ( int i = 0; i < size; i++ ) {
//         if ( bcod[i] != '0' ) count++;
//         if ( bcod[i] == 'M' ) cm++;
//         if ( bcod[i] == 'V' ) cv++;
//         if ( bcod[i] == 'L' ) cl++;
//         if ( bcod[i] == 'I' ) ci++;
//     }
//     if ( verbose ) std::cout << "start:" << count << "\n";
//     if ( verbose ) std::cout << "M:" << cm << "\tV:" << cv << "\tL:" << cl << "\tI:" << ci << "\n";
    int count = 0;
    for ( int i = 0; i < size; i++ )
        if ( ecod[i] == true ) count++;
    if ( verbose ) std::cout << "stop:" << count << "\n";
}


int getSize()
{
    int size = 0;
    //fstream in;
    string line;
    std::ifstream fstrm(fgs_out.c_str(), std::ios_base::in | std::ios_base::binary);
    if ( !fstrm ) {
        std::cout << "Can't open " << fgs_out << "\n";
        exit (1);
    }
    
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    if ( fgs_out.substr( fgs_out.find_last_of(".") + 1 ) == "gz" )
        in.push(boost::iostreams::gzip_decompressor());
    in.push(fstrm);
    std::istream istrm(&in);

//     io::openFile( in, fgs_out.c_str(), ios::in );
//     getline(in, line);
//     while ( !in.eof() ) {
    while ( std::getline(istrm, line ) ) {
        if ( line[0] != '>' ) size++;
        //getline(in, line);
    }
    //in.close();
    fstrm.close();
    return size;
}

int getCount( string &indel )
{
    int count = 0;
    for ( size_t i = 0; i < indel.length(); i++ ) 
        if ( indel[i] == ',' ) count++;
    return count;
}


void parse( string &line, int &spos, int &epos, int &fram, int &inum, int &dnum, bool &pstr )
{
    sList cols;
    typedef tokenizer<char_separator<char> > tokenizer;
    char_separator<char> sep("\t");
    tokenizer tokens(line, sep);
    for (tokenizer::iterator it = tokens.begin(); it != tokens.end(); ++it) 
        cols.push_back(*it);
    
    int i = 0;
    for ( sList::iterator it = cols.begin(); it != cols.end(); ++it ) {
        ++i;
        if ( i == 1 ) spos = atoi( (*it).c_str() );
        if ( i == 2 ) epos = atoi( (*it).c_str() );        
        if ( i == 3 ) 
            *it == "+" ? pstr = true : pstr = false;
        if ( i == 4 ) fram = atoi( (*it).c_str() );        
        if ( i == 6 ) inum = getCount(*it);        
        if ( i == 7 ) dnum = getCount(*it);
    }
    //cout << spos << "\t" << epos << "\t" << pstr << "\t" << inum << "\t" << dnum << "\n";
}

int findTag( string &tag, char **tags, int nreads, int last )
{
    for ( int i = last; i < nreads; i++ )
        if ( tag == tags[i] ) 
            return i;
    return -1;
}
 
void loadout( FgsGene *gene,
              char **tags,
              int nreads )
{
    int s, e, f, i, d;
    bool p;
    string line;
    int c = 0; // curr gene call
    int r = 0; // read id

    std::ifstream fstrm(fgs_out.c_str(), std::ios_base::in | std::ios_base::binary);
    if ( !fstrm ) {
        std::cout << "Can't open " << fgs_out << "\n";
        exit (1);
    }
    
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    if ( fgs_out.substr( fgs_out.find_last_of(".") + 1 ) == "gz" )
        in.push(boost::iostreams::gzip_decompressor());
    in.push(fstrm);
    std::istream istrm(&in);
    
//     fstream in;
//     io::openFile( in, fgs_out.c_str(), ios::in );
//     getline(in, line);
//     while ( !in.eof() ) {
    while ( std::getline(istrm, line ) ) {
        if ( line[0] == '>' ) {
            line.erase(0,1);
            r = findTag(line, tags, nreads, r );
            if ( r == -1 ) {
                if ( verbose ) cout << "[Erorr] " << line << " not found\n"; exit(1);
            }
        } else {
            parse(line, s, e, f, i, d, p);
            gene[c] = FgsGene(r,s,e, f,i,d,p);
            c++;
        }
        //getline(in, line);
    }
    //in.close();
    fstrm.close();
}

//void loadffn( const char *file )
//{
//
//}
//
//void loadfaa( const char *file )
//{
//
//}

void attach( fstream &out,
             char **tags,
             FgsGene *gene,
             std::string &tag,
             std::string &seq,
             char *bcod,
             bool *ecod,
             int index )
{
    char strand;
    int beg, end;
    parseTag(tag, beg, end, strand);
    if ( !verify( tags, index, tag, beg, end, strand, gene  ) ) {
        if ( verbose ) std::cerr << "[Error] Tag verification error\n"; 
        exit(1);
    }
    
    if ( bcod[index] != '0'  ) seq = bcod[index] + seq;
    if ( ecod[index] == true ) seq += "*";
    
    out << ">" << tag << "_" << beg << "_" << end << "_" << strand << "\n";
    out << seq << "\n";
}
void writeNewFaa(char *bcod, bool *ecod, FgsGene *gene, int size, char **tags)
{
    //string nfaa = fgs_faa + ".post";
    string nfaa;
    fgs_faa.substr( fgs_faa.find_last_of(".") + 1 ) == "gz" ?  nfaa = dropGzipSuffix(fgs_faa) + ".post" : nfaa = fgs_faa + ".post";
    fstream out;
    io::openFile(out, nfaa.c_str(), ios::out);

    int curr = 0;
    string line, tag, seq;

    std::ifstream fstrm(fgs_faa.c_str(), std::ios_base::in | std::ios_base::binary);
    if ( !fstrm ) {
        std::cerr << "Can't open " << fgs_out << "\n";
        exit (1);
    }

    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    if ( fgs_faa.substr( fgs_faa.find_last_of(".") + 1 ) == "gz" )
        in.push(boost::iostreams::gzip_decompressor());
    in.push(fstrm);
    std::istream istrm(&in);


//     fstream fgs;
//     io::openFile(fgs, fgs_faa.c_str(), ios::in);
//     while ( !fgs.eof() ) {
//         getline(fgs, line);
    while ( std::getline(istrm, line ) ) {
        if ( line[0] == '>' ) {
            if (  seq != "" ) 
                attach(out, tags, gene, tag,  seq, bcod, ecod, curr++ );
             tag = line;  seq = ""; 
        }
        else  seq += line;
    }
    if (  tag != "" &&  seq != "") 
        attach(out, tags, gene, tag,  seq, bcod, ecod, curr++ );

    fstrm.close();
}

void readParams(int argc, char **argv)
{
    args::parse_cmd_args_postfgs(argc, argv, input_files, fgs_out, fgs_ffn, fgs_faa, tmp_dir, verbose);
}
