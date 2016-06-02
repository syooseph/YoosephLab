/** 
 * @file validate.cpp
 * @date Thu 2012-01-12 01:14:50 PM
 * @author Youngik Yang
 * @version 0.001
 * @brief Valication of read placement
 * @details 
 * This validates path prediction by calling orf predictor.
 */

#include "rtran.h"


int main( int argc, char* argv[] ) 
{   
    readParams( argc, argv );

    setbuf(stdout, NULL); // no buffering

    //std::cout << "faagz:" << faagz << "\n";
    int norfs;
    //std::cout << faa_file << "\n";
    norfs = seq::getSequenceCount( faa_file.c_str() );//, faagz );
    //norfs = seq::getGzipSequenceCount( faa_file );
    std::cout << "#orfs:" << norfs << "\n";

    char **orfs = new char *[norfs];
    int count = 0;
    seq::loadSequenceFile( faa_file, NULL, orfs, count, SEQONLY );//, faagz );
    std::cout << "ORFs loaded\n";
    assert(count == norfs);

    count = 0;
    char **dnas = new char *[norfs];
    seq::loadSequenceFile( ffn_file, NULL, dnas, count, SEQONLY );//, ffngz );
    std::cout << "DNA read loaded\n";
    assert(count == norfs);

    trimReference( dnas, orfs, norfs );
    std::cout << "DNA read trimmed\n";

    generateDnaSeqs( orfs, dnas );
    //validate(orfs, dnas);

    delete[] orfs;
    delete[] dnas;
    
}

void trimReference( char **dnas, char **orfs, int norfs )
{
    setbuf(stdout, NULL); // no buffering

    int count = 0;
    for ( int i = 0; i < norfs; i++ ) {

        if ( strlen(orfs[i]) == 0 ) {
            std::cout << "[Warning] emtry faa:" << i << "\n";
            continue;
        }

        if ( strlen(dnas[i]) == 0 ) {
            std::cout << "[Warning] emtry ffn:" << i << "\n";
            continue;
        }

        std::string faa = orfs[i];
        std::string ffn = dnas[i];
        std::string orf;
        size_t found = std::string::npos;


        if ( ffn.size() < faa.size()*3 ) {
            std::cout << "[Warning] ffn length:" << ffn.size() << "\tfaa length:" << faa.size() << "\n";
            continue;
        }


        //std::cout << i << "\tfaa:" << faa << "\tffn:" << ffn << "\t";
        //std::cout << ffn.size() << ", " << faa.size() << "\t";
        //std::cout << i << "\tffn:" << ffn.size() << "\tfaa:" << faa.size() << "\t";
        while ( ffn.size() >= faa.size()*3 ) {
            //std::cout << "[" << ffn.size() << "," << faa.size() << " ";
            //std::cout << "\t" << "ffn:" << ffn << "\n";
            orf = tran::translate(ffn);            
            //std::cout << "orf:" << orf << " ";
            if ( orf.size() < faa.size() ) break;
            //std::cout << "\t" << "orf:" << orf << "\n";
            found = orf.find(faa);
            if (found != std::string::npos) break;
            ffn.erase(0, 1);
            if ( ffn.size() < 3 ) break; // not a codon
        }
        //std::cout << "\n";
        //std::cout << i << "orf:" << orf.size() << "\tafter:\tffn:" << ffn.size() << "\tfaa:" << faa.size() << "\n";

        if ( found != std::string::npos ) {     
            //std::cout << "\tFound\n";
            if ( ffn == dnas[i] ) continue;
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

        ++count;
        if ( count % 100000 == 0 ) std::cout << count << " trimming processed\n";
    }
}


void filterGzip( std::istream *istrm,
                 std::fstream &fstrm )
{
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(fstrm);
    std::istream bis(&in);
    istrm = &bis;
}


void generateDnaSeqs( char **orfs,
                      char **dnas )
{
    setbuf(stdout, NULL); // no buffering

    Placement place;
    std::string line;

    std::fstream out;
    io::openFile(out, out_file.c_str(), std::ios::out);

    //std::cout << "Validating\n";
    std::fstream fstrm;
    io::openFile(fstrm, place_file.c_str(), std::ios::in | std::ios::binary);    

    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    //if ( placegz ) 
    if ( io::getFileExtension( place_file ) == "gz" ) 
        in.push(boost::iostreams::gzip_decompressor());
    in.push(fstrm);
    std::istream istrm(&in);


    int count = 0;
    while ( getline( istrm, line ) ) {
        //std::cout << line << "\n";
        if ( line == "" ) continue;
        char c = line[0];
       
        switch(c) {
        case '/' : // placement block ends
            //std::cout << "verification\n";
            rtranslate( place, orfs, dnas, out );
            place.clear();
            ++count;
            if ( count % 100000 == 0 ) std::cout << count << " processed\n";
            break;
        case 'I': // ID
            //std::cout << "id\n";
            //char *num = line.substr(3, line.length()-3).c_str();
            place.id = (unsigned)atoi(line.substr(3, line.length()-3).c_str());
            break;
        case 'R': // Raw sequences
            //std::cout << "raw sequence\n";
            place.seq = biostr::stripGap( line.substr(4, line.length()-4) );
            break;
        case 'C': // Read count
            //std::cout << "count\n";
            break;
        default:
            //std::cout << "adding read\n";
            addRead( line, place );
            break;
        }
    }
//     fstrm->close();
//     delete fstrm;
    fstrm.close();
}

void setAlignPositions( std::map<int, std::list<ReadId> > &aln_pos,
                        Placement &place )
{
    std::vector<ReadId> reads = std::vector<ReadId>( place.reads.begin(), place.reads.end() );
    std::vector<int> starts = std::vector<int>( place.starts.begin(), place.starts.end() );

    for ( size_t i = 0; i < reads.size(); i++ ) {
        aln_pos[starts[i]].push_back( reads[i] );
    }
        //aln_pos.insert( std::pair<int, ReadId>( starts[i], reads[i] ) );
}


// int getRead( std::map<int, ReadId> &aln_pos,
//              std::string &orf,
//              int &lpos,
//              int &rpos )
// {
//     int cpos = rpos;
//     while ( cpos >= lpos ) {
//         int off = rpos-cpos;
//         if ( aln_pos.find(npos) != aln_pos.end() ) {
//             for ( size_t i = 0; i < aln_pos[npos].size(); i++ ) {
//                 std::string pep = orfs[ aln_pos[npos][i] ];
//                 std::string dna = dnas[ aln_pos[npos][i] ];

//                 std::string tra = tran::translate(dna);
//                 if ( tra != pep ) {
//                     std::cout << "Error: two sequences are not identical\n";
//                     std::cout << pep << "\t" << tra << "\n";
//                     exit(-1);
//                 }
//                 //if ( pep[off] == seq[rpos] ) {
//                 if ( pep.substr(0, off+1) == seq.substr(rpos-off, off+1) ) {
//                     lpos = cpos;
//                     return aln_pos[npos][i];
//                 }
//             }
//         }
//         cpos--;
//     }
//     return -1;
// }

std::pair<int, ReadId> searchRead( std::map<int, std::list<ReadId> > &aln_pos,
                                   char **orfs,
                                   char **dnas,
                                   std::string &seq,
                                   int lpos,
                                   int rpos )
{
    int cpos = rpos;
    std::string pep, dna, tra;
    int off;
    while ( cpos >= lpos ) {
        off = rpos-cpos;
        std::cout << "cpos:" << cpos << "\n";
        if ( aln_pos.find(cpos) != aln_pos.end() ) {
            std::vector<ReadId> reads = std::vector<ReadId>( aln_pos[cpos].begin(), aln_pos[cpos].end() );
            
            for ( size_t i = 0; i < reads.size(); i++ ) {
                std::cout << "size:" << reads.size() << "\t" << "i:" << i << "\n";
                pep = orfs[ reads[i] ];
                dna = dnas[ reads[i] ];
                tra = tran::translate(dna);
                std::cout << "sub:" << seq.substr(cpos) << "\n";
                std::cout << "pep:" << pep << "\n";
                std::cout << "dna:" << dna << "\n";
                std::cout << "tra:" << tra << "\n";
                if ( tra != pep ) {
                    std::cerr << "Error: two sequences are not identical\n";
                    std::cerr << pep << "\t" << tra << "\n";
                    exit(-1);
                }
                std::cout << "Pep sub:" << "0\t" << pep.substr(0, off+1) << "\n";
                std::cout << "Seq sub:" << rpos-off << "\t" << seq.substr(rpos-off, off+1) << "\n";
                if ( pep.substr(0, off+1) == seq.substr(cpos, off+1) ) {
                    return std::pair<int, ReadId>( cpos, reads[i] );
                }
            }
        }
        cpos--;
    }
    return std::pair<int, ReadId>( BADPOS, 0);
}

// void reverseTranslate( std::map<int, std::list<ReadId> > &aln_pos,
//                        Placement &place,
//                        char **orfs,
//                        char **dnas )
// {
//     std::string con = place.seq;
//     if ( con == "" ) return;

//     std::string seq = std::string(3*con.size(), '-');
    
//     int lpos = 0;
//     int rpos = 0;
//     //int dpos = 0;
    
//     while ( rpos < (int)con.size() ) {
//         //std::cout << rpos << "\t" << con.size() << "\n";
//         std::cout << "l:" << lpos << "\tr:" << rpos << "\n";
//         std::pair<int, ReadId> srch = searchRead( aln_pos, orfs, dnas, con, lpos, rpos );
//         if ( srch.first == BADPOS ) {
//             std::cout << "Bad pos\n";
//             lpos = rpos; rpos++;
//         }
//         else {
//             std::string pep = orfs[ srch.second ];
//             std::string dna = dnas[ srch.second ];
//             int match = 0;            
//             for ( size_t i = 0; i < pep.size(); i++ ) {
//                 if ( pep[i] == con[srch.first+i] ) match++;
//             }
//             std::cout << "#match:" << match << "\n";
//             std::cout << con.substr(srch.first, match) << "\n";
//             std::cout << "dna:" << dna.substr(0, 3*match) << "\n";
//             seq.replace( 3*srch.first, 3*match, dna.substr(0, 3*match) );
//             std::cout << "s:" << 3*srch.first << "\tn:" << 3*match << "\n";
//             std::cout << "so far:" << seq << "\n";
//             lpos = srch.first;
//             rpos = srch.first + match;
//         }
//     }

//     std::cout << "ORF:" << con << "\n";
//     std::cout << "DNA:" << seq << "\n\n";

//     std::cout << ">" << place.id << "\n";
//     std::cout << seq << "\n";
// }

void setPaddedReads( std::vector<std::string> &DnaReads,
                     Placement &place,
                     char **dnas )
{
    std::vector<ReadId> rids = std::vector<ReadId>( place.reads.begin(), place.reads.end() );
    std::vector<int>  starts = std::vector<int>( place.starts.begin(), place.starts.end() );
    for ( size_t i = 0; i < rids.size(); i++ ) {
        std::string padded = std::string( place.seq.size()*3, 'X' );
        std::string nread = dnas[rids[i]];
        //std::cout << rids[i] << "\t" << nread << "\n";
        std::vector<int> dels, inss;
        if ( place.dels.find(rids[i]) != place.dels.end() ) dels = place.dels[rids[i]];
        if ( place.inss.find(rids[i]) != place.inss.end() ) inss = place.inss[rids[i]];
        

        std::multimap<int, int> indels;
        for ( int j = 0; j < (int)inss.size(); j++ ) 
            indels.insert(std::pair<int,int>(inss[j], 1));
        for ( int j = 0; j < (int)dels.size(); j++ ) 
            indels.insert(std::pair<int,int>(dels[j], 0));

        std::multimap<int, int>::reverse_iterator it;
        for ( it = indels.rbegin(); it != indels.rend(); ++it ) {
            assert( (int)nread.size() > it->first*3 );
            if ( it->second == 1 ) nread.erase( it->first*3, 3 );      // insertion -> drop 3 DNAs
            if ( it->second == 0 ) nread.insert( it->first*3, "NNN" ); // deletion -> insert 3 Ns
        }


//         //if ( starts[i] < 0 ) {
//         //    nread.erase( 0, -1*starts[i]*3 );
//         for ( int j = (int)inss.size()-1; j >=0 ; j-- ) {
//             if ( nread.size() <= inss[j]*3 ) {
//                 std::cout << "Error:" << rids[i] << "\tsize:" << nread.size() << "\tins:" << inss[j] << "\n";
//             }
//             assert( (int)nread.size() > inss[j]*3 );
//             //nread.erase( (inss[j]+starts[i])*3, 3);
//             nread.erase( inss[j]*3, 3 );
//         }
//         for ( int j = (int)dels.size()-1; j >=0 ; j-- ) {
//             if ( nread.size() <= dels[j]*3 )
//                 std::cout << "Error:" << rids[i] << "\tsize:" << nread.size() << "\tdel:" << dels[j] << "\n";
//             assert( (int)nread.size() > dels[j]*3 );
//             //nread.insert( (dels[j]+starts[i])*3, "NNN");
//             nread.insert( dels[j]*3, "NNN" );
//         }
//         //std::cout << rids[i] << "\t" << nread << "\n";

        //starts[i] = 0;
        //} else {
//         for ( int j = (int)inss.size()-1; j >=0 ; j-- ) {
//             assert( (int)nread.size() > inss[j]*3 );
//             nread.erase( inss[j]*3, 3 );
//         }
//             for ( int j = (int)dels.size()-1; j >=0 ; j-- ) {
//                 assert( (int)nread.size() > dels[j]*3 );
//                 nread.insert( dels[j]*3, "NNN");
//             }
//         }
        if ( starts[i] < 0 ) {
            nread.erase( 0, -1*starts[i]*3 );
            starts[i] = 0;
            //std::cout << rids[i] << "\t" << nread << "\n";
        }

        assert( starts[i]*3 < (int)padded.size() );
        size_t nbase = nread.size();
        if ( starts[i]*3 + nread.size() > padded.size() ) 
            nbase = padded.size() - starts[i]*3;
        padded.replace( starts[i]*3, nbase, nread, 0, nbase );

        DnaReads.push_back(padded);
        //std::cout << rids[i] << "\t" << padded << "\n";
    }
}

void setProfile( Profile &profile,
                 std::vector<std::string> &DnaReads )
{
    for ( size_t i = 0; i < profile.ncol; i++ ) 
        profile.matrix[i] = new unsigned[NDNA];
    
    for ( size_t i = 0; i < profile.ncol; i++ ) {
        for ( size_t j = 0; j < NDNA; j++ ) {
            profile.matrix[i][j] = 0;
        }
    }

    for ( size_t i = 0; i < DnaReads.size(); i++ ) {
        for ( size_t j = 0; j < profile.ncol; j++ ) {
            char dna = DnaReads[i][j];
            if ( dna == 'X' ) continue;
            
            switch(dna) {
            case 'A':
                profile.matrix[j][0]++;
                break;
            case 'C':
                profile.matrix[j][1]++;
                break;
            case 'G':
                profile.matrix[j][2]++;
                break;
            case 'T':
                profile.matrix[j][3]++;
                break;
            default :
                profile.matrix[j][4]++;
                break;
            }
        }
    }
}

std::string getConsensus(Profile &profile)
{
    std::string consensus = std::string(profile.ncol, '-');
    for ( size_t i = 0; i < profile.ncol; i++ ) {
        std::multimap<int, int> cmap;
        for ( size_t j = 0; j < NDNA; j++ ) 
            cmap.insert( std::pair<int, int>(profile.matrix[i][j], j) );

        std::multimap<int, int>::reverse_iterator it = cmap.rbegin();
        consensus[i] = DNA[it->second];
    }
    return consensus;
}

void rtranslate( Placement &place,
                 char **orfs,
                 char **dnas,
                 std::fstream &out )
{
    std::vector<std::string> DnaReads;
    setPaddedReads( DnaReads, place, dnas );

    Profile profile;
    profile.ncol = place.seq.size()*3;
    profile.matrix = new unsigned*[profile.ncol];
    setProfile(profile, DnaReads);

    std::string consensus = getConsensus(profile);    
    std::string ntra = tran::translate(consensus);
 
    std::cout << "ID:" << place.id << "\n";
    std::cout << "ORF:" << place.seq << "\n";
    std::cout << "DNA:" << consensus << "\n";
    std::cout << "NEW:" << ntra << "\n";

    double positive = 0;
    double identity = 0;
    if ( ntra == place.seq ) {
        positive = identity = 100;
        //std::cout << "[Warning] ORF != NEW\n";
        //std::cout << "Identity:"
    } else {
        int cmat = 0;
        int cpos = 0;
        for ( size_t i = 0; i < place.seq.size(); i++ ) {
            if ( place.seq[i] == ntra[i] ) cmat++;
            if ( scoring::getScore(place.seq[i], ntra[i], BLOSUM62) > 0 ) cpos++;
        }
        identity = (double)cmat/place.seq.size() * 100;
        positive = (double)cpos/place.seq.size() * 100;
    }

    char buff[100];
    sprintf(buff, "%.2f", identity);
    std::cout << "identity:" << buff << "\n";
    sprintf(buff, "%.2f", positive);
    std::cout << "positive:" << buff << "\n";
    std::cout << "//\n";

    out << ">" << place.id << "\n" << consensus << "\n";
    profile.clean();
}

void addRead( std::string &line,
              Placement &place )
{
    using namespace boost;

    std::list<std::string> cols;
    typedef tokenizer<char_separator<char> > tokenizer;
    char_separator<char> sep("\t");
    tokenizer tokens(line, sep);
    for (tokenizer::iterator it = tokens.begin(); it != tokens.end(); ++it)
        cols.push_back(*it);

    ReadId rid; 
    int i = 0;
    for ( std::list<std::string>::iterator it = cols.begin(); it != cols.end(); ++it ) {
        ++i;
        if ( i == 2 ) {
            rid = (ReadId) atoi( (*it).c_str() );
            //std::cout << rid << "\n";
            place.reads.push_back( rid );
        }
        if ( i == 3 ) place.starts.push_back( atoi( (*it).c_str() ) );
        if ( i == 4 ) { 
            std::vector<int> ins = getIndel(*it);
            if ( ins.size() ) place.inss.insert( std::pair<ReadId, std::vector<int> >( rid, ins ) );
        }
        if ( i == 5 ) {
            std::vector<int> del = getIndel(*it);
            if ( del.size() ) place.dels.insert( std::pair<ReadId, std::vector<int> >( rid, del ) );
        }
    }
}

std::vector<int> getIndel( std::string &col )
{
    using namespace boost;

    std::vector<int> pos;
    col.erase(0, 2);
    if ( col == "" ) return pos;
 

    std::list<std::string> locs;
    typedef tokenizer<char_separator<char> > tokenizer;
    char_separator<char> sep(",");
    tokenizer tokens(col, sep);
    for (tokenizer::iterator it = tokens.begin(); it != tokens.end(); ++it)
        locs.push_back(*it);
    
    for ( std::list<std::string>::iterator it = locs.begin(); it != locs.end(); ++it ) {
        if ( *it == "" )  continue;
        pos.push_back( atoi( (*it).c_str() ) );
    }
    //for ( size_t i = 0; i < pos.size(); i++ )
    //std::cout << pos[i] << "\n";
    return pos;
}

void readParams(int argc, char **argv)
{
    args::printCommand(argc, argv);
    
    args::parse_cmd_args_postspa( argc, 
                                  argv, 
                                  faa_file,
                                  ffn_file,
                                  place_file,
                                  out_file );

}
