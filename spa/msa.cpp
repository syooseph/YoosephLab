#include "msa.h"

/**
 * Default constructor
 */
MSA::MSA()
{
    init();
}

/** 
 * Constructor 
 */
MSA::MSA(SpaPath &spath, PathId pid, PathId *path_reads, BitString *bstrs, char *strands, ReadId *pairs, int mode, int &error, Param &param )
{
    init();
	build(spath, pid, path_reads, bstrs, strands, pairs, mode, error, param);
}

/** 
 * Destructor 
 */
MSA::~MSA()
{
    profile.clean();
}

void MSA::init()
{
    profile.ncol = 0;
    aligned = false;
}

/**
 * @return consensus string
 */
std::string MSA::getConsensus()
{
    //return makeConsensus();
    return consensus;
}

/**
 * @return alignment profile
 */
Profile MSA::getProfile()
{
    return profile;
}

/**
 * Consensus module.
 * Based on profile, report single most frequet AA in each column.
 * If there is a tie involing stop codon, use stop as a consensus.
 */
void MSA::makeConsensus()
{
    consensus = std::string(profile.ncol, '.');
    for ( size_t i = 0; i < profile.ncol; i++ ) {
        std::multimap<int, int> cmap;
        for ( size_t j = 0; j < (size_t)naas; j++ ) 
            cmap.insert( std::pair<int, int>(profile.matrix[i][j], j) );

        std::multimap<int, int>::reverse_iterator it = cmap.rbegin();
        int max = it->first;
        if ( max == 0 ) continue;

        consensus[i] = alpha::AminoAcid[(it->second)+1];
        if ( consensus[i] != '*' ) {
            ++it;
            for ( ; it != cmap.rend(); ++it ) {
                if ( it->first < max ) break;
                if ( alpha::AminoAcid[(it->second)+1] == '*' ) {
                    consensus[i] = alpha::AminoAcid[(it->second)+1];
                    break;
                }
            }
        }
    }
}



//====================================================================
// Generation of 2D profile matrix
// Row: reads
// Col: AAs
//====================================================================
void MSA::makeProfile()
{
    if ( profile.ncol ) profile.clean();
    
    profile.init(consensus.size());
    profile.update(nreads);

//     __initializeProfile();
//    __updateProfile();
}

//====================================================================
// 2D profile matrix filled with 0 count
//====================================================================
void MSA::__initializeProfile()
{
    profile.ncol = consensus.size();
    profile.matrix = new unsigned*[profile.ncol];
    for ( size_t i = 0; i < profile.ncol; i++ ) 
        profile.matrix[i] = new unsigned[naas];
}

//====================================================================
// Count filling in profile matrix
//====================================================================
void MSA::__updateProfile()
{
    __resetProfile();

    for ( size_t i = 0; i < nreads.size(); i++ ) {
        for ( size_t j = 0; j < profile.ncol; j++ ) {
            char aa = nreads[i][j];
            if ( aa == '.' ) continue;
            int na = alpha::AsciiToAA[(int)aa];
            if ( na < 1 || na > 27 ) 
                std::cout << "Invalid aa:" << aa << "\tnum:" << na << "\tcol:" << j << "\t" << nreads[i] << "\n";
            profile.matrix[j][na-1]++;
        }
    }
}

void MSA::__resetProfile()
{
    for ( size_t i = 0; i < profile.ncol; i++ ) {
        for ( size_t j = 0; j < (size_t)naas; j++ ) {
            profile.matrix[i][j] = 0;
        }
    }
}



//====================================================================
// Profile display
//====================================================================
void MSA::printProfile(std::ostream &out)
{
    char buf[25];

    //------------
    // amino acids
    //------------
    out << "\t";    
    for ( int i = 0; i < naas; i++ ) {
        sprintf(buf, "%3c ", alpha::AminoAcid[i+1]);
        out << buf;
    }
    out << "\n";

    //-------
    // matrix
    //------- 
    for ( size_t i = 0; i < profile.ncol; i++ ) {
        out << i << "\t";
        for ( int j = 0; j < naas; j++ ) {
            sprintf(buf, "%3d ", profile.matrix[i][j]);
            out << buf;
        }
        out << "\n";
    }
}

Task MSA::parseMode( int mode )
{
    Task flag;
    if ( mode >= CODON ) {
        flag.codon = true;
        mode -= CODON;
    } 
    if ( mode >= TRIM ) {
        flag.trim = true;
        mode -= TRIM;
    }
    if ( mode >= FLAT ) {
        flag.flat = true;
        mode -= FLAT;
    }
    if ( mode >= PILE ) {
        flag.pile = true;
    }
    return flag;
}

void MSA::__update(SpaPath &spath, int kmer_size, bool verbose)
{
    makeProfile();
    updateConsensus(spath, kmer_size); 
    updateAllPosFromConsensus(spath, verbose);
    //consensus = biostr::stripGap(consensus);
    aligned = true;
}

void MSA::build(SpaPath &spath, PathId pid, PathId *path_reads, BitString *bstrs, char *strands, ReadId *pairs, int mode, int &error, Param &param )
{
    Task task = parseMode( mode );

    /* reset path id membership */
    joinUsedReads( NOT_PATH, spath.getReads(), spath.getReadCount(), path_reads );

    //=====================================================
    // Multiple sequence alignment with temporary consensus
    //=====================================================
    if ( task.pile ) {
        if ( !pileUpReads(spath, bstrs, param.kmer_size, param.verbose) ) {
            error = PILE; return; 
        }
    }
    //=================================
    // Trim end bases with low coverage
    //=================================
    if ( task.trim ) {
        while (true) {
            size_t csize = consensus.size();
            if ( csize == 0 ) { error = TRIM; return; }

            trimPartialMatches( spath, bstrs, param.merge_score, param.verbose );
            if ( spath.getReadCount() == 0 ) {
                error = TRIM; return; 
            }
            FlagPair flags = trimEnds(spath, param.base_depth, param.verbose);
            if ( flags.second ) { // error occurred
                error = TRIM; return; 
            }
            if ( flags.first ) break; // no need to trim
            if (!aligned) __update(spath, param.kmer_size, param.verbose);
            //if ( csize == consensus.size() ) break; // no change made
        }
    }

    //=======================
    // Fix zero coverage hole
    //=======================
    if ( task.flat ) {
		if ( ! __flatten(bstrs, strands, pairs, spath, param ) ) {
            error = FLAT; return; 
        }
    }

    //=============================================
    // Trim sequence after stop condon in consensus
    //=============================================
    if ( task.codon ) {
        if ( !aligned ) __update(spath, param.kmer_size, param.verbose);
        while ( true ) {
            size_t csize = consensus.size();
            if ( csize == 0 ) { error = CODON; return; }
            
            if ( !trimAfterStopCodon( spath, param.verbose ) ) break;

            //updateConsensus(spath, param.kmer_size); 
            //__update(spath, param.kmer_size);
            trimPartialMatches( spath, bstrs, param.merge_score, param.verbose );
            if ( spath.getReadCount() == 0 ) {
                error = TRIM; return; 
            }
            //bool flag = trimEnds(spath, param.base_depth, param.verbose);
            FlagPair flags = trimEnds(spath, param.base_depth, param.verbose);
            if ( flags.second ) {
                error = CODON; return; 
            }
            if ( flags.first ) break;
            
            //updateConsensus(spath, param.kmer_size); 
            if (!aligned) __update(spath, param.kmer_size, param.verbose);
            if ( csize == consensus.size() ) {
                size_t pos = consensus.find("*");
                if ( pos == std::string::npos || pos == consensus.size()-1 ) 
                    break; // no change made
            }
        }
    }
     
    ////updateConsensus(spath, param.kmer_size);    

    if (!aligned) __update(spath, param.kmer_size, param.verbose);

    spath.resetProfile();
    spath.setProfile(profile);

    //updateAllPosFromConsensus(spath);

    joinUsedReads( pid, spath.getReads(), spath.getReadCount(), path_reads );

}

bool MSA::trimPartialMatches( SpaPath &spath, BitString *bstrs, double merge_score, bool verbose )
{
    if ( verbose ) std::cout << "# reads before trimming partial matches:" << spath.getReadCount() << "\n";
    std::list<int> partials = spath.getPartialMatches( bstrs, merge_score, verbose);
    if ( verbose && partials.size() ) 
        std::cout << "# Partial aligned reads:" << partials.size() << "\n";

    if ( partials.size() == 0 ) return false;

    dropReads(partials);
    spath.dropReads(partials);
    if ( verbose ) std::cout << "# reads after trimming partial matches:" << spath.getReadCount() << "\n";
    return true;
}


int MSA::maxEndPosition( SpaPath &spath, BitString *bstrs, int k, bool verbose)
{
    int maxe = (int)consensus.size()-1;
    size_t nr  = spath.getReadCount();
    ReadId* rids = spath.getReads();
    int *inits = spath.getInits();
    for ( size_t i = 0; i < nr; i++ ) {
        std::string str = bstrs[rids[i]].toString();
        int epos = inits[i] + str.size() - 1;
        if ( epos > maxe ) {
            maxe = epos;
        }
    }
    return maxe;
}

int MSA::minStartPosition( SpaPath  &spath )
{
    int mins = 0;
    size_t nr  = spath.getReadCount();
    int *inits = spath.getInits();
    for ( size_t i = 0; i < nr; i++ )
        if ( inits[i] < mins ) 
            mins = inits[i];
    return mins;
}


void MSA::adjustStartPositions( SpaPath &spath, int k, bool verbose )
{
    /* Find minimum start position alignment */
    int mins = minStartPosition( spath );
    if ( mins >= 0 ) return;

    if ( verbose ) std::cout << "** Minimum start:" << mins << "\n";
    spath.adjustPositions(-1*mins);
    std::string ostr = biostr::getSequenceString( spath.getKmers(), spath.getKmerCount(), k );
    std::string nstr = std::string( -1*mins, 'X' ) + ostr;
    if ( verbose ) {
        std::cout << "ostr:" << ostr << "\n";
        std::cout << "nstr:" << nstr << "\n";
    }
    //KmerArray kmers = biostr::getKmers( biostr::stripGap(nstr), k );
    KmerArray kmers = biostr::getKmers( nstr, k );
    spath.updateKmers( &kmers[0], kmers.size() );
}

bool MSA::pileUpReads( SpaPath &spath, BitString *bstrs, int k, bool verbose)
{
    if ( spath.getReadCount() == 0 ) return false;   


    /* In case of negative start position */
    adjustStartPositions( spath, k, verbose );

    /* Find deletion locations in each read and insert gaps */
    insertGapsToReads(spath, bstrs, verbose);

    //---------------------------------------------------------------------------
    // Find insertion locations in each read and corresponding reference location 
    // Key: read-id, Value(seq-pos, ref-pos)
    //---------------------------------------------------------------------------
    ReadPosPairListMap pmap = loadInsertions(spath, verbose);
    
    //---------------------------------------------------
    // Count no of gaps of reads in reference position
    // Key: reference position
    // Value: count map - key:ReadId, value:#gaps
    //---------------------------------------------------
    PosReadCountsMap pcmap = countReadPositions(pmap);

    //----------------------------------------------------------------
    // Find count of maximum insertion size in each sbjct gap position
    // Key: reference position
    // Value: maximum gaps
    //---------------------------------------------------------------
    PosCountsMap cmap = getMaxCounts(pcmap);

    /* Insert gaps to reference sequence */
    insertGapsToRef(cmap, spath, k, verbose);

    /* Stretch consensus sequence and reads */
    stretch(spath, bstrs, k, verbose);

    /* Insert gaps to reads due to sbjct insertions */
    adjustInsertions(spath, pmap, pcmap, cmap, verbose);

    fitLength(verbose);

    /* Update profile, consensus, and positions */

    __update(spath, k, verbose);

//     makeProfile();
//     makeConsensus();

//     updateAllPosFromConsensus(spath);

//     //     std::cout << "After pile up alignment\n";
//     //     printAlignment(std::cout, spath, 100);

//     //=======================================================================
//     //  Do not make consensus yet, which makes zero coverage regions all gaps
//     //======================================================================= 
    return true;
}

// int countLeadGap(int s)
// {
//     int lgap = 0;
//     for ( size_t i = 0; i < s; i++ )
//         if ( consensus[i] == '-') lgap++;
//     return lgap;
// }

void MSA::updateInit( SpaPath &spath, bool verbose )
{
    ReadId* rids = spath.getReads();
    int *inits   = spath.getInits();
    
    for ( size_t i = 0; i < nreads.size(); i++ ) {
        int sbjct_gaps = 0;
        size_t j = 0;
        while ( nreads[i][j] == '.' ) {
            if ( j == nreads[i].size() ) break;
            if ( consensus[j] == '-' ) sbjct_gaps++;
            j++;
        }
        
        if ( j == nreads[i].size() ) {
            if (verbose) {
                std::cout << "Something wrong\n";
                std::cout << rids[i] << "\tinit:" << inits[i] << "\tj:" << j << "\t" << nreads[i] << "\n";
            }
            continue;
        }
        
        
        if ( j > 0 ) {
            //if ( verbose ) 
            //std::cout << "Start\t" << rids[i] << "\t" << inits[i] << " -> " << j-sbjct_gaps << "\n";
            spath.setInit(i, j-sbjct_gaps);
        }
    }
}

// void MSA::updateAllPosFromConsensus( SpaPath &spath, bool verbose )
// {
//     if ( verbose ) std::cout << "Updating all positions\n";

//     ReadId* rids = spath.getReads();
//     int *inits   = spath.getInits();

//     spath.resetMismatches();
//     std::list<Mismatch>ilist, dlist;

//     updateInit(spath);

//     size_t c = 0;
//     for ( ; ; ) {
//         if ( c >= consensus.size() ) break;
//         if ( consensus.size() == 0 ) break;
//         /* Insertion handling */
//         if ( consensus[c] == '-' ) {
//             for ( size_t i = 0; i < nreads.size(); i++ ) {
//                 if ( nreads[i][c] != '-' && nreads[i][c] != '.' ) {
//                     ilist.push_back( Mismatch( rids[i], c-inits[i], c) );
//                     std::cout << rids[i] << "\tinsertion\t" << c-inits[i] << "\t" << c << "\n";
//                 }
//                 nreads[i].erase(c, 1);
//             }                                   
//             consensus.erase(c,1);
//         }
//         else {
//             for ( size_t i = 0; i < nreads.size(); i++ ) {
//                 if ( nreads[i][c] == '-' ) {
//                     dlist.push_back( Mismatch( rids[i], c-inits[i], c) );
//                     std::cout << rids[i] << "\tdeletion\t" << c-inits[i] << "\t" << c << "\n";
//                 }
//             }
//             c++;
//         }
//     }
//     spath.setMismatches(ilist, INSERTION);
//     spath.setMismatches(dlist, DELETION);
// }


void MSA::updateAllPosFromConsensus( SpaPath &spath, bool verbose )
{
    if ( verbose ) std::cout << "Updating all positions\n";

    ReadId* rids = spath.getReads();
    int *inits   = spath.getInits();

    spath.resetMismatches();

    std::list<Mismatch>ilist, dlist;

    for ( size_t i = 0; i < nreads.size(); i++ ) {
        int sbjct_gaps = 0;
        int query_gaps = 0;
        size_t j = 0;
        while ( nreads[i][j] == '.' ) {
            if ( j == nreads[i].size() ) break;
            if ( consensus[j] == '-' ) sbjct_gaps++;
            j++;
        }

        if ( j == nreads[i].size() ) {
            if (verbose) {
                std::cout << "Something wrong\n";
                std::cout << rids[i] << "\tinit:" << inits[i] << "\tj:" << j << "\t" << nreads[i] << "\n";
            }
            continue;
        }

        //spath.setInit(i, j); // update alignment start position
        
        //if ( inits[i] > 0 && sbjct_gaps ) {
        if ( j > 0 ) {
//             if ( inits[i] != int(j-sbjct_gaps) ) {
//                 if ( verbose ) std::cout << "New start\t" << rids[i] << "\t" << inits[i] << " -> " << j-sbjct_gaps << "\n";
//             }
            //if ( verbose ) std::cout << "Start\t" << rids[i] << "\t" << inits[i] << " -> " << j-sbjct_gaps << "\n";
            spath.setInit(i, j-sbjct_gaps);
        }
        
        int lead_gaps = sbjct_gaps;
        while ( nreads[i][j] != '.' ) {
            if ( j == nreads[i].size() ) break;

            //     .    :    .    :
            // ABCD--GHIJK-MN--QRST
            // ---------J---NOP-R
            //           123 456
            // init = 7
            // 1. 10 - 7 - 2 = 1      del
            // 2. skip
            // 3. 12 - 7 - 2 - 2 = 1  del
            // 4. 14 - 7 - 2 - 3 = 2  ins
            // 5. 15 - 7 - 2 - 3 = 3  ins
            // 6. 16 - 7 - 2 - 3 = 4  ins

            if ( consensus[j] != '-' ) {
                if ( nreads[i][j] == '-' ) {
                    //inits[i] >= 0 ?
                    //dlist.push_back( Mismatch( rids[i], j-query_gaps, j-sbjct_gaps) );
                    if ( verbose ) std::cout << rids[i] << "\tdeletion\t" << j-inits[i]-lead_gaps-query_gaps << "\t" << j-sbjct_gaps << "\t" << j << "\n";
                    dlist.push_back( Mismatch( rids[i], j-inits[i]-lead_gaps-query_gaps, j-sbjct_gaps) );
                }
            }

            if ( consensus[j] == '-' ) {
                if ( nreads[i][j] != '-' ) {
                    //inits[i] >= 0 ?
                    //ilist.push_back( Mismatch( rids[i], j-query_gaps, j-sbjct_gaps) ) ;
                    if ( verbose ) std::cout << rids[i] << "\tinsertion\t" << j-inits[i]-lead_gaps-query_gaps << "\t" << j-sbjct_gaps << "\t" << j << "\n";
                    ilist.push_back( Mismatch( rids[i], j-inits[i]-lead_gaps-query_gaps, j-sbjct_gaps) ) ;
                }
            }
            
            // increment subject gaps count
            if ( consensus[j] == '-' ) sbjct_gaps++; 

            // increment query gaps count
            if ( nreads[i][j] == '-' ) query_gaps++;
            
            j++;
            
        }
    }
    spath.setMismatches(ilist, INSERTION);
    spath.setMismatches(dlist, DELETION);
}


int MSA::__getMinAlignStart( SpaPath &spath )
{
    int nreads = spath.getReadCount();
    int *inits = spath.getInits();

    int min = 1000000;
    for ( int i = 0; i < nreads; i++ )
        if ( inits[i] < min ) 
            min = inits[i];

    if ( min < 0 ) return min;
    return 0;
}

FlagPair MSA::trimEnds( SpaPath &spath, int base_depth, bool verbose )
{
    bool skip = true;
    bool error = false;

    int lgap = __countLeadingGaps(base_depth);
    if ( verbose && lgap )
        std::cout << "Leading low coverage base:" << lgap << "\n";
    if ( lgap ) {
        __trimReads(lgap, dLEFT);
        __trimProfile(lgap, dLEFT);
        consensus.erase(0, lgap);
        spath.adjustPositions(-1*lgap);
        aligned = false;
        skip = false;
    }
    if ( consensus.size() == 0 ) return FlagPair(skip, true);

    int tgap = __countTrailingGaps(base_depth);
    if ( verbose && tgap )
        std::cout << "Trailing low coverage base:" << tgap << "\n";
    if ( tgap ) {
        __trimReads(tgap, dRIGHT);
        __trimProfile(tgap, dRIGHT);
        assert(tgap <= (int)consensus.size());
        consensus.erase(consensus.size()-tgap, tgap);
        spath.trimIndels(consensus.size());
        aligned = false;
        skip = false;
    }
    if ( consensus.size() == 0 ) return FlagPair(skip, true);
    
    if ( lgap == 0 && tgap == 0 ) skip = true;
    return FlagPair(skip, error );
}

void MSA::__trimReads(int gap, int direction)
{
    assert(gap <= (int)consensus.size());
    for ( size_t i = 0; i < nreads.size(); i++ ) {
        if ( direction == dLEFT )
            nreads[i].erase(0, gap);
        else 
            nreads[i].erase(consensus.size()-gap, gap);
    }
}

void MSA::__trimProfile(int gap, int direction)
{
    //     Profile old;
    //     old.matrix = profile.matrix;


    size_t osize = profile.ncol;
    assert(gap <= (int)osize);

    std::vector< std::vector<unsigned> > old( profile.ncol, std::vector<unsigned>(naas) );
    for ( size_t i = 0; i < profile.ncol; i++ ) {
        for ( size_t j = 0; j < (size_t)naas; j++ ) {
            old[i][j] = profile.matrix[i][j];
        }
    }

    profile.clean();

    profile.ncol = osize - gap;
    //profile.ncol = consensus.size();
    profile.matrix = new unsigned*[profile.ncol];
    for ( size_t i = 0; i < profile.ncol; i++ ) 
        profile.matrix[i] = new unsigned[naas];
    
    for ( size_t i = 0; i < profile.ncol; i++ ) {
        for ( size_t j = 0; j < (size_t)naas; j++ ) {
            if ( direction == dLEFT ) 
                profile.matrix[i][j] = old[i+gap][j];
            else
                profile.matrix[i][j] = old[i][j];
        }
    }
    //old.clean();
}

//==============================================================================
// 
//==============================================================================
void MSA::updateConsensus(SpaPath &spath, int k)
{
    //__resetProfile();
    //__updateProfile();
    makeConsensus();
    KmerArray kmers = biostr::getKmers( biostr::stripGap(consensus), k );
    //KmerArray kmers = biostr::getKmers( consensus, k );
    /* Kmers from ungapped consensus */
    spath.updateKmers( &kmers[0], kmers.size() );

    //spath.updateConsensusSequence( consensus.c_str() );
    spath.updateConsensusSequence( biostr::stripGap(consensus).c_str() );
}

void MSA::printAlignment(std::ostream &out, SpaPath &spath, int csize)
{
    if ( consensus.size() == 0 ) return;

    ReadId* rids = spath.getReads();
    int* inits = spath.getInits();
    unsigned nread = spath.getReadCount();

    std::multimap<int, unsigned> imap;
    for ( unsigned i = 0; i < nread; i++ )
        imap.insert( std::pair<int, unsigned>( inits[i], i ) );

    int from = 0;
    char buf[25];
    while ( from < (int)consensus.size()-1 ) {
        if ( from+csize > (int)consensus.length() )
            csize = consensus.length()-from;

        assert(from+csize <= (int)consensus.length());
        std::string rstr = consensus.substr(from, csize);

        sprintf(buf, "\n%10d", from);
        out << buf << "\t";
        //out << "     " << "\t";
        for ( size_t i = 0; i < rstr.size(); i++ ) {
            if ( (i+1)%5 == 0 && (i+1)%10 != 0 ) out << ".";
            else if ( (i+1)%10 == 0 ) out << ":";
            else out << " ";
        }
        out << "\n";

        sprintf(buf, "%10s", "consensus");
        out << buf << "\t" << rstr << "\n";
        //out << buf << "\t";
        //out << "start\t";
        //out << rstr << "\n";
        for ( std::multimap<int, unsigned>::iterator it = imap.begin(); it != imap.end(); ++it ) {
            std::string read = nreads[it->second];
            assert(from+csize <= (int)read.size());
            std::string sstr = read.substr(from, csize);
            if ( sstr != std::string(sstr.size(), '.') ) {
                sprintf(buf, "%10d", rids[it->second]);
                //out << buf << "\t";
                //char spos[5];
                //sprintf(spos, "%5d", inits[it->second]);
                //out << spos << "\t" << sstr << "\n";
                
                out << buf << "\t" << sstr << "\n"; 
            }
        }
        from += csize;
    }
}


void MSA::fitLength(bool verbose)
{
    size_t len = consensus.size();
    for ( size_t i = 0; i < nreads.size(); i++ ) {
        if ( nreads[i].size() > len ) {
            if ( nreads[i][len] != '.' ) {
                if ( verbose ) std::cout << "Tail trimming read:" << i << "\tAA:" << nreads[i][len] << "\n";
            }
            nreads[i] = nreads[i].substr(0, len);
        } else if ( nreads[i].size() < len ) {
            nreads[i] += std::string( len-nreads[i].size(), '.');
        }
    }
}

bool __blank(std::string &seq, int beg, int end)
{
    for ( int i = beg; i <= end; i++ )
        if ( seq[i] != '.' && seq[i] != '-' ) return false;
    return true;
}

bool __end_right(std::string &seq, int pos, int gap) 
{
    return __blank(seq, pos+gap, seq.size()-1);
}

bool __end_left(std::string &seq, int pos, int gap) 
{
    return __blank(seq, 0, pos-1);
}

void MSA::adjustInsertions( SpaPath &spath,
                            ReadPosPairListMap &rpmap,
                            PosReadCountsMap &prmap,
                            PosCountsMap &cmap, 
                            bool verbose )
{
    if (verbose) std::cout << "Adjusting insertions to reads\n";

    std::tr1::unordered_map<ReadId, std::string> sequence_map;
    int nr = spath.getReadCount();
    ReadId* rids = spath.getReads();
    for ( int i = 0; i < nr; i++ ) 
        sequence_map.insert(std::pair<ReadId, std::string>(rids[i], nreads[i]));

    for ( PosCountsMap::reverse_iterator it = cmap.rbegin(); it != cmap.rend(); ++it ) {
        int rpos = it->first;  // reference position of insertion
        int cins = it->second; // count of insertions

        if ( verbose ) std::cout << "Ref pos:" << rpos << "\tCount:" << cins << "\n";

        std::tr1::unordered_map<ReadId, bool> checked;
        for ( int i = 0; i < nr; i++ ) 
            checked.insert(std::pair<ReadId, bool>(rids[i], false));
        
        ReadCountsMap rmap = prmap[rpos];
        for ( ReadCountsMap::iterator rt = rmap.begin(); rt != rmap.end(); ++rt ) {
            ReadId rid = rt->first;
            int    cin = rt->second;
            
            int diff = cins - cin;
            assert(diff >= 0);
            
            if ( diff > 0 ) {
                if ( verbose ) {
                    std::cout << "Read:" << rid << "\t" << cin << "\t" << diff << "\n";
//                     std::cout << "Sbject:" << consensus << "\n";
//                     std::cout << "Query0:" << sequence_map[rid] << "\n";
                }
                if ( __end_right(sequence_map[rid], rpos, cins) )
                    sequence_map[rid].insert(rpos+cin, diff, '-');
                else if ( __end_left(sequence_map[rid], rpos, cins) )
                    sequence_map[rid].insert(rpos, diff, '-');
                else
                    sequence_map[rid].insert(rpos, diff, '-');
//                 if ( verbose ) std::cout << "Query1:" << sequence_map[rid] << "\n";
            }
            checked[rid] = true;
        }
        
        if ( verbose ) std::cout << "Filling non-insert sequences\n";
        std::tr1::unordered_map<ReadId, bool>::iterator ct;
        for ( ct = checked.begin(); ct != checked.end(); ++ct ) {
            if ( ct->second == false ) {
//                 if ( verbose ) {
//                     std::cout << "Read:" << ct->first << "\n";
//                     std::cout << "Sbject:" << consensus << "\n";
//                     std::cout << "Query0:" << sequence_map[ct->first] << "\n";
//                 }
                sequence_map[ct->first].insert(rpos, cins, '.');
//                 if ( verbose ) std::cout << "Query1:" << sequence_map[ct->first] << "\n";

            }
        }
    }
    
    if ( verbose ) std::cout << "Sequence cleaned\n";
    for ( size_t i = 0; i < nreads.size(); i++ ) {
        nreads[i] = sequence_map[rids[i]];
        int s = 0;
        while ( true ) {
            if ( s >= (int)nreads[i].size() ) break;
            if ( nreads[i][s] != '.' && nreads[i][s] != '-' ) break;
            if ( nreads[i][s] == '-' ) nreads[i][s] = '.';
            s++;
        }
        int e = nreads[i].size()-1;
        while ( true ) {
            if ( e < 0 ) break;
            if ( nreads[i][e] != '.' && nreads[i][e] != '-' ) break;
            if ( nreads[i][e] == '-' ) nreads[i][e] = '.';
            e--;
        }

        for ( int j = s; j <= e; j++ )
            if ( nreads[i][j] == '.' ) nreads[i][j] = '-';
        //if ( verbose )std::cout << rids[i] << "\t" << nreads[i] << "\n";
    }
}

// void MSA::adjustInsertions( SpaPath &spath,
//                                ReadPosPairListMap &rpmap,
//                                PosReadCountsMap &prmap,
//                                PosCountsMap &cmap)
// {
//     std::cout << "Adjusting insertions to reads\n";

//     std::vector<std::string> ostrs = nreads;


//     int nr = spath.getReadCount();
//     ReadId* rids = spath.getReads();
//     int *inits = spath.getInits();

//     std::tr1::unordered_map<ReadId, int> start_poss;
//     std::tr1::unordered_map<ReadId, std::string> read_strs;
    
//     // ---------------------------------------------------------------------
//     // First, temporarily erase all successive insertion letters from reads.
//     // ----- /----------------------------------------------------------------
//     std::tr1::unordered_map<int, std::list<InChar> > maskMap;
//     for ( ReadPosPairListMap::iterator it = rpmap.begin(); it != rpmap.end(); ++it ) {
//         ReadId rid = it->first;
//         PosPairList plist = it->second;
//         for ( int i = 0; i < nr; i++ ) {
//             if ( rid != rids[i] ) continue;

//             std::cout << rid << ":";
//             //for ( PosPairList::iterator jt = plist.begin(); jt != plist.end(); ++jt ) {
//             for ( PosPairList::reverse_iterator jt = plist.rbegin(); jt != plist.rend(); ++jt ) {
//                 int rpos = jt->second;
//                 std::cout << rpos << ",";
//                 assert(rpos >= 0 && rpos < (int)nreads[i].size());
//                 char inch = nreads[i][rpos];
//                 std::cout << inch << " ";
//                 nreads[i].erase(rpos, 1);
//                 //maskMap[i].push_back( InChar(rpos, inch) );
//                 maskMap[i].push_front( InChar(rpos, inch) );
//             }
//             std::cout << "\n"; 
//         }
//     }

//     //-----------------------------------------------------------
//     // Insert gaps to reads
//     // Replace original letters to temporarily masked out letters
//     //-----------------------------------------------------------
//     for ( int i = 0; i < nr; i++ ) {
//         for ( PosCountsMap::reverse_iterator it = cmap.rbegin(); it != cmap.rend(); ++it ) {
//             int rpos = it->first;  // reference position of insertion
//             int cins = it->second; // count of insertions

//             std::cout << rids[i] << "\trpos:" << rpos << "\tcins:" << cins << "\n";

//             assert( rpos >= 0 && rpos < (int)consensus.size() );
//             assert( cins >= 0 );

//             assert( prmap.find(rpos) != prmap.end() );
//             ReadCountsMap rmap = prmap[rpos];

//             std::cout << "\t" << ostrs[i] << "\n";
//             std::cout << "\t" << nreads[i] << "\n";
//             if ( nreads[i][rpos] == '.' )
//                 nreads[i].insert(rpos, cins, '.');
//             else
//                 nreads[i].insert(rpos, cins, '-');
//             std::cout << "\t" << nreads[i] << "\n";

// //             /* Put original letters back to reads */
// //             if ( maskMap.find(i) != maskMap.end() ) {
                
// //                 if ( rpos-cins <= inits[i] ) {
// //                     int pos = cins-1;
// //                     for ( std::list<InChar>::reverse_iterator jt = maskMap[i].rbegin(); jt != maskMap[i].rend(); ++jt ) {
// //                         if ( jt->rpos == rpos ) {
// //                             assert(rpos+pos <= (int)consensus.size());
// //                             nreads[i][rpos+pos] = jt->ch;
// //                             --pos;
// //                         }
// //                     }
// // //                     for ( int j = 0; j <= pos; j++ )
// // //                         if ( nreads[i][rpos+j] == '-' ) nreads[i][j] = '.';
// //                 }
// //                 else {
// //                     int pos = 0;
// //                     for ( std::list<InChar>::iterator jt = maskMap[i].begin(); jt != maskMap[i].end(); ++jt ) {
// //                         if ( jt->rpos == rpos ) {
// //                             assert(rpos+pos <= (int)consensus.size());
// //                             nreads[i][rpos+pos] = jt->ch;
// //                             ++pos;
// //                         }
// //                     }
// // //                     for ( int j = pos; j < cins; j++ )
// // //                         if ( nreads[i][rpos+j] == '-' ) nreads[i][j] = '.';
// //                 }
// //             }
            
// //             /* delete leading gaps */
// //             int j = rpos;
// //             while ( nreads[i][j] == '-' ) {
// //                 if ( j >= rpos+cins ) break; // something not right
// //                 nreads[i][j] = '.';
// //                 j++;
                
// //             }
//             std::cout << "\t" << nreads[i] << "\n";
            
//         }
//     }
// }

void MSA::stretchConsensus( SpaPath &spath, BitString *bstrs, unsigned k, bool verbose)
{
    int maxe = maxEndPosition( spath, bstrs, k, verbose);
    if ( maxe >= (int)consensus.size() ) {
        if ( verbose ) {
            std::cout << "** Stretching consensus:" << maxe-consensus.size()+1 << "\n";
            std::cout << "Old:" << consensus << "\n";
        }
        int off = maxe - consensus.size() + 1;
        consensus += std::string( off, 'X' );
        if ( verbose ) std::cout << "New:" << consensus << "\n";
    }
}

void MSA::stretch( SpaPath &spath, BitString *bstrs, unsigned k, bool verbose )
{
    stretchConsensus( spath, bstrs, k, verbose );

    int* inits = spath.getInits();
    int size   = consensus.size();

    int i = 0;
    for ( SequenceArray::iterator it = nreads.begin(); it != nreads.end(); ++it ) {
        if ( inits[i] > 0 ) {
            assert(inits[i]<size);
            (*it).insert(0, inits[i], '.');
        }
        int diff = size - (*it).size();
        if ( diff > 0 ) {
            assert(diff<size);
            (*it) += std::string(diff, '.');
        }
        i++;
    }
}

// ReadPosPairListMap MSA::loadMismatches( SpaPath &spath, int type )
// {
//     ReadPosPairListMap pmap;

//     Mismatch *miss;
//     int cmis;
//     if ( type == INSERTION ) {
//         miss = spath.getInsertions();
//         cmis = spath.countInsertions();
//     } else {
//         miss = spath.getDeletions();
//         cmis = spath.countDeletions();
//     }
//     for ( int i = 0; i < cmis; i++ ) {
//         Mismatch mis = miss[i];
//         pmap[mis.read].push_back( PosPair(mis.spos, mis.rpos) );
//     }
//     return pmap;
// }

ReadPosPairListMap MSA::loadInsertions( SpaPath &spath, bool verbose )
{
    ReadPosPairListMap pmap;

    Mismatch* inss = spath.getInsertions();
    int       cins = spath.countInsertions();
    for ( int i = 0; i < cins; i++ ) {
        Mismatch ins = inss[i];
        //if ( verbose ) std::cout << "\tInsertion:" << ins.read << "\t" << ins.spos << "\t" << ins.rpos << "\n";
        pmap[ins.read].push_back( PosPair(ins.spos, ins.rpos) );
    }
    return pmap;
}


PosReadCountsMap MSA::countReadPositions(ReadPosPairListMap &pmap)
{
    PosReadCountsMap prmap;
    for ( ReadPosPairListMap::iterator it = pmap.begin(); it != pmap.end(); ++it ) {
        ReadId rid = it->first;
        PosPairList plist = it->second;
        for ( PosPairList::iterator jt = plist.begin(); jt != plist.end(); ++jt ) {
            int rp = (*jt).second; // ref pos
            
            if ( prmap.find(rp) == prmap.end() ) {
                ReadCountsMap rmap;;
                rmap[rid] = 1;
                prmap[rp] = rmap;
            }
            else {
                ReadCountsMap rmap = prmap[rp];
                if ( rmap.find(rid) == rmap.end() ) rmap[rid] = 1;
                else rmap[rid]++;
                prmap[rp] = rmap;
            }
        }
    }

    return prmap;
}

PosCountsMap MSA::getMaxCounts(PosReadCountsMap &prmap)
{
    PosCountsMap cmap;
    for ( PosReadCountsMap::iterator it = prmap.begin(); it != prmap.end(); ++it ) {
        int pos = it->first;
        ReadCountsMap rmap = it->second;
        unsigned max = 0;
        for ( ReadCountsMap::iterator jt = rmap.begin(); jt != rmap.end(); ++jt ) 
            if ( jt->second > max ) max = jt->second;

        cmap[pos] = max;
    }
    return cmap;
}
 
void MSA::insertGapsToReads( SpaPath &spath, BitString *bstrs, bool verbose )
{
    if (verbose) std::cout << "Inserting gaps to reads\n";

    /* deletion map */
    std::tr1::unordered_map<ReadId, std::list<SeqPos> > delmap;
    
    Mismatch *dels = spath.getDeletions();
    int       cdel = spath.countDeletions();
    
    for ( int i = 0; i < cdel; i++ ) {
        Mismatch del = dels[i];
        if ( delmap.find(del.read) == delmap.end()) 
            delmap[del.read] = std::list<SeqPos>();
        delmap[del.read].push_back(del.spos);
    }
    
    SequenceList lreads;
    unsigned nr = spath.getReadCount();
    ReadId* rids = spath.getReads();
    for ( unsigned i = 0; i < nr; i++ ) {
        std::string nread = bstrs[rids[i]].toString();
        if ( delmap.find(rids[i]) == delmap.end() ) {
            lreads.push_back(nread);
            continue;
        }
        
        std::vector<SeqPos> sorted = std::vector<SeqPos>(delmap[rids[i]].begin(), delmap[rids[i]].end());
        std::sort( sorted.begin(), sorted.end() ); 
        delmap[rids[i]] = std::list<SeqPos>(sorted.begin(), sorted.end() );
        //std::sort( delmap[rids[i]].begin(), delmap[rids[i]].end() );

        std::list<SeqPos>::iterator it;
        for ( it = delmap[rids[i]].begin(); it != delmap[rids[i]].end(); ++it ) {
            if ( *it < 0 || *it >= (int)nread.size() ) {
                    std::cerr << "Invalid insert position\t" << *it << "\n";
                    std::cerr << "Read:" << rids[i] << "\t" << nread << "\n";
                    exit(1);
            }
            nread.insert(*it, 1, '-');
            if ( verbose ) std::cout << rids[i] << "\t" << *it << "\t" << nread << "\n";
        }
        lreads.push_back(nread);
    }
    nreads = SequenceArray(lreads.begin(), lreads.end());
}
    


void MSA::insertGapsToRef(PosCountsMap &cmap, SpaPath &spath, int k, bool verbose)
{
    if ( verbose ) std::cout << "Insert gaps to reference\n";

    unsigned nk = spath.getKmerCount();
    KmerId* kmers = spath.getKmers();

    if ( verbose )  {
        if ( cmap.size() ) {
            std::cout << "Insertion positions:" << cmap.size() << "\n";
            for ( PosCountsMap::reverse_iterator it = cmap.rbegin(); it != cmap.rend(); ++it ) 
            //for ( PosCountsMap::iterator it = cmap.begin(); it != cmap.end(); ++it ) 
                std::cout << it->first << ":" << it->second << "\t";
            std::cout << "\n";
        }
    }
    consensus = biostr::getSequenceString( kmers, nk, k );
    if ( verbose ) std::cout << "\nReference sequence\n";
    if ( verbose ) std::cout << "old:" << consensus << "\n";
    for ( PosCountsMap::reverse_iterator it = cmap.rbegin(); it != cmap.rend(); ++it ) {
        //for ( PosCountsMap::iterator it = cmap.begin(); it != cmap.end(); ++it ) {
        if ( it->first < 0 || it->first >= (int)consensus.size() ) {
                std::cerr << "Invalid insertion position:" << it->first << "\n";
                std::cerr << "Consensus:" << consensus << "\n";
                continue; //exit(1);
        } 
        assert(it->second > 0);
        consensus.insert( it->first, it->second, '-');
    }
    if ( verbose ) std::cout << "new:" << consensus << "\n";
}

int MSA::__sumColumn( int col )
{
    int sum = 0;
    for ( int i = 0; i < naas; i++ )
        sum += profile.matrix[col][i];
    return sum;
}

int MSA::__countLeadingGaps(int base_depth)
{
    if ( consensus.size() == 0 ) return 0;

    int count = 0;
    size_t i = 0; 
    while (true) {
        if ( i == consensus.size() ) break;
        if (  __sumColumn(i) >= base_depth ) break;
        count++; i++;
    }
    return count;
}

int MSA::__countTrailingGaps(int base_depth)
{
    if ( consensus.size() == 0 ) return 0;

    int count = 0;
    int i = consensus.size()-1;
    while (true) {
        if ( i < 0 ) break;
        if (  __sumColumn(i) >= base_depth ) break;
        count++; i--;
    }
    return count;
}

std::vector<IntPair> MSA::getZeroCoverageRegions()
{
    std::vector<IntPair> poors;
    size_t i = 0;
    while(true) {
        if ( i == consensus.size() ) break;
        while( __sumColumn(i) ) {        
            i++;
            if ( i == consensus.size() ) return poors;
        }
        size_t s = i;
        i++;
        if ( i == consensus.size() ) return poors;
        while ( __sumColumn(i) == 0 ) {
            i++;
            if ( i == consensus.size() ) return poors;
        }
        size_t e = i;
        poors.push_back( IntPair(s,e) );
    }
    return poors;
}

bool MSA::__flatten( BitString *bstrs, char *strands, ReadId *pairs, SpaPath &spath, Param &param )
{
    if ( !aligned ) __update(spath, param.kmer_size, param.verbose);

    std::vector<IntPair> poors = getZeroCoverageRegions();
    if ( poors.size() == 0 ) return true;

    for ( size_t i = 0; i < poors.size(); i++ ) {
        if ( param.verbose ) std::cout << "Zero coverage:" << poors[i].first << "\t" << poors[i].second << "\n";
    }
    bool success;
    // !pair_flag ? success =__flattenSingleEnds(poors, spath, bstrs, merge_score, offset, seed, verbose ) :
    //     success =  __flattenPairedEnds(poors, spath, bstrs, strands, pairs, merge_score, offset, platform, insert, stddev, verbose);
    !param.pair_flag ? success =__flattenSingleEnds(poors, spath, bstrs, param ) :
        success =  __flattenPairedEnds(poors, spath, bstrs, strands, pairs, param);
    if ( !success ) return false;

    //=========================================================
    // Do MSA again because of possible rearrangement w/ indels
    //=========================================================
    if ( profile.ncol > 0 )  profile.clean();  // clean profile if already exist
    success = pileUpReads(spath, bstrs, param.kmer_size, param.verbose);
    if ( !success ) return false;

    return true;
}

bool MSA::__flattenSingleEnds( std::vector<IntPair> &poors, SpaPath &spath, BitString *bstrs, Param &param )
{
    for ( size_t i = 0; i < poors.size(); i++ ) {
        std::vector<IntPair> sim = __searchSimilarRegion( poors[i], param.verbose );
        if ( param.verbose ) std::cout << "# Similar regions:" << sim.size() << "\n";

        if ( sim.size() == 0 ) {
            if ( param.verbose ) std::cout << "Can't flat reads: 0 similar regions\n"; 
            return false;
        }

        __trimPoorRegion( sim ); 
        if ( sim.size() == 0 ) {
            if ( param.verbose ) std::cout << "Can't flat reads: all region trimmed\n"; 
            return false;
        }

        if ( param.verbose ) std::cout << "Similar region:" << sim.size() << "\n";

        // distribute
        __borrow( poors[i], sim, spath, bstrs, param.merge_score, 2*param.read_spur, param.seed, param.verbose );
    }
    return true;
}

void MSA::__borrow( IntPair bad, std::vector<IntPair> &sim, SpaPath &spath, BitString *bstrs, double merge_score, int offset, unsigned seed, bool verbose )
{
    size_t nsim = sim.size();
    for ( size_t i = 0; i < nsim; i++ ) {
        if ( verbose ) std::cout << "Sim region:" << sim[i].first << "\t" << sim[i].second << "\n";
        std::vector<int> index  = __getReadsInRegion( spath, bstrs, sim[i], offset );
        if ( verbose ) std::cout << "# reads:" << index.size() << "\n";
        std::vector<int> subset = __getRandomSubset( index, index.size()/(nsim+1), seed );
        if ( verbose ) std::cout << "# subset:" << subset.size() << "\n";
        __placeToRegion(bad, subset, spath, bstrs, offset, merge_score, verbose);
    }
}

std::vector<int> MSA::__getRandomSubset( std::vector<int> &index, size_t size, unsigned seed )
{
    srand(seed);
    std::set<int> subset;
    while( subset.size() < size ) {
        int rnum = rand()%index.size();
        subset.insert( index[rnum] );
    }
    return std::vector<int>( subset.begin(), subset.end() );
}


std::vector<int> MSA::__getReadsInRegion( SpaPath &aln, BitString *bstrs, IntPair region, int offset )
{
    int s = region.first - offset;
    int e = region.second - offset;// + offset;
    if ( e < s ) e = s + offset;
    
    //std::cout << "s:" << s << "\te:" << e << "\n";
    //if ( s < 0 ) s = 0;
    //if ( e >= (int)consensus.size() ) e = consensus.size()-1;
    //std::cout << "s:" << s << "\te:" << e << "\n";

    std::vector<int> index;
    int *inits = aln.getInits();
    ReadId *reads = aln.getReads();
    unsigned nread = aln.getReadCount();
    for ( unsigned i = 0; i < nread; i++ ) {
        if ( inits[i] < s ) continue;
        if ( inits[i] > e ) continue;
      
        std::string rstr = bstrs[reads[i]].toString();
        //if ( inits[i] + (int)rstr.size() > e ) continue;
        
        index.push_back(i);
    }
    return index;
}

void MSA::__placeToRegion( IntPair region, std::vector<int> &subset, SpaPath &spath, BitString *bstrs, int offset, double merge_score, bool verbose )
{
    int s = region.first - offset;
    int e = region.second + offset;
    if ( s < 0 ) s = 0;
    if ( e >=(int) consensus.size() ) e = consensus.size()-1;
    
    assert(s <= e);
    std::string sbjct = consensus.substr(s, e-s+1);
    
    ReadIdArray revise;
    std::list<Mismatch>ilist, dlist; 

    int    *inits = spath.getInits();
    ReadId *reads = spath.getReads();
    for ( size_t i = 0; i < subset.size(); i++ ) {
        if ( verbose ) std::cout << "read:" << reads[subset[i]] << "\told init:" << inits[subset[i]] << "\n";
        std::string query = bstrs[reads[subset[i]]].toString();

        
        /* substring match */
        size_t pos = sbjct.find(query);
        if ( pos != std::string::npos ) {
            if ( verbose ) std::cout << "\tSubstring match at:" << pos << "\t";

            inits[subset[i]] = pos + s;
            spath.setInit( subset[i], inits[subset[i]] );
            if ( verbose ) std::cout << "\tnew init:" << inits[subset[i]] << "\n";
            
            __resetRead(subset[i]);
            __assignRead(subset[i], inits[subset[i]], query);
            revise.push_back( reads[subset[i]] );

            continue;
        }


        GlobalAlignPair paln  = GlobalAlignPair(sbjct, query);
        AlignSummary summary = paln.getSummary();
        if ( verbose ) {
            std::cout << "\tScore:" << summary.posrate << "\n";
            paln.printAlignment(std::cout);
        }
        if ( summary.posrate < merge_score ) {
            if ( verbose ) std::cout << "\tWeak replacement:" << summary.posrate << "\n";
            continue;
        }
        inits[subset[i]] = summary.range.first + s;
        spath.setInit( subset[i], inits[subset[i]] );
        if ( verbose ) std::cout << "\tnew init:" << inits[subset[i]] << "\n";
        //std::cout << "new init:" << spath.getInit(subset[i]) << "\n";
        
        __resetRead(subset[i]);
        __assignRead(subset[i], inits[subset[i]], query);
        revise.push_back( reads[subset[i]] );
        
        for ( AlignPosList::iterator it = summary.ilist.begin(); it != summary.ilist.end(); ++it )
            ilist.push_back( Mismatch(reads[subset[i]], it->seq_pos, it->ref_pos + s ) );
        for ( AlignPosList::iterator it = summary.dlist.begin(); it != summary.dlist.end(); ++it )
            dlist.push_back( Mismatch(reads[subset[i]], it->seq_pos, it->ref_pos + s ) );
    }

    spath.dropMismatches( revise );
    spath.appendIndels( ilist, INSERTION );
    spath.appendIndels( dlist, DELETION );
}

void MSA::__resetRead(int index)
{
    for ( size_t i = 0; i < nreads[index].length(); i++ ) 
        nreads[index][i] = '.';
}

void MSA::__assignRead(int index, int start, std::string seq)
{
    nreads[index].replace(start, seq.length(), seq);
}

/**
 * 
 */
bool MSA::__flattenPairedEnds( std::vector<IntPair> &poors, SpaPath &spath, BitString *bstrs, char *strands, ReadId *pairs, Param &param )
{
    if ( param.verbose ) printAlignment( std::cout, spath, 100 );

    // find incorrectly placed read pairs
    // try to place in right place
    // check profile
    // still gap, not paired reads -> distribute
    
    //return __correctMisplacement(spath, bstrs, strands, pairs, merge_score, offset, platform, insert, stddev, verbose);
	return __correctMisplacement(spath, bstrs, strands, pairs, param);
}

bool MSA::__correctMisplacement(SpaPath &spath, BitString *bstrs, char *strands, ReadId *pairs, Param &param )

{
    /* Only Illumina is supported for now. */
    if ( param.platform != ILLUMINA ) return false;

    if ( param.verbose ) std::cout << "Correcting misplaced reads ...\n";

    double min_insert = (param.insert_size - 2*param.insert_sd)/3;
    double max_insert = (param.insert_size + 2*param.insert_sd)/3;
    //double avg_insert = insert/3;

	int offset = param.read_spur * 2;

    std::tr1::unordered_map<int, bool> processed;


    ReadId *reads = spath.getReads();
    int    nread  = spath.getReadCount();
    int    *inits = spath.getInits();
    for ( int i = 0; i < nread; i++ ) {
        if ( processed.find(i) != processed.end() ) continue;
        if ( pairs[reads[i]] == NOT_PAIR ) continue;
        
        int j = findPair( reads, nread, pairs[reads[i]] );
        if ( j == -1 ) continue;
        
        char str1 = strands[reads[i]];
        char str2 = strands[reads[j]];
        
        if ( param.verbose ) {
            std::cout << "\nRead Pairs:\t" << reads[i] << "\t" << reads[j] << "\n";
            std::cout << "Strands:" << str1 << "\t" << str2 << "\n";
            std::cout << "StartPos:" << inits[i] << "\t" << inits[j] << "\n";
        }

        if ( str1 == str2 ) {
            if ( param.verbose ) std::cout << "** SAME STRANDS:\n";
            processed.insert( std::pair<int, bool>( i, true ) );
            processed.insert( std::pair<int, bool>( j, true ) );
            continue;
        } 
        
        size_t l, r;
        if ( str1 == '+' ) {
            l = i; r = j;
        } else {
            l = j; r = i;
        }
        
        int linit = inits[l];
        int rinit = inits[r];     

        //if ( verbose ) std::cout << "left init:" << linit << "\tright init:" << rinit << "\n";
        /*
         * Correct
         ooooooooooooooooooo
         +++++     -----

         * Wrong 1
         ooooooooooooooooooo
         +++++     
         -----         
         * Wrong 2
         ooooooooooooooooooo
         +++++     
         -----         
         * Wrong 3
         ooooooooooooooooooo
         -----    +++++                       
        */


        std::string lstr = bstrs[reads[l]].toString();
        std::string rstr = bstrs[reads[r]].toString();
        size_t lslen = lstr.length();
        size_t rslen = rstr.length();

        bool good = false;
        if ( linit < rinit ) {
            //if ( linit + lslen + min_insert - offset <= rinit && linit + lslen + max_insert + offset >= rinit ) {
            if ( linit + lslen + min_insert <= rinit && linit + lslen + max_insert >= rinit ) {
                good = true;
                if ( param.verbose ) std::cout << "Good placement\n";
            }
            // too close
            else if ( linit + lslen + min_insert > rinit ) {
                if ( param.verbose ) std::cout << "Too close\n";
                bool success = false;

                // Move read left
                if ( linit >= (int)lslen ) {
                    // try to relocate left
                    if ( param.verbose ) std::cout << "Relocating left...\n";
                    assert(linit <= (int)consensus.size());
                    std::string sbjct = consensus.substr(0, linit);
                    GlobalAlignPair paln = GlobalAlignPair( sbjct, lstr );
                    AlignSummary summary = paln.getSummary();
                    if ( summary.posrate >= param.merge_score ) {
                        success = true;
                        if ( param.verbose ) {
                            std::cout << "Left replacement success\n";
                            paln.printAlignment(std::cout);
                            std::cout << "New start:0\n";
                        }
                        spath.updateReadPlacement(l, 0, summary, param.verbose);
                    }
                } 
                // Mover read right
                if ( !success ) {
                    if ( param.verbose ) std::cout << "Relocating right...\n";
                    int beg = linit + lslen + min_insert - offset;
                    int end = linit + lslen + max_insert + offset + rslen;

                    if ( param.verbose ) std::cout << "Subject range\tstart:" << beg << "\tend:" << end << "\n";

                    if ( end < (int)consensus.size() ) {
                        assert(beg>=0);
                        std::string sbjct = consensus.substr(beg, end-beg+1);
                        GlobalAlignPair paln = GlobalAlignPair( sbjct, rstr );
                        AlignSummary summary = paln.getSummary();
                        if ( summary.posrate >= param.merge_score ) {
                            if ( param.verbose ) {
                                std::cout << "Right replacement success\n";
                                paln.printAlignment(std::cout);
                                std::cout << "New start:" << beg << "\n";
                            }
                            spath.updateReadPlacement(r, beg, summary, param.verbose);
                        }
                    }
                    else {
                        if ( param.verbose ) std::cout << "Wrong end:" << end << "\n";
                    }
                }
            } 
            // too far
            else if ( linit + lslen + max_insert < rinit ) {
                if ( param.verbose ) std::cout << "Too far\n";

                // try to shift left read to right
                bool success = false;
                int beg = linit + lslen;
                int end = rinit;
                
                if ( end < (int)consensus.size() ) {
                    assert(beg>=0);
                    assert(beg <= end);
                    std::string sbjct = consensus.substr(beg, end-beg+1);
                    GlobalAlignPair paln = GlobalAlignPair( sbjct, lstr );
                    AlignSummary summary = paln.getSummary();
                    if ( summary.posrate >= param.merge_score ) {
                        success = true;
                        if ( param.verbose ) {
                            std::cout << "Left replacement success\n";
                            std::cout << "New start:" << beg << "\n";
                        }
                        spath.updateReadPlacement(l, beg, summary, param.verbose);
                    }
                }
                else {
                    if ( param.verbose ) std::cout << "Wrong end:" << end << "\n";
                }
                
                // try to shift right read to left
                if ( !success ) {
                    int beg = linit + lslen;
                    int end = rinit;
                    if ( param.verbose ) std::cout << "Subject range\tstart:" << beg << "\tend:" << end << "\n";
                    assert(beg>=0);
                    assert(beg <= end);
                    std::string sbjct = consensus.substr(beg, end-beg+1);
                    GlobalAlignPair paln = GlobalAlignPair( sbjct, rstr );
                    AlignSummary summary = paln.getSummary();
                    if ( summary.posrate >= param.merge_score ) {
                        if ( param.verbose ) {
                            std::cout << "Right replacement success\n";
                            std::cout << "New start:" << beg << "\n";
                        }
                        spath.updateReadPlacement(r, beg, summary, param.verbose);
                    }
                }
            }
        }
        else {
            if ( param.verbose ) std::cout << "Inverse start positions -> skip\n";
            //             int lbeg = linit + lslen;
            //             int lend = rinit;
            //             std::string sbjct = consensus.substr(beg, end-beg+1);
            //             GlobalAlignPair paln = GlobalAlignPair( sbjct, rstr );
            //             AlignSummary summary = paln.getSummary();
            
        }
        processed.insert( std::pair<int, bool>( i, true ) );
        processed.insert( std::pair<int, bool>( j, true ) );
    }
    return true;
}

int MSA::findPair( ReadId *reads, size_t nread, ReadId pread )
{
    for ( int i = 0; i < (int)nread; i++ ) 
        if ( reads[i] == pread ) 
            return i;
    return -1;
}

void MSA::__redistribute()
{

}

void MSA::__trimPoorRegion( std::vector<IntPair> &regions )
{
    std::list<intPair> good;
    std::vector<IntPair>::iterator it;
    for ( it = regions.begin(); it != regions.end(); ++it ) {
        bool zero_depth = false; 
        std::cout << "Region range:" << it->first << "\t" << it->second << "\n";
        for ( int i = (*it).first; i <= (*it).second; i++ ) {
            if (  __sumColumn(i) == 0 ) {
                std::cout << "Zero depth at " << i << "\n";
                zero_depth = true;
                break;
            }
        }
        if ( !zero_depth ) good.push_back(*it);//regions.erase(it++);
        //else ++it;
    }
    regions = std::vector<IntPair>( good.begin(), good.end() );
}

std::vector<IntPair> MSA::__searchSimilarRegion( IntPair qr, bool verbose )
{
    int s = qr.first;
    int l = qr.second-qr.first+1;
    assert(s+l <= (int)consensus.size());
    assert(s>=0);
    std::string query = consensus.substr(s,l);
    if (verbose) std::cout << "Zero coverage sequence:" << query << "\n";
    std::vector<IntPair> r1 = __getSimilarRegions(query, qr.first, qr.second, dLEFT, verbose);
    std::vector<IntPair> r2 = __getSimilarRegions(query, qr.first, qr.second, dRIGHT,verbose);
    //std::vector<IntPair> comb;
    //std::merge( r1.begin(), r1.end(),  r2.begin(), r2.end(), comb.begin() );
    //std::cout << "merged << "\n";
    for ( size_t i = 0; i < r2.size(); i++ ) 
        r1.push_back(r2[i]);
    return r1;
}

std::vector<IntPair> MSA::__getSimilarRegions( std::string &query, int qs, int qe, int direction, bool verbose )
{
    std::vector<IntPair> regions;

    int ss, se, sl;
    if ( direction == dLEFT )  { 
        ss = 0; se = qs-1; sl = qs; 
    }
    else { 
        ss = qe+1; se = consensus.size()-1; sl = se-ss+1; 
    }
    if ( sl < (int)query.length() ) return regions;

    int offset = ss;
    while ( true ) {
        assert(ss>=0);
        assert(ss+sl <= (int)consensus.size());
        std::string sbjct = consensus.substr(ss, sl);
        if ( verbose ) std::cout << "ss:" << ss << "\tsl:" << sl << "\tsbjct:" << sbjct << "\n";
        GlobalAlignPair aln  = GlobalAlignPair(sbjct, query);
        AlignSummary summary = aln.getSummary();
        
        int rs = offset + summary.range.first;
        int re = offset + summary.range.second;
        if ( re >= (int)consensus.length() ) break;

        if ( summary.posrate == 1 ) {
            regions.push_back( IntPair(rs, re) );
        }

        offset += ss;
        ss = rs+1;
        sl = se-ss+1;
        if ( direction == dLEFT  && ss >  qs ) break;
        if ( direction == dRIGHT && ss >= se ) break;
        if ( sl < (int)sbjct.length() ) break;
    }
    return regions;
}


bool MSA::trimAfterStopCodon( SpaPath &spath, bool verbose )
{
    size_t pos = consensus.find("*");
    if ( pos == std::string::npos || pos == consensus.size()-1 ) 
        return false;

    while (true) {
        size_t pos = consensus.find("*");
        if ( pos == std::string::npos || pos == consensus.size()-1 ) break;

        aligned = false; // need to MSA 

        if ( verbose ) 
            std::cout << "Stop Pos:" << pos << "/" << consensus.size() << "\n";
        
        pos += 1;
        int len  = consensus.length()-pos;
        assert(len > 0);

        std::list<int> bad;
        std::string gaps(len, '.');
        for ( size_t i = 0; i < nreads.size(); i++ ) {
            assert( nreads[i].size() >= pos+len );
            assert(pos>=0);
            std::string right = nreads[i].substr(pos, len);
            if ( right != gaps ) bad.push_back(i);
        }    

        __trimReads(len, dRIGHT);
        __trimProfile(len, dRIGHT);        
        makeConsensus();
        spath.trimIndels(consensus.size());
        if ( bad.size() > 0 ) {
            dropReads(bad);
            spath.dropReads(bad);
        }
    }
    return true;
}

void MSA::dropReads( std::list<int> &bad )
{
    //std::cout << "Old reads:" << nreads.size() << "\n";
    std::map<int, bool> good;
    for ( size_t i = 0; i < nreads.size(); i++ ) 
        good.insert( std::pair<int, bool>(i, true) );

    for ( std::list<int>::iterator it = bad.begin(); it != bad.end(); ++it )
        good[*it] = false;

    SequenceList ngood;
    for ( std::map<int, bool>::iterator it = good.begin(); it != good.end(); ++it ) {
        if ( it->second == true ) 
            ngood.push_back( nreads[it->first] );
    }
    nreads = SequenceArray(ngood.begin(), ngood.end());
    //std::cout << "New reads:" << nreads.size() << "\n";
}

void MSA::joinUsedReads( PathId pid,
						 ReadId *reads,
						 int nread,
						 PathId *used_reads )
{
    for ( int i = 0; i < nread; i++ ) 
        used_reads[reads[i]] = pid;
}

void MSA::joinUsedReads( PathId pid,
						 ReadIdArray &path_rids,
						 PathId *used_reads )
{
    joinUsedReads( pid, &path_rids[0], path_rids.size(), used_reads );
}
