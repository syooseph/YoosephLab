#include "msa.h"

MSA::MSA()
{
}

MSA::MSA(PathAligner *raln, char **r)//, Param &p )
{
    reads = r;
    // use reference sequence not the consensus
    pivot  = raln->getSequence();
    places = raln->getPlacements();
    nreads.reserve((*places).size());
    nstart.reserve((*places).size());
    rids.reserve((*places).size());

    run();
}

MSA::~MSA()
{
}

void MSA::run()
{
    ReadPlacementList::iterator it;
    for ( it = (*places).begin(); it != (*places).end(); ++it ) {
        rids.push_back(it->rid);
        nstart.push_back( it->path_pos );
    }


    /* Find deletion locations in each read and insert gaps */
    insertGapsToReads();

    //---------------------------------------------------------------------------
    // Find insertion locations in each read and corresponding reference location 
    // Key: read-id, Value(seq-pos, ref-pos)
    //---------------------------------------------------------------------------
    ReadPosPairListMap pmap = loadInsertions();

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
    insertGapsToPivot(cmap);

    /* Stretch consensus sequence and reads */
    stretch();

    /* Insert gaps to reads due to sbjct insertions */
    adjustInsertions(pmap, pcmap, cmap);
    
    equalize();
}

void MSA::insertGapsToReads()
{
    ReadPlacementList::iterator it;
    for ( it = (*places).begin(); it != (*places).end(); ++it ) {
        ReadPlacement place = *it;
        ReadId rid = place.rid;

        std::list<Mismatch> dlist = place.dlist;
        std::vector<int> poss; poss.reserve(dlist.size());
        for ( std::list<Mismatch>::iterator mt = dlist.begin(); mt != dlist.end(); ++mt ) {
            assert( mt->read == rid );
            poss.push_back(mt->qry_pos);
        }
        std::sort( poss.begin(), poss.end() );

        std::string str = reads[rid];
        std::vector<int>::reverse_iterator pt;
        for ( pt = poss.rbegin(); pt != poss.rend(); ++pt ) {
            assert(*pt < (int)str.size());
            str.insert(*pt, 1, '-');
        }

        nreads.push_back(str);
    }
}

ReadPosPairListMap MSA::loadInsertions()
{
    ReadPosPairListMap pmap;
    ReadPlacementList::iterator it;
    for ( it = (*places).begin(); it != (*places).end(); ++it ) {
        ReadPlacement place = *it;
        ReadId rid = place.rid;
        std::list<Mismatch> ilist = place.ilist;
        std::vector<int> poss; poss.reserve(ilist.size());
        for ( std::list<Mismatch>::iterator mt = ilist.begin(); mt != ilist.end(); ++mt ) {
            assert( mt->read == rid );
            poss.push_back(mt->ref_pos);
            pmap[rid].push_back( PosPair(mt->qry_pos, mt->ref_pos) );
        }
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

void MSA::insertGapsToPivot(PosCountsMap &cmap)
{
    if ( Param::verbose > 1 ) { 
        std::cout << "Insert gaps to reference\n";
        if ( cmap.size() ) {
            std::cout << "Insertion positions:" << cmap.size() << "\n";
            for ( PosCountsMap::reverse_iterator it = cmap.rbegin(); it != cmap.rend(); ++it ) 
                std::cout << it->first << ":" << it->second << "\t";
            std::cout << "\n";
        }
    }
    
    if ( Param::verbose > 1 ) {
        std::cout << "\nReference sequence\n";
        std::cout << "old:" << pivot << "\n";
    }
    for ( PosCountsMap::reverse_iterator it = cmap.rbegin(); it != cmap.rend(); ++it ) {
        if ( it->first < 0 || it->first >= (int)pivot.size() ) {
            std::cerr << "Invalid insertion position:" << it->first << "\tmax:" << (int)pivot.size()-1 << "\n";
                std::cerr << "Pivot:" << pivot << "\n";
                exit(1);
                //continue;
        } 
        assert(it->second > 0);
        assert(it->first < (int)pivot.size());
        pivot.insert( it->first, it->second, '-');
    }
    if ( Param::verbose > 1 ) std::cout << "new:" << pivot << "\n";
}

void MSA::stretch()
{
    if ( Param::verbose > 1 ) std::cout << "Stretching sequences\n";
    size_t ncolumn = pivot.size();
    if ( Param::verbose > 1 ) std::cout << "Ref size:" << ncolumn << "\n";
    int i;
    if ( Param::verbose > 1 ) std::cout << "#seqs:" << nreads.size() << "\t#places:" << places->size() << "\n";
    ReadPlacementList::iterator it;
    for ( it = (*places).begin(), i=0; it != (*places).end(); ++it, ++i ) {
        ReadPlacement place = *it;
        ReadId rid = place.rid;
        int read_pos = place.read_pos;
        int path_pos = place.path_pos;
        int aln_len  = place.length;
        if ( Param::verbose > 1 ) std::cout << "i:" << i << "\trp:" << read_pos << "\tpp:" << path_pos << "\trid:" << rid << "\taln-len:" << aln_len << "\n";

        if ( Param::verbose > 1 ) std::cout << "Seq:" << reads[rid] << "\n";
        if ( nreads[i].size() > pivot.size() ) {
            if ( Param::verbose > 1 ) std::cout << "Short path\n";
            pivot += std::string(nreads[i].size()-pivot.size(), '.');
            ncolumn = pivot.size();
        }
    }


    for ( it = (*places).begin(), i=0; it != (*places).end(); ++it, ++i ) {
        ReadPlacement place = *it;
        ReadId rid = place.rid;
        int read_pos = place.read_pos;
        int path_pos = place.path_pos;

        if ( Param::verbose > 1 ) std::cout << "i:" << i << "\trp:" << read_pos << "\tpp:" << path_pos << "\trid:" << rid << "\n";

        std::string str = nreads[i];
        if ( Param::verbose > 1 ) std::cout << str << "\n";

        assert( read_pos < (int)str.size() );
        assert( path_pos < (int)ncolumn );

        // read: xxooooooo
        // path:   oooooooxxxxxx
        if ( read_pos > path_pos ) {
            if ( path_pos > 0 ) {
                int diff  = read_pos - path_pos;
                if ( Param::verbose > 1 ) std::cout << "rp:" << read_pos << "\tpp:" << path_pos << "\tdf:" << diff << "\n";
                assert(diff < (int)str.size());
                str = str.substr(diff);
            } else str = str.substr(read_pos);
            int diff = ncolumn-str.size();
            if ( diff > 0 ) str += std::string( ncolumn-str.size(), '.');
            nstart[i] = 0;
        }
        // read:    ooooo
        // path: xxxoooooxxxxx
        else {
            if ( read_pos > 0 ) {
                int diff = path_pos - read_pos;
                assert( diff >= 0 );
                str = std::string(diff, '.') + str;
                nstart[i] = diff;
            } else str = std::string(path_pos, '.') + str;

            int diff = ncolumn-str.size();
            if ( Param::verbose > 1 ) {
                std::cout << "ref:" << pivot << "\n";
                std::cout << "str:" << str << "\n";
                std::cout << "ncolumn:" << ncolumn << "\tstr-len:" << str.size() << "\n";
            }
            if ( diff > 0 ) str += std::string( ncolumn-str.size(), '.');
        } 
        nreads[i] = str;
    }
}

void MSA::adjustInsertions(ReadPosPairListMap &rpmap,
                           PosReadCountsMap &prmap,
                           PosCountsMap &cmap )
{
    if (Param::verbose > 1) std::cout << "Adjusting insertions to reads\n";

    for ( PosCountsMap::reverse_iterator it = cmap.rbegin(); it != cmap.rend(); ++it ) {
        int rpos = it->first;  // reference position of insertion
        int cins = it->second; // count of insertions

        if ( Param::verbose > 1 ) std::cout << "Ref pos:" << rpos << "\tCount:" << cins << "\n";

        size_t i;
        ReadPlacementList::iterator pt;
        for ( pt = (*places).begin(), i = 0; pt != (*places).end(); ++pt, ++i ) {
            ReadPlacement place = *pt;
            ReadId crid = place.rid;
        
            int diff = cins;
            ReadCountsMap rmap = prmap[rpos];
            for ( ReadCountsMap::iterator rt = rmap.begin(); rt != rmap.end(); ++rt ) {
                ReadId rid = rt->first;
                int    cin = rt->second;
                if ( crid == rid ) diff = cins - cin; 
            }
            assert(diff >= 0);
            if ( diff > 0 ) {
                if ( Param::verbose > 1 ) 
                    std::cout << "Read:" << crid << "\t#ins:" << diff << "\n";
                if ( rpos >= (int)nreads[i].size() ) {
                    std::cout << "Invalid rpos:" << rpos << "\tmax:" << (int)nreads[i].size()-1 << "\trid:" << crid << "\n";
                }
                assert(rpos < (int)nreads[i].size());
                char ch = (nreads[i][rpos] == '.') ? '.' : '-';
                nreads[i].insert(rpos, diff, ch);

                // adjust start
                if ( nstart[i] >= rpos ) nstart[i]+=diff;
            }
        }
    }

    for ( size_t i = 0; i < nreads.size(); i++ ) {
        if ( nreads[i].size() == pivot.size() ) continue;
        if ( nreads[i].size() >  pivot.size() ) 
            nreads[i] = nreads[i].substr(0,pivot.size());
        else
            nreads[i] += std::string(pivot.size()-nreads[i].size(), '.');
    }
}

void MSA::print(std::ostream &out, int csize)
{
    std::multimap<int, size_t> imap;
    for ( size_t i = 0; i < nreads.size(); i++ )
        imap.insert( std::pair<int, unsigned>( nstart[i], i ) );

    char buf[25];
    size_t ncol = csize;
    for ( size_t from = 0; from < pivot.length(); from+=csize ) {
        if ( from+ncol > pivot.length() )
            ncol = pivot.length()-from;

        assert(from+ncol <= pivot.length());
        std::string rstr = pivot.substr(from, ncol);

        from == 0 ? sprintf(buf, "%10zu", from) : sprintf(buf, "\n%10zu", from) ;
        out << buf << "\t";
        for ( size_t i = 0; i < rstr.size(); i++ ) {
            if ( (i+1)%5 == 0 && (i+1)%10 != 0 ) out << ".";
            else if ( (i+1)%10 == 0 ) out << ":";
            else out << " ";
        }
        out << "\n";

        sprintf(buf, "%10s", "reference");
        out << buf << "\t" << rstr << "\n";
        std::multimap<int, size_t>::iterator it;
        for (  it = imap.begin(); it != imap.end(); ++it ) {
            std::string read = nreads[it->second];
            assert(from+ncol <= read.size());
            std::string sstr = read.substr(from, ncol);
            
            if ( sstr != std::string(sstr.size(), '.') ) {
                sprintf(buf, "%10d", rids[it->second]);
                out << buf << "\t" << sstr << "\n"; 
            }
        }
    }
}

void MSA::equalize()
{
    size_t plen = pivot.size();
    while(pivot[plen-1]=='.'){
        plen--;
        if ( plen == 0 ) break;
    }
    pivot = pivot.substr(0, plen);
    size_t max = plen;

    for ( size_t i = 0; i < nreads.size(); i++ ) {
        size_t rlen = nreads[i].size();
        if ( nreads[i].size() > max ) {
            for ( size_t j = nreads[i].size()-1; j >= max; j-- ) 
                if ( nreads[i][j] == '.' ) {
                    rlen--;
                } else break;
            if ( rlen < nreads[i].size() ) 
                nreads[i] = nreads[i].substr(0, rlen);
        }
        if ( rlen > max ) max = rlen;
    }
    
    if ( max > plen ) pivot += std::string(max-plen, '.');
    
    for ( size_t i = 0; i < nreads.size(); i++ ) {
        size_t rlen = nreads[i].size();
        if ( rlen < max ) 
            nreads[i] += std::string(max-rlen, '.');
        else if ( rlen > max ) 
            nreads[i] = nreads[i].substr(0, max);
    }
}
