#include "ralign.h"

ReadAligner::ReadAligner()
{
    gpath  = NULL;
    gsa    = NULL;
    seqs   = NULL;
    preads = NULL;
    log    = NULL;
}

ReadAligner::ReadAligner( GraphPath *gp, 
                          GSA *g, 
                          ReadStartCount *c,
                          char **r, 
                          int n, 
                          PathId *u, 
                          ReadAlignLog *l )
{
    double t0 = mytime();

    init(gp,g,c, r,n,u,l);
    align();
    count();
    trim();

    // To do.
    // Can't do the following correctly at this moment.
    updateUsedCount(gp->used_begs, 1);
    updateUsedCount(gp->used_ends, 0);

    log->t_total = mytime()-t0;
}

ReadAligner::~ReadAligner()
{
    if ( gpath != NULL ) gpath = NULL;
	if ( gsa   != NULL ) gsa   = NULL;
	if ( seqs  != NULL ) seqs  = NULL;
	if ( preads  != NULL ) preads  = NULL;
    if ( log   != NULL ) log   = NULL;
    if ( read_starts != NULL ) read_starts = NULL;
}

void ReadAligner::init( GraphPath *gp, 
                        GSA *g, 
                        ReadStartCount *c,
                        char **r, 
                        int n, 
                        PathId *u, 
                        ReadAlignLog *l )
{
    double tic = mytime();

    nreads = n;
    gpath = gp; 
    gsa = g; 
    seqs = r; 
    preads = u;
    log = l;
    read_starts = c;

    ltrim = rtrim = 0;
    size = 0;
    status = false;

    lstop = gp->getLStop();
    rstop = gp->getRStop();

    reference = gpath->toString(Param::kmer_size);
    depths = std::vector<size_t>(reference.size(), 0);


    counts.init( reference.size(), Param::kmer_size+2, Param::back_trace );

    log->t_start = mytime()-tic;
}

void ReadAligner::release()
{
    if ( ranges.size() ) ranges.clear();
    if ( depths.size() ) depths.clear();
    if ( used_begs.size() ) used_begs.clear();
    if ( used_ends.size() ) used_ends.clear();
    if ( dropped.size() ) dropped.clear();
}

void ReadAligner::dropCount()
{
    counts.clear();
}


void ReadAligner::align()
{
    double tic = mytime();
    alignReads();
    log->t_align = mytime()-tic;
}

void ReadAligner::alignReads()
{
    double t0 = mytime();

    extractReads();
    placeReads();

    log->t_reads = mytime()-t0;
}

void ReadAligner::extractReads()
{
    double t0 = mytime();


    if ( Param::verbose ) {
        std::cout << "Seed:" << alpha::IntegerToAminoAcid(gpath->seed, Param::kmer_size)
                  << "\tSeed pos:" << gpath->spos
                  << "\t#nodes:" << gpath->kmers.size()
                  << "\t#SFA bounds:" << gpath->bounds.size()
                  << "\t#traces:" << gpath->traces.size() << "\n";
    }

    SuffixBounds::iterator st;
    TraceSizeList::iterator tt;
    int ki;

    for ( st=gpath->bounds.begin(), tt=gpath->traces.begin(),ki=0;
          st!=gpath->bounds.end(),  tt!=gpath->traces.end();
          ++st, ++tt, ++ki ) {
        
        updateAlignRange( ki, (int)*tt, *st );
    }
    if ( Param::verbose ) std::cout << "Read range extracted:" << mytime()-t0 << " sec\n";

    log->ct_range = ranges.size();
    log->t_range = mytime()-t0;
}

void ReadAligner::updateAlignRange( int ki, int l, BoundArray ba)
{
    if ( Param::verbose > 1 ) std::cout << "kmer-pos:" << ki << "\tspos:" << gpath->spos << "\tlength:" << l << "\n";
    
    //-----------------------------
    // positions and length in path
    //-----------------------------
    int ps = ki, pe = ki, pl = l;
    
    //---------------
    // left extension
    //---------------
    if ( ki <= gpath->spos ) {
        pe = ki+pl-1;
    }
    //----------------
    // right extension
    //----------------
    else {
        pe = ki + Param::kmer_size - 1;
        ps = pe - pl + 1;
    }
    
    assert( pe < (int)reference.size() );
    
    if ( Param::verbose > 1 ) {
        std::cout << "Range start:" << ps << "\tlength:" << pl << "\n";
        std::cout << "Back trace:" << reference.substr(ps, pl) << "\n";
    }

    for ( int i = 0; i < Param::nparts; i++ ) {
        //--------------------------------------------
        // start read number in multiple suffix arrays
        //-------------------------------------------- 
        int sr = i*int(nreads/Param::nparts);                
        
        if ( Param::verbose > 1 ) std::cout << "#part:" << i << " sr:" << sr << " range:(" << ba[i].first << "\t" << ba[i].second << ")\n";
        for ( int j = ba[i].first; j <= ba[i].second; j++ ) {
            
            GsaType s = gsa[i].getAt(j);
            ReadId rid = sr + s.doc;
            int    pos = (int)s.pos;
            assert( (int)rid < nreads );
            
            if ( preads[rid] != NOT_PATH ) continue;
            if ( preads[rid] == BAD_READ ) continue;
            
            //----------------------------------------------
            // Read range which aligned to the path
            //----------------------------------------------
            int rs = pos;
            int re = rs + (abs(pl)-1);
            if ( Param::verbose > 1 ) {
                std::cout << "Read:" << rid << "\tbeg:" << rs << "\tend:" << re 
                          << "\tpos-len:" << (re-rs+1) << "\tseq-len:" << strlen(seqs[rid]) << "\t"
                          << seqs[rid] << "\n";
            }
            
            ReadRangesMap::iterator rit = ranges.find(rid);
            if ( rit == ranges.end() ) {
                ranges.insert(std::pair<ReadId,AlignRange>(rid, AlignRange(ps,pe,rs,re)) );
            }
            else {
                if ( rit->second.path_beg > ps ) rit->second.path_beg = ps;
                if ( rit->second.path_end < pe ) rit->second.path_end = pe;
                if ( rit->second.read_beg > rs ) rit->second.read_beg = rs;
                if ( rit->second.read_end < re ) rit->second.read_end = re;
            }        
        }
    }
}


void ReadAligner::placeReads()
{
    double t0 = mytime();

    if ( Param::verbose ) std::cout << "#reads:" << ranges.size() << "\n";
    for ( ReadRangesMap::iterator it = ranges.begin(); it != ranges.end();  ) {
        ReadId rid = it->first;
        AlignRange pos = it->second;
        std::string rseq = seqs[rid];

        ReadPlacement place; place.rid = rid;
        if ( Param::verbose > 1 ) std::cout << "Read:" << rid << "\n";
        if ( !placeSingleRead( rseq, pos, place ) ) {
            if ( Param::verbose > 1 ) std::cout << "rid:" << rid << " skipped\n";
            ranges.erase(it++);
            dropped.insert( std::pair<ReadId, bool>(rid, true) );
            continue;
        } 
        else {
            if ( Param::verbose > 1 ) std::cout << "rid:" << rid << " placed\n";
            ++it;
        }
        places.push_back(place);
        size++;
    }
    if ( Param::verbose ) std::cout << "Read placed:" << mytime()-t0 << " sec\n";

    log->ct_place = places.size();
    log->t_place = mytime()-t0;
    
    if ( Param::verbose ) std::cout << "#placed reads:" << places.size() << "\t" << mytime()-t0 << "\n";
}

bool ReadAligner::placeSingleRead( //std::string &pseq,
                                   std::string &rseq,
                                   AlignRange  &poss,
                                   ReadPlacement &place )
{
    std::string *pseq = &reference;

    int sr = poss.read_beg, sp = poss.path_beg;
    int er = poss.read_end, ep = poss.path_end;

    if ( Param::verbose > 1 ) {
        std::cout << "sr:" << sr << " sp:" << sp << " er:" << er << " ep:" << ep << " rlen:" << rseq.size() << " plen:" << pseq->size() << "\n";
    }

    assert( sr >= 0 );
    assert( er < (int)rseq.size() );
    assert( er >= sr );
    assert( ep >= sp );
    assert( ep < (int)pseq->size() );

    size_t raln_len = er-sr+1;
    size_t paln_len = ep-sp+1;

    if ( Param::verbose > 1 ) {
        std::cout << "sr:" << sr << " sp:" << sp << " er:" << er << " ep:" << ep
                  << "\trlen:" << raln_len << "\tplen:" << paln_len
                  << "\tRead:" << rseq.substr(sr, raln_len)
                  << "\tPath:" << pseq->substr(sp, paln_len) << "\n";            
    }

    if ( raln_len == rseq.size() && raln_len == paln_len )
        log->ct_same_len++;
    else if ( (double)raln_len/rseq.size() >= Param::read_align_ratio &&
              (double)paln_len/rseq.size() >= Param::read_align_ratio ) 
        log->ct_good_len++;
    else
        log->ct_poor_len++;

    std::string sbjct, query;
    if ( Param::exact_match_only ) {
        if ( raln_len != paln_len ) {
            log->ct_diff_aln++;
            return false;
        }
        
        if ( sr != 0 ) {
            log->ct_read_off++;
            return false;
        }
        
        if ( raln_len != rseq.size() ) {
            log->ct_diff_len++;
            return false;
        }
        sbjct = pseq->substr(sp, paln_len);
        query = rseq;
    }
    else {
        if ( (double)raln_len/rseq.size() < Param::read_align_ratio ) {
            log->ct_short++;
            if ( Param::verbose > 1 ) std::cout << "Short read alignment:" << raln_len << "\n";
            //if ( Param::verbose > 1 ) std::cout << "read:" << rseq << "\n";
            return false;
        }
        if ( (double)paln_len/rseq.size() < Param::read_align_ratio ) {
            log->ct_short++;
            if ( Param::verbose > 1 ) std::cout << "Short path range: (" <<  sp << " " << ep << ")\n";
            return false;
        }
        
        // too short read
        if ( raln_len < (size_t)Param::kmer_size+2 ) {
            log->ct_tiny++;
            if ( Param::verbose > 1 ) std::cout << "Too short alignment length\n";
            return false; 
        }
        sbjct = pseq->substr(sp, paln_len);
        query = rseq.substr(sr, raln_len);
    }

    // std::string sbjct = pseq.substr(sp, paln_len);
    // std::string query = rseq.substr(sr, raln_len);

    // std::string sbjct = pseq->substr(sp, paln_len);
    // std::string query = rseq;

    place.length = raln_len;
    place.read_pos = sr;
    place.path_pos = sp;
    
    size_t pos = sbjct.find(query);
    if ( pos != std::string::npos ) {
        if ( Param::verbose > 1 ) std::cout << "Exact match\n";
        log->ct_substr++;
        return true; 
    }

    // Exact match only
    if ( Param::exact_match_only ) return false;

    AlignSummary summary;
    if ( sbjct.size() == query.size()  ) {
        summary = compareBase(query, sbjct);
        if ( summary.posrate < Param::read_align_score ) {
            log->ct_weak++;
            if ( Param::verbose > 1 ) std::cout << "Weak alignment:" << summary.posrate << "\n";
            return false;
        } else {
            if ( Param::verbose > 1 ) std::cout << "Strong score:" << summary.posrate << "\n";
            log->ct_basecmp++;
            return true;
        }
    }

    // probably not better allow alignment here 
    else {
        doAlignment( summary, query, sbjct );
        if ( summary.posrate < Param::read_align_score ) return false;            
        if ( summary.ilist.size() > 0 || summary.dlist.size() > 0 ) {
            if ( Param::verbose > 1 ) std::cout << "gappy alignment\n";
            return false;
        }
        if ( Param::verbose > 1 ) std::cout << "good alignment:" << query << "\n";

        place.length = summary.length;
        if ( summary.lgap.first > 0 && summary.lgap.second == 0 ) 
            place.read_pos += summary.lgap.first;
        else if ( summary.lgap.first == 0 && summary.lgap.second > 0 ) 
            place.path_pos += summary.lgap.second;
        else return false;

        log->ct_align++;
        return true;
    }
}

bool ReadAligner::doAlignment(AlignSummary &summary, 
                              std::string &query, 
                              std::string &sbjct )
{
    int band_size = int ( Param::band_ratio * sbjct.size() );
    SemiGlobalAlign paln = !Param::banded_align ? 
        SemiGlobalAlign(sbjct, query, ANCHOR_CENTER) :
        SemiGlobalAlign(sbjct, query, ANCHOR_CENTER, Param::gap_ext, Param::gap_open, -1*band_size, band_size) ;
    
    summary = paln.getSummary();

    size_t aln_beg = summary.range.first;
    size_t aln_end = summary.range.second;
    size_t aln_len = aln_end - aln_beg + 1;
                    
    if ( Param::verbose ) {
        std::cout << paln.getAlignment();
        std::cout << "\t#ins:" << summary.ilist.size() << "\t";
        std::cout << "\t#del:" << summary.dlist.size() << "\t";
        std::cout << "\t#mat:" << summary.match << "\t";
        std::cout << "\t#pos:" << summary.positive << "\t";
        std::cout << "\t#mms:" << summary.mismatch << "\n";
        std::cout << "\t#out:" << summary.outer.first << "\t" << summary.outer.second << "\n";
        std::cout << "\t#in:" << summary.range.first << "\t" << summary.range.second << "\n";
    }

    if ( aln_len < size_t(Param::read_align_ratio * query.size()) ) {
        if ( Param::verbose > 1 ) std::cout << "\talignment - short\n";
        return false;
    }
    else if ( aln_len * Param::read_align_ratio  >  query.size() ) {
        if ( Param::verbose > 1 ) std::cout << "\talignment - long\n";
        return false;
    }

    if ( (double)summary.positive <  aln_len * Param::read_align_score ) {
        if ( Param::verbose > 1 ) std::cout << "\talignment - weak\n";
        return false;
    }

    return true;

}

void ReadAligner::updateIndel( ReadId rid,
                               int init,
                               AlignSummary &summary,
                               std::list<Mismatch> &ilist,
                               std::list<Mismatch> &dlist)
{
    AlignPosList::iterator it;
    for ( it = summary.ilist.begin(); it != summary.ilist.end(); ++it ) {
        ilist.push_back( Mismatch(rid, it->qry_pos, init + it->ref_pos) );
    }
    for ( it = summary.dlist.begin(); it != summary.dlist.end(); ++it ) {
        dlist.push_back( Mismatch(rid, it->qry_pos, init + it->ref_pos) );
    }
}

AlignSummary ReadAligner::compareBase(std::string &query, std::string &sbjct)
{
    AlignSummary summary;
    summary.length = query.length();
    summary.positive = scoring::countPositive(query, sbjct, BLOSUM62);
    summary.posrate = (double)summary.positive/summary.length;
    return summary;
}

void ReadAligner::printPlacements()
{
    std::multimap<int, ReadPlacement> omap;
    for ( ReadPlacementList::iterator it = places.begin(); it != places.end(); ++it )
        omap.insert(std::pair<int, ReadPlacement>( it->path_pos, *it ) );

    std::multimap<int, ReadPlacement>::iterator ot;
    for ( ot = omap.begin(); ot != omap.end(); ++ot ) {
        int pos = ot->first;
        ReadPlacement place = ot->second;
        ReadId rid = place.rid;
        int beg = place.read_pos;
        int len = place.length;
        std::string seq = seqs[rid];
        assert( beg+len <= (int)seq.size() );
        std::string sub = seq.substr(beg, len);
        std::cout << ot->second.rid << "\tpos:" << pos 
                  << "\tbeg:" << beg << "\tlen:" << len 
                  << "\t#ins:" << place.ilist.size() 
                  << "\t#del:" << place.dlist.size()
                  << "\tseq:" << sub << "\n";
    }
}

void ReadAligner::computeBaseDepth()
{
    double t0 = mytime();

    assert(reference.size()>0);
    int plen = (int)reference.size();

    ReadPlacementList::iterator it;
    for ( auto place : places ) {
        ReadId rid = place.rid;
        int rpos = place.read_pos;
        int ppos = place.path_pos;

        std::string rseq = seqs[rid];
        int rlen = (int)rseq.size();
        int i,k;

        // include non-aligned region of read bases 
        // because it was passed of alignment threshold
        for ( i=rpos,k=0; i<rlen && k<plen; i++,k++ ) {
            if ( i >= rlen ) break;
            if ( k+ppos >= plen ) break;
            depths[k+ppos]++;
        }
    }
    if ( Param::verbose > 1 ) {
        std::cout << "Base depths\n";
        for ( size_t i = 0; i < depths.size(); i++ ) 
            std::cout << i << ":" << depths[i] << " ";
        std::cout << "\n";
    }
    if ( Param::verbose ) std::cout << "Read coverage computed:" << mytime()-t0 << " sec\n";
}

void ReadAligner::trim()
{
    double tic = mytime();
    status = false;

    if ( places.size() == 0 ) {
        log->t_clean = mytime()-tic;
        return;
    }

    computeBaseDepth();
 
    double t0 = mytime();

    int plen = (int)reference.size();
    int ns,ne;
    for ( ns = 0; ns < plen; ns++ )
        if ( depths[ns] != 0 ) break;
    for ( ne = plen-1; ne > ns; ne-- )
        if ( depths[ne] != 0 ) break;

    if ( Param::verbose > 1 ) 
        printf("trimmed\tns:%d\tne:%d\n", ns, ne);

    if ( ns >= ne ) {
        if ( Param::verbose > 1 ) std::cout << "Invalid range after trimming\n";
        return;
    }

    // find holes
    for ( int i = ns; i <= ne; i++ ) {
        if ( depths[i] == 0 ) {
            if ( Param::verbose ) std::cout << "Found a hole at " << i << "\n";
            log->t_clean = mytime()-tic;
            return;
        }
    }

    rtrim = (plen-1)-ne;
    ltrim = ns;

    if ( ne < plen-1 ) {
        int l = ne+1;
        reference = reference.substr(0, l);
        depths    = std::vector<size_t>(depths.begin(), depths.begin()+l);
        if ( Param::verbose ) std::cout << "Reference trimmed:" << reference << "\n";
    }
    if ( ns > 0 ) {
        reference = reference.substr(ns);
        depths    = std::vector<size_t>(depths.begin()+ns, depths.end());
        adjust(-1*ns);
        if ( Param::verbose ) std::cout << "Reference trimmed:" << reference << "\n";
    }
    
    
    if ( ns != 0 || ne != plen-1 )
        counts.trim(ns,ne);

    if ( Param::verbose ) std::cout << "Read trimmed:" << mytime()-t0 << " sec\n";

    if ( reference.size() == 0 ) {
        if ( Param::verbose ) std::cout << "Zero length reference\n";
        //return false;
        log->t_clean = mytime()-tic;
        return;
    }

    status = true;
    log->t_clean = mytime()-tic;
}


void ReadAligner::adjust( int pos )
{
    for ( ReadPlacementList::iterator it = places.begin(); it != places.end(); ++it ) {
        it->path_pos += pos;
        std::list<Mismatch>::iterator mt;
        for ( mt = it->ilist.begin(); mt != it->ilist.end(); ++mt ) 
            mt->ref_pos += pos;
        for ( mt = it->dlist.begin(); mt != it->dlist.end(); ++mt ) 
            mt->ref_pos += pos;
    }
}

void ReadAligner::dump( std::fstream &out )
{
    size_t len = reference.size();
    out.write((char*)&len, sizeof(size_t));
    if ( len ) {
        const char *cstr = reference.c_str();
        out.write(cstr, len);
    }

    size_t nplace = places.size();
    out.write((char*)&nplace, sizeof(size_t));
    
    ReadPlacementList::iterator it;
    for ( it = places.begin(); it != places.end(); ++it )
        it->dump(out);

    out.write((char*)&lstop, sizeof(int));
    out.write((char*)&rstop, sizeof(int));

    out.write((char*)&ltrim, sizeof(int));
    out.write((char*)&rtrim, sizeof(int));
}

void ReadAligner::load( std::fstream &in )
{
    size_t len;
    in.read((char*)&len, sizeof(size_t));
    if ( len ) {
        char *seq = new char[len+1];
        in.read(seq, len);
        seq[len]='\0';
        reference = std::string(seq);
        delete[] seq;
    }

    in.read((char*)&len, sizeof(size_t));
    for ( size_t i = 0; i < len; i++ ) {
        ReadPlacement p;
        p.load(in);
        places.push_back(p);
    }

    in.read((char*)&lstop, sizeof(int));
    in.read((char*)&rstop, sizeof(int));

    in.read((char*)&ltrim, sizeof(int));
    in.read((char*)&rtrim, sizeof(int));
}


void ReadAligner::getWordFrequency( WordFreqMap &word_freqs,
                                    int wsize,
                                    GSA *g, 
                                    char **r )
{
    if ( places.size() >= (size_t)Param::ncpus ) return getWordFrequencyMP(word_freqs, wsize, g, r, Param::ncpus);

    gsa = g; 
    seqs = r;

    for ( ReadPlacementList::iterator it = places.begin(); it != places.end(); ++it ) {
        ReadId rid = it->rid;
        
        std::string rseq = seqs[rid];
        int rlen = (int)rseq.size();

        int nword = rlen-wsize+1;
        for ( int i = 0; i<nword; i++ ){
            std::string sstr = rseq.substr(i, wsize);
            if ( word_freqs.find(sstr) == word_freqs.end() ) 
                word_freqs.insert(std::pair<std::string, size_t>(sstr,0));
            word_freqs[sstr]++;
        }
    }
}

void ReadAligner::getWordFrequencyMP( WordFreqMap &word_freqs,
                                      int wsize,
                                      GSA *g, 
                                      char **r,
                                      int ncores )
{
    gsa = g; 
    seqs = r;

    WordFreqMap *freqs = new WordFreqMap[ncores];

    typedef std::vector<ReadPlacement> Placements;
    Placements ps = Placements( places.begin(), places.end() );
    

    size_t k, n = ps.size();
#pragma omp parallel for schedule(runtime) if(ncores>1) private(k) num_threads(ncores)
    for ( k = 0; k < n; k++ ) {
        int tid = omp_get_thread_num();
        
        ReadPlacement p = ps[k];
        ReadId rid = p.rid;
        
        std::string rseq = seqs[rid];
        int rlen = (int)rseq.size();
        
        int nword = rlen-wsize+1;
        for ( int i = 0; i<nword; i++ ){
            std::string sstr = rseq.substr(i, wsize);
            WordFreqMap::iterator it = freqs[tid].find(sstr);
            if ( it == freqs[tid].end() ) 
                freqs[tid].insert(std::pair<std::string, size_t>(sstr,1));
            else it->second++;
            // if ( freqs[tid].find(sstr) == freqs[tid].end() ) 
            //     freqs[tid].insert(std::pair<std::string, size_t>(sstr,0));
            // freqs[tid][sstr]++;
        }
    }

    WordFreqMap::iterator it, jt;
    for ( k = 0; k < (size_t)ncores; k++ ) {
        for ( it = freqs[k].begin(); it != freqs[k].end(); ++it ) {
            jt = word_freqs.find(it->first);
            if ( jt == word_freqs.end() ) 
                word_freqs.insert( std::pair<std::string,size_t>(it->first, it->second) );
            else
                jt->second += it->second;
        }
    }
    delete[] freqs;

}

void ReadAligner::count()
{
    if ( places.size() == 0 ) return;
    if ( reference.size() == 0 ) return;

    double tic = mytime();

    counts.init(reference.size(), Param::kmer_size+2, Param::back_trace);
    for ( ReadRangesMap::iterator it = ranges.begin(); it != ranges.end(); ++it ) {
        ReadId rid = it->first;
        AlignRange pos = it->second;

        int sp = pos.path_beg;
        int ep = pos.path_end;
        assert( sp >= 0 );
        assert( ep >= sp );
        assert( ep < (int)reference.size() );

        if ( ep-sp+1 < Param::kmer_size+2) {
            if ( Param::verbose ) {
                printf("[Warning] invalid length for read:%u (beg:%d end:%d len:%d)\t",
                       rid, sp, ep, ep-sp+1);
                std::cout << "str:" << seqs[rid] << "\n";
            }
            continue;
        }
        counts.increment(sp, ep-sp+1, Param::kmer_size+2, Param::back_trace);
    }

    log->t_count = mytime()-tic;
}
 
int ReadAligner::getReadCount( int s, int e )
{
    return counts.get(s,e, Param::kmer_size+2, Param::back_trace);
}


// End bound is indeed the bound of all reads end with '\0'
// Start bound is not correct.
// It need to have start array to be passed
//-------------------------------------------------------------
void ReadAligner::updateUsedCount( UsedBoundInfoLists &used_lists, bool start )
{
    double tic = mytime();
    
    UsedBoundMap *used_bounds = start ? &used_begs : &used_ends;

    // UsedBoundInfoLists
    for ( auto it = used_lists.begin(); it != used_lists.end(); ++it ) {
        // UsedBoundInfoList
        for ( auto jt = it->begin(); jt != it->end(); ++jt ) {
            // UsedBoundInfo
            if ( used_bounds->find(jt->first) != used_bounds->end() )
                continue;
            else
                used_bounds->insert( std::pair<UsedBound, size_t>(jt->first, jt->second ));
        }
    }


    size_t k, n = used_bounds->size();
    std::vector< size_t > Used( n, 0 );
    std::vector< size_t > Nums; Nums.reserve( n );
    std::vector< UsedBound > Bounds; Bounds.reserve( n );
    for ( auto it = used_bounds->begin(); it != used_bounds->end(); ++it ) {
        Bounds.push_back( it->first );
        Nums.push_back( it->second );
    }

    bool par_run = ( Param::ncpus > 1 && n > 1 ) ? true : false;
#pragma omp parallel for schedule(dynamic, 1) if (par_run) private(k) num_threads(Param::ncpus)
    for ( k = 0; k < n; k++ ) {
    
        UsedBound bound = Bounds[k];

    // for ( auto it = used_bounds->begin(); it != used_bounds->end(); ++it ) {
    //     if ( it->first.min > it->first.max ) continue;

        if ( bound.min > bound.max ) continue;

        // start: count fro read_starts
        // end: count from getEndBound
        //int num = it->second; // count
        int num = Nums[k];
        int bad = 0;

        if ( Param::verbose > 1 ) 
            std::cout << "bound: (" << bound.min << ", " << bound.max 
                      << ") count:" << num << "\n";

        for ( int i = bound.min; i <= bound.max; i++ ) {
            LcpType gid  = bound.gsa;
            GsaType item = gsa[gid].getAt(i);
            ReadId rid = (ReadId)item.doc;
            assert( rid >= 0 && rid < (size_t)nreads );

            if ( start ) {
                if ( i == 0 ) {
                    if ( read_starts[gid].get(i) == 0 ) continue;
                }
                else if (i>0) {
                    // Not a starting read at this suffix 
                    if ( read_starts[gid].get(i-1) == read_starts[gid].get(i) ) 
                        continue;
                }
            }

            if ( preads[rid] != NOT_PATH && preads[rid] != BAD_READ ) bad++;
            else if ( dropped.find(rid) != dropped.end() ) bad++;
        }

        if ( Param::verbose > 1 ) std::cout << "num:" << num << "\tbad:" << bad << "\n";
        
        //num -= bad;
        int use = num - bad;
        //it->second = num >= 0 ? num : 0;
        Used[k] = use >= 0 ? use : 0;
    }

    for ( k = 0; k < n; k++ )
        (*used_bounds)[Bounds[k]] = Used[k];

    log->t_update += (mytime()-tic);
}

