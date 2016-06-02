#include "palign.h"

PathAligner::PathAligner()
{
    dirty = true;
}

PathAligner::PathAligner( ReadAligner &raln )
{
    sequence = raln.getReference();
    ReadPlacementList *p = raln.getPlacements();
    for ( ReadPlacementList::iterator it = p->begin(); it != p->end(); ++it ) {
        if ( Param::verbose > 1 ) it->print(std::cout);
        places.push_back(*it);
    }

    ltrim = raln.getLTrim();
    rtrim = raln.getRTrim();

    lstop = raln.getLStop();
    rstop = raln.getRStop();

    dirty = true;
}

void PathAligner::reset()
{
    sequence = consensus = "";
    places.clear();
    ltrim = rtrim = lstop = rstop = -1;
    
    dirty = true;
}

void PathAligner::addGap( int npos )
{
    for ( ReadPlacementList::iterator it = places.begin(); it != places.end(); ++it ) {
        ReadPlacement place = *it;
        ReadId rid = place.rid;
        int ppos = place.path_pos;
        int rpos = place.read_pos;
        //int len = place.length;
        
        if ( ppos > npos ) continue;
        //if ( param.verbose ) std::cout << "Adding gap to rid:" << rid << "\n";

        int diff = npos-ppos;
        Mismatch ins;
        ins.read = rid;
        ins.ref_pos = rpos+diff;
        ins.qry_pos = diff;
        place.ilist.push_back(ins);
    }    
}

bool PathAligner::dropGaps( int end )
{
    bool found = false;
    for ( auto it = places.begin(); it != places.end(); ++it ) {
        ReadPlacement place = *it;

        if ( place.ilist.size() ) {
            std::list<Mismatch> ins = place.ilist;
            size_t bad = 0;
            for ( auto jt = ins.rbegin(); jt != ins.rend(); ++jt ) {
                if ( jt->ref_pos > end ) bad++;
            }
            if ( bad > 0 ) {
                found = true;
                if ( Param::verbose > 1 ) {
                    std::cout << "# bad gaps:" << bad << "\n";
                }
                while( bad > 0 ) {
                    Mismatch mm = ins.back();
                    if ( Param::verbose > 1 ) 
                        std::cout << "read:" << mm.read << "\tqpos:" << mm.qry_pos << "\trpos:" << mm.ref_pos << "\tend:" << end << "\n";
                    ins.pop_back();
                    bad--;
                }
                place.ilist = ins;
            }
        }
        
        if ( place.dlist.size() ) {
            std::list<Mismatch> del = place.dlist;
            size_t bad = 0;
            for ( auto jt = del.rbegin(); jt != del.rend(); ++jt ) {
                if ( jt->ref_pos > end ) bad++;
            }
            if ( bad > 0 ) {
                found = true;
                while( bad > 0 ) {
                    Mismatch mm = del.back();
                    if ( Param::verbose > 1 ) 
                        std::cout << "read:" << mm.read << "\tqpos:" << mm.qry_pos << "\trpos:" << mm.ref_pos << "\tend:" << end << "\n";
                    del.pop_back();
                    bad--;
                }
                place.dlist = del;
            }
        }
        // update
        if ( found ) *it = place;
    }    
    return found;
}

void PathAligner::adjust( int pos )
{
    //assert( pos >= 0 );
    if ( pos == 0 ) return;

    for ( ReadPlacementList::iterator it = places.begin(); it != places.end(); ++it ) {
        it->path_pos += pos;
        std::list<Mismatch>::iterator mt;
        for ( mt = it->ilist.begin(); mt != it->ilist.end(); ++mt ) 
            mt->ref_pos += pos;
        for ( mt = it->dlist.begin(); mt != it->dlist.end(); ++mt ) 
            mt->ref_pos += pos;
    }
}


void PathAligner::merge( PathAligner &other )
{
    ReadPlacementList::iterator it;
    for ( it = other.places.begin(); it != other.places.end(); ++it ) 
        places.push_back(*it);
}

void PathAligner::realign(char **seq)
{
    int plen = sequence.size();
    if ( Param::verbose > 1 ) {
        std::cout << "Realigning\n";
        std::cout << "Pivot:" << sequence << "\n";
    }
                              
    ReadPlacementList::iterator it;
    for ( it = places.begin(); it != places.end(); ) {
        bool success = false;

        it->ilist.clear();
        it->dlist.clear();

        std::string rstr = seq[it->rid];
        int rlen = rstr.size();

        if ( Param::verbose > 1 ) {
            std::cout << "\nread:" << it->rid << "\t" << rstr << "\n";
            std::cout << "ppos:" << it->path_pos << "\trpos:" << it->read_pos << "\tshift:" << Param::shift_length << "\n";
        }
        
        if ( it->path_pos >= (int)sequence.size() ) {
            places.erase(it++); 
            if ( Param::verbose > 1 ) std::cout << "Wrong path position\n";
            continue;
        }

        // rpos > 0
        if ( it->read_pos > 0 ) {
            if ( it->path_pos < it->read_pos ) {
                it->path_pos = 0;
                it->read_pos -= it->path_pos;
            } else {
                it->path_pos -= it->read_pos;
                it->read_pos = 0;
            }
        }

        if ( Param::verbose > 1 ) {
            std::cout << "ppos:" << it->path_pos << "\trpos:" << it->read_pos << "\n";
        }

        
        // still > 0
        assert( it->read_pos >= 0 && it->read_pos < (int)rstr.size() );
        if ( it->read_pos > 0 ) {
            rstr = rstr.substr(it->read_pos);
            if ( Param::verbose > 1 ) std::cout << "sub read:" << rstr << "\n";
        }

        int s = it->path_pos > Param::shift_length ? it->path_pos-Param::shift_length : 0;
        int e = it->path_pos+rlen+Param::shift_length < plen ? it->path_pos+rlen+Param::shift_length : plen-1;
        if ( Param::verbose > 1 ) {
            std::cout << "s:" << s << "\te:" << e << "\n";
        }

        if ( s < 0 || s >= (int)sequence.size() ) {
            //std::cout << "s:" << s << "\tsequence-len:" << sequence.size() << "\n";
            std::cout << "ppos:" << it->path_pos << "\trpos:" << it->read_pos << "\tlen:" << it->length << "\t#ins:" << it->ilist.size() << "\t#del:" << it->dlist.size() << "\n";
            std::cout << "plen:" << sequence.size() << "\t" << sequence << "\n";
            std::cout << "rlen:" << rstr.size() << "\t" << rstr << "\n";
            std::cout << "s:" <<s << "\n";
            std::cout << "e:" <<e << "\n";
        }

        assert(s >= 0 && s < (int)sequence.size());
        assert(e >= s && e < (int)sequence.size());
        std::string pstr = sequence.substr(s,e+s-1);


        if ( Param::verbose > 1 ) {
            std::cout << "\nrid:" << it->rid << "\n";
            std::cout << "rstr:" << rstr << "\n";
            std::cout << "pstr:" << pstr << "\n";
        }
        std::size_t found = pstr.find(rstr);
        if ( found != std::string::npos ) {
            it->path_pos = s + found;
            //it->read_pos = 0;
            it->length = rlen;
            success = true;
            if ( Param::verbose > 1 ) std::cout << "Substring match at:" << it->path_pos << "\trpos:" << it->read_pos << "\n";
        }
        else {
            AlignSummary sum;
            if ( Param::verbose > 1 ) std::cout << "old:" << it->path_pos << "\n";
            success = alignReadToPath( sum, pstr, rstr, it, s, it->read_pos );
            if ( Param::verbose > 1 ) {
                std::cout << "Recruited?" << success << "\n";
                if ( success ) std::cout << "Aligned at:" << it->path_pos << "\trpos:" << it->read_pos << "\n";
            }
            
        }
        
        success ? ++it : places.erase(it++);
    }
}


bool PathAligner::alignReadToPath( AlignSummary &sum,
                                   std::string &pstr,
                                   std::string &rstr,
                                   //ReadPlacement &p,
                                   ReadPlacementList::iterator &it,
                                   int rbeg,
                                   int qbeg )
{
    int band_size = int ( Param::band_ratio * pstr.size() );
    SemiGlobalAlign aln = ! Param::banded_align ?
        SemiGlobalAlign(pstr, rstr, ANCHOR_CENTER, Param::gap_ext, Param::gap_open) :
        SemiGlobalAlign(pstr, rstr, ANCHOR_CENTER, Param::gap_ext, Param::gap_open, -1*band_size, band_size) ;
    //SemiGlobalAlign aln;
    //aln.setVerbose(Param::verbose);
    //aln = SemiGlobalAlign(pstr, rstr, ANCHOR_CENTER, Param::gap_ext, Param::gap_open);
    sum = aln.getSummary();
    
    if ( Param::verbose > 1 ) aln.printAlignment(std::cout);
    if ( Param::verbose > 1 ) sum.print(std::cout);
    
    if ( sum.posrate < Param::read_align_score ) {
        if ( Param::verbose > 1 ) std::cout << "Weak score:" << sum.posrate << "\n";
        return false;
    }
    if ( (double)sum.length/rstr.size() < Param::read_align_ratio ) {
        if ( Param::verbose > 1 ) std::cout << "Short alignment:" << sum.length << "\n";
        return false;
    }
    
    if ( sum.lgap.first > 0 && sum.lgap.second > 0 ) return false;
    it->path_pos = rbeg;
    it->read_pos = qbeg;
    it->length = sum.length;
    if ( sum.lgap.first > 0 ) {
        it->read_pos += sum.lgap.first;
    }
    if ( sum.lgap.second > 0 ) {
        it->path_pos += sum.lgap.second;
    }
    
    AlignPosList::iterator mt;
    for ( mt = sum.ilist.begin(); mt != sum.ilist.end(); ++mt ) {
        //it->ilist.push_back( Mismatch(it->rid, mt->qry_pos, start+mt->ref_pos) );
        it->ilist.push_back( Mismatch(it->rid, qbeg+mt->qry_pos, rbeg+mt->ref_pos) );
        if ( Param::verbose > 1 ) std::cout << "Ins:\t" << it->rid << "\t" << qbeg+mt->qry_pos << "\t" << rbeg+mt->ref_pos << "\n";
    }
    for ( mt = sum.dlist.begin(); mt != sum.dlist.end(); ++mt ) {
        //it->dlist.push_back( Mismatch(it->rid, mt->qry_pos, start+mt->ref_pos) );
        it->dlist.push_back( Mismatch(it->rid, qbeg+mt->qry_pos, rbeg+mt->ref_pos) );
        if ( Param::verbose > 1 ) std::cout << "Del:\t" << it->rid << "\t" << qbeg+mt->qry_pos << "\t" << rbeg+mt->ref_pos << "\n";
    }
    
    return true;
}


void PathAligner::printPlacement( std::ostream &out )
{
    // out << "ID:" << count << "\n";
    // out << "Raw:" << spath->getConsensus() << "\n";
    out << "Count:" << places.size() << "\n";

    size_t i;
    ReadPlacementList::iterator it;
    std::list<Mismatch>::iterator mt;
    for ( it = places.begin(), i=0; it != places.end(); ++it, ++i ) {
        out << i << "\t" << it->rid << "\t" << it->read_pos << "\t" << it->path_pos;
        
        out << "\tI:";
        for ( mt = it->ilist.begin(); mt != it->ilist.end(); ++mt )
            out << mt->qry_pos << ",";

        out << "\tD:";
        for ( mt = it->dlist.begin(); mt != it->dlist.end(); ++mt )
            out << mt->qry_pos << ",";

        out << "\n";
    }
}

void PathAligner::dump( std::fstream &out )
{
    size_t len = sequence.size();
    out.write((char*)&len, sizeof(size_t));
    if ( len ) {
        const char *cstr = sequence.c_str();
        out.write(cstr, len);
    }

    len = consensus.size();
    out.write((char*)&len, sizeof(size_t));
    if ( len ) {
        const char *cstr = consensus.c_str();
        out.write(cstr, len);
    }

    size_t count = places.size();
    out.write((char*)&count, sizeof(size_t));
    for ( ReadPlacementList::iterator it = places.begin(); it != places.end(); ++it )
        it->dump(out);

    out.write((char*)&lstop, sizeof(int));
    out.write((char*)&rstop, sizeof(int));

    out.write((char*)&ltrim, sizeof(int));
    out.write((char*)&rtrim, sizeof(int));

    out.write((char*)&dirty, sizeof(char));
}

void PathAligner::load( std::fstream &in )
{
    size_t len;
    in.read((char*)&len, sizeof(size_t));
    if ( len ) {
        char *seq = new char[len+1];
        in.read(seq, len);
        seq[len]='\0';
        sequence = std::string(seq);
        delete[] seq;
    }
    
    in.read((char*)&len, sizeof(size_t));
    if ( len ) {
        char *seq = new char[len+1];
        in.read(seq, len);
        seq[len]='\0';
        consensus = std::string(seq);
        delete[] seq;
    }

    size_t count;
    in.read((char*)&count, sizeof(size_t));
    for ( size_t i = 0; i < count; i++ ) {
        ReadPlacement p;
        p.load(in);
        places.push_back(p);
    }

    in.read((char*)&lstop, sizeof(int));
    in.read((char*)&rstop, sizeof(int));

    in.read((char*)&ltrim, sizeof(int));
    in.read((char*)&rtrim, sizeof(int));

    char d;
    in.read((char*)&d, sizeof(char));
    dirty = (bool)d;
}


std::vector<size_t> PathAligner::computeBaseDepth( char **seqs )
{
    int plen = (int)sequence.size();
    std::vector<size_t> depths(plen,0);

    //std::cout << "plen:" << plen << "\n";
    //std::cout << sequence << "\n";
    ReadPlacementList::iterator it;
    for ( it = places.begin(); it != places.end(); ++it ) {
        ReadPlacement place = *it;
        ReadId rid = place.rid;
        int rpos = place.read_pos;
        int ppos = place.path_pos;
        //int alen = place.length;
        
        //std::cout << "rid:" << rid << " rpos:" << rpos << " ppos:" << ppos << " alen:" << alen << "\n";
        std::string rseq = seqs[rid];
        //std::cout << "seq:" << rseq << "\n";
        int rlen = (int)rseq.size();
        int i,k;
        //for ( i=rpos,k=0; i<alen && k<plen; i++,k++ ) {
        for ( i=rpos,k=0; i<rlen && k<plen; i++,k++ ) {
            //if ( i >= rlen ) break;
            if ( k+ppos >= plen ) break;
            depths[k+ppos]++;
        }
    }
    return depths;
}
