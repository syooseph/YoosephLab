#include "cluster.h"

//====================================================================
// Add a new member to a cluster
//====================================================================
void Cluster::add(PathId &p, AlignSummary &a )
{
    members.push_back(p);
    summarys.push_back(a);
}

//====================================================================
// Find maximum leading or trailing gap by looking into alignment
// summaries.
// This maximum gap will be used to prepend/append new extra bases
// to pivot sequence.
//====================================================================
void Cluster::findMaxGap( size_t &max_gap, PathId &max_pid, bool left )
{
    max_gap = 0;

    PathIdList::iterator it;
    AlignSummaryList::iterator jt;
    for ( it  = members.begin(), jt  = summarys.begin();
          it != members.end(),   jt != summarys.end();
          ++it, ++jt ) {

        PathId pid = *it;

        if ( left ) {
            if ( jt->lgap.first > (int)max_gap ) {
                max_gap = jt->lgap.first;
                max_pid = pid;
            }
        } else {
            if ( jt->egap.first > (int)max_gap ) {
                max_gap = jt->egap.first;
                max_pid = pid;
            }
        }
    }
}

//====================================================================
// For each alisummaries, update leading gap position and indel 
// positions.
//    oooooooooooooooooooo center sequence (==pivot)
// xxxooooo                member 1
//   xooooo                member 2 
//--------------------------------------------------------------------
// After cluster center sequences becomes longer with extra bases 
// (see in merger), it needs to update leading gap position of pivot
// and INDEL positions BECAUSE pivot path is untouched to trace back.
// Only cluster center sequence is changed.
//--------------------------------------------------------------------
// xxxoooooooooooooooooooo center
//    oooooooooooooooooooo pivot
// xxxooooo                member 1
// --xooooo                member 2 
//====================================================================
void Cluster::updatePositions( size_t &max_lgap )
{
    for ( auto jt = summarys.begin(); jt != summarys.end(); ++jt ) {
        
        int lgap = jt->lgap.first;
        assert( lgap >= 0 && lgap <= (int)max_lgap );
        int diff = max_lgap - lgap;
     
        if ( diff > 0 ) {
            jt->lgap.first += diff;
            jt->shift( diff, true );
        }
    }
}

//====================================================================
// Write a cluster object to a file
//====================================================================
void Cluster::dump( std::fstream &out ) 
{
    size_t len = sequence.size();
    const char *cstr = sequence.c_str();
    out.write((char*)&len, sizeof(size_t));
    out.write(cstr, len);
	
    size_t ctmem = members.size();
    out.write((char*)&ctmem, sizeof(size_t));
    for ( auto pid : members ) 
        out.write((char*)&(pid), sizeof(PathId));

    size_t ctsum = summarys.size();
    assert( ctmem == ctsum );
    out.write((char*)&ctsum, sizeof(size_t));
    for ( auto item : summarys ) 
        item.dump(out);

    out.write((char*)&lstop, sizeof(int));
    out.write((char*)&rstop, sizeof(int));
    out.write((char*)&ltrim, sizeof(int));
    out.write((char*)&rtrim, sizeof(int));
    out.write((char*)&lbase, sizeof(int));
    out.write((char*)&rbase, sizeof(int));
}

//====================================================================
// Read a cluster object from a file
//====================================================================
void Cluster::load( std::fstream &in ) 
{
    size_t len;
    in.read((char*)&len, sizeof(size_t));
    assert(len>0);
    char *seq = new char[len+1];
    in.read(seq, len);
    seq[len]='\0';
    sequence = std::string(seq);
    delete[] seq;

    size_t ctmem;
    in.read((char*)&ctmem, sizeof(size_t));
    for ( size_t i = 0; i < ctmem; i++ ) {
        PathId pid;
        in.read((char*)&pid, sizeof(PathId));
        members.push_back(pid);
    }

    size_t ctsum;
    in.read((char*)&ctsum, sizeof(size_t));
    for ( size_t i = 0; i < ctsum; i++ ) {
        AlignSummary summary;
        summary.load(in);
        summarys.push_back( summary );
    }

    in.read((char*)&lstop, sizeof(int));
    in.read((char*)&rstop, sizeof(int));
    in.read((char*)&ltrim, sizeof(int));
    in.read((char*)&rtrim, sizeof(int));
    in.read((char*)&lbase, sizeof(int));
    in.read((char*)&rbase, sizeof(int));
}

