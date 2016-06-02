#include "path.h"

std::string GraphPath::toString(int k) 
{
    return biostr::getSequenceString( &kmers[0], kmers.size(), k );
}

void GraphPath::dump(std::fstream &out) 
{
    out.write( (char*)&seed, sizeof(KmerId) );
    out.write( (char*)&spos, sizeof(unsigned) );

    size_t size = kmers.size();
    out.write( (char*)&size, sizeof(size_t) );
    for ( KmerArray::iterator it = kmers.begin(); it != kmers.end(); ++it )
        out.write( (char*)&(*it), sizeof(KmerId) );
 
    out.write( (char*)&(lstop), sizeof(int) );
    out.write( (char*)&(rstop), sizeof(int) );

    // SuffixBounds::iterator bt;
    // for (  bt = bounds.begin(); bt != bounds.end(); ++bt ) {
    //     BoundArray b = *bt;
    //     for ( size_t i = 0; i < b.size(); i++ )  
    //         out.write( (char*)&(b[i]), sizeof(BoundType) );            
    // }
}
void GraphPath::load(std::fstream &in, size_t nparts) 
{
    in.read( (char*)&seed, sizeof(KmerId) );
    in.read( (char*)&spos, sizeof(unsigned) );

    size_t size;
    in.read( (char*)&size, sizeof(size_t) );
    kmers.reserve(size);
    for ( size_t i = 0; i < size; i++ ) {
        KmerId kid;
        in.read( (char*)&(kid), sizeof(KmerId) );
        kmers.push_back(kid);
    }
    in.read( (char*)&(lstop), sizeof(int) );
    in.read( (char*)&(rstop), sizeof(int) );

    // for ( size_t i = 0; i < size; i++ ) {
    //     BoundArray b = BoundArray(nparts, BoundType());
    //     for ( size_t j = 0; j < nparts; j++ ) {
    //         BoundType t;
    //         in.read( (char*)&t, sizeof(BoundType) );
    //         b[j] = t;
    //     }
    //     bounds.push_back(b);
    // }
}

void GraphPath::clear()
{
    path.clear();
    kmers.clear();
    bounds.clear();
    traces.clear();
    used_begs.clear();
    used_ends.clear();

    spos = -1;
    lstop = rstop = 0;
}

void GraphPath::trim( int size, bool left )
{
    NodeArray t_nodes = NodeArray( path.begin(), path.end() );
    assert(size < (int)t_nodes.size());
    assert( kmers.size() == t_nodes.size() );

    std::vector<BoundArray> t_bounds = std::vector<BoundArray>( bounds.begin(), bounds.end());
    assert( t_bounds.size() == t_nodes.size() );
    std::vector<size_t>     t_traces = std::vector<size_t>( traces.begin(), traces.end() );
    assert( t_traces.size() == t_nodes.size() );

    std::vector<UsedBoundInfoList> t_begs = std::vector<UsedBoundInfoList>( used_begs.begin(), used_begs.end() );
    std::vector<UsedBoundInfoList> t_ends = std::vector<UsedBoundInfoList>( used_ends.begin(), used_ends.end() );

    if ( left ) {
        kmers    = KmerArray( kmers.begin()+size, kmers.end() );
        t_nodes  = NodeArray( t_nodes.begin()+size, t_nodes.end() );
        t_bounds = std::vector<BoundArray>( t_bounds.begin()+size, t_bounds.end() );
        t_traces = std::vector<size_t>( t_traces.begin()+size, t_traces.end() );
        t_begs   = std::vector<UsedBoundInfoList>( t_begs.begin()+size, t_begs.end() );
        t_ends   = std::vector<UsedBoundInfoList>( t_ends.begin()+size, t_ends.end() );
    }
    else {
        kmers    = KmerArray( kmers.begin(), kmers.end()-size );
        t_nodes  = NodeArray( t_nodes.begin(), t_nodes.end()-size );
        t_bounds = std::vector<BoundArray>( t_bounds.begin(), t_bounds.end()-size );
        t_traces = std::vector<size_t>( t_traces.begin(), t_traces.end()-size );
        t_begs   = std::vector<UsedBoundInfoList>( t_begs.begin(), t_begs.end()-size );
        t_ends   = std::vector<UsedBoundInfoList>( t_ends.begin(), t_ends.end()-size );
    }
 
    path   = PathType( t_nodes.begin(), t_nodes.end() );
    bounds = SuffixBounds( t_bounds.begin(), t_bounds.end() );
    traces = TraceSizeList( t_traces.begin(), t_traces.end() );
    used_begs = UsedBoundInfoLists( t_begs.begin(), t_begs.end() );
    used_ends = UsedBoundInfoLists( t_ends.begin(), t_ends.end() );
}

void GraphPath::release()
{
    used_begs.clear();
    used_ends.clear();
    bounds.clear();
    traces.clear();
}
