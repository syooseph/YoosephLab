#include "seed.h"

void Seed::insert( Vertex node, BinType bin )
{
    assert( has(node, bin) == false );

    BinMemberMap::iterator it = bin2nodes.find(bin);
    if ( it == bin2nodes.end() ) {
        NodeFlagMap nfmap; 
        nfmap.insert( std::pair<Vertex, bool>(node, true) );
        bin2nodes.insert( std::pair<BinType, NodeFlagMap>(bin, nfmap) );
    } else {
        bin2nodes[bin].insert( std::pair<Vertex, bool>(node, true) );
    }

    node2bins.insert( std::pair<Vertex,BinType>(node,bin) );

}

bool Seed::has( Vertex node, BinType &bin )
{
    VertexToBinMap::iterator it = node2bins.find(node);
    if ( it == node2bins.end() ) return false;

    bin = it->second;
    NodeFlagMap::iterator jt = bin2nodes[bin].find(node);
    assert( jt != bin2nodes[bin].end() );

    return true;
}

void Seed::update( Vertex node, BinType bin )
{
    BinType obin;
    assert ( has(node, obin) );
    
    // If graph is updated with deleted edge, 
    // it is possible that bin>=obin
    //assert( obin >= bin );
    if ( obin == bin ) return;
    else {
        erase( node );
        insert( node, bin);
    }
}

void Seed::erase( Vertex node )
{
    BinType bin;
    assert( has(node, bin) == true );

    bin2nodes[bin].erase(node);
    if ( bin2nodes[bin].size() == 0 ) 
        bin2nodes.erase(bin);

    node2bins.erase(node);
}

Vertex Seed::getSeed()
{
    if ( node2bins.size() == 0 ) return NULL;
    
    BinMemberMap::reverse_iterator rit = bin2nodes.rbegin();
    assert( rit->second.size() > 0 );
    
    NodeFlagMap::iterator nit = rit->second.begin();
    return nit->first;
}

Vertex Seed::pop()
{
    if ( node2bins.size() == 0 ) return NULL;
    
    BinMemberMap::reverse_iterator rit = bin2nodes.rbegin();
    assert( rit->second.size() > 0 );
    NodeFlagMap::iterator nit = rit->second.begin();
    Vertex node = nit->first;
    
    if ( Param::verbose > 1 ) std::cout << "Bin:" << getBin(node);
    erase(node);
    
    return node;
}

void Seed::clear()
{
    node2bins.clear();
    bin2nodes.clear();
}

size_t Seed::getSize()
{
    return node2bins.size();
}

BinType Seed::getBin(Vertex node)
{
    BinType bin;
    assert( has(node, bin) == true );
    return bin;
}

size_t Seed::getBinSize()
{
    return bin2nodes.size();
}
