#include "PathKmerPosIndex.h"

PathKmerPosMap PathKmerPosIndex::getMap() const
{
    return index;
}

void PathKmerPosIndex::clear() 
{
    index.clear();
}

void PathKmerPosIndex::add( const PathId pid, const KmerPosMap &kmer_map )
{
    if ( has(pid) ) index[pid] = kmer_map;
    else index.insert( std::pair<PathId, KmerPosMap>( pid, kmer_map ) );
}

void PathKmerPosIndex::erase( const PathId pid )
{
    PathKmerPosMap::const_iterator it = index.find(pid);
    assert( it != index.end() );
    index.erase(it);
}

void PathKmerPosIndex::update( const PathId pid, const KmerPosMap &kmer_map )
{
    erase(pid);
    add(pid, kmer_map);
}


KmerPosMap PathKmerPosIndex::getValue( const PathId pid ) const
{
    PathKmerPosMap::const_iterator it = index.find(pid);
    assert( it != index.end() );
    return it->second;
}

const KmerPosMap* PathKmerPosIndex::getPointer( const PathId pid ) const
{
    PathKmerPosMap::const_iterator it = index.find(pid);
    assert( it != index.end() );
    return &it->second;
}

bool PathKmerPosIndex::has(const PathId pid) const
{
    return index.find(pid) != index.end();
}
