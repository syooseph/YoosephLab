#include "PathKmerIndex.h"

PathKmersMap PathKmerIndex::getMap() const
{
    return index;
}

void PathKmerIndex::clear() 
{
    index.clear();
}

void PathKmerIndex::add( const PathId pid, const KmerArray &kmers )
{
    if ( has(pid) ) index[pid] = kmers;
    else index.insert( std::pair<PathId, KmerArray>( pid, kmers ) );
}

void PathKmerIndex::erase( const PathId pid )
{
    PathKmersMap::const_iterator it = index.find(pid);
    assert( it != index.end() );
    index.erase(it);
}

void PathKmerIndex::update( const PathId pid, const KmerArray &kmers )
{
    erase(pid);
    add(pid, kmers);
}

KmerArray PathKmerIndex::getValue( const PathId pid ) const
{
    PathKmersMap::const_iterator it = index.find(pid);
    assert( it != index.end() );
    return it->second;
}

const KmerArray* PathKmerIndex::getPointer( const PathId pid )
{
    PathKmersMap::const_iterator it = index.find(pid);
    assert( it != index.end() );
    return &it->second;
}

bool PathKmerIndex::has(const PathId pid) const
{
    return index.find(pid) != index.end();
}
