#include "KmerPathIndex.h"

/** Getter: inverted index */
KmerPathIndexMap KmerPathIndex::getKmerPathIndex()
{
    return index;
}

/** Reclaim memory */
void KmerPathIndex::clear() 
{
    index.clear();
}


KmerFreqMap KmerPathIndex::getKmerFreqMap(const KmerArray &kmers)
{
    KmerFreqMap count_map;
    for ( KmerArray::const_iterator it = kmers.begin(); it != kmers.end(); ++it ) {
        if ( count_map.find(*it) == count_map.end() ) 
            count_map.insert(std::pair<KmerId, size_t>( *it, 0 ) );
        count_map[*it]++;
    }
    return count_map;
}

void KmerPathIndex::setKmerFreqMap( KmerFreqMap &count_map, KmerArray &kmers)
{
    for ( KmerArray::iterator it = kmers.begin(); it != kmers.end(); ++it ) {
        KmerFreqMap::iterator mt = count_map.find(*it);
        // if ( count_map.find(*it) == count_map.end() ) 
        //     count_map.insert(std::pair<KmerId, size_t>( *it, 0 ) );
        // count_map[*it]++;
        if ( mt == count_map.end() ) 
            count_map.insert(std::pair<KmerId, size_t>( *it, 1 ) );
        else 
            mt->second ++;
    }
}
	
void KmerPathIndex::add( const KmerArray &kmers, PathId &pid )
{
    KmerFreqMap count_map = getKmerFreqMap(kmers);
    for ( KmerFreqMap::iterator it = count_map.begin(); it != count_map.end(); ++it ) {
        KmerId kid = it->first;
        size_t num = it->second;
        // if ( index.find(kid) == index.end() ) {
        //     index.insert( std::pair<KmerId, PathFreqArray>( kid, PathFreqArray() ) );
        //     index[kid].reserve( 1000 );
        // }
        //index[kid].push_back( std::pair<PathId, size_t>( pid, num ) );
        index[kid].push_back( PathFreq( pid, num ) );
        //index[kid].insert(pid);
    }
}


void KmerPathIndex::setPathFreqMap( PathFreqMap &path_freq, KmerFreqMap &count_map )
{    
    for ( KmerFreqMap::iterator it = count_map.begin(); it != count_map.end(); ++it ) {
        KmerId kmer = it->first;
        //size_t freq = it->second;

        KmerPathIndexMap::iterator jt = index.find(kmer);
        if ( jt == index.end() ) continue;
		
        PathFreqList::iterator kt;
        //PathFreqArray::iterator kt;
        //std::set<PathId>::iterator kt;
        for ( kt = jt->second.begin(); kt != jt->second.end(); ++kt ) {
            PathId pid = kt->first;
            size_t num = kt->second;
            // PathId pid = *kt;
            // size_t num = 1;

            //size_t min = freq > num ? num : freq;
            PathFreqMap::iterator mt = path_freq.find(pid);
            // if ( path_freq.find(pid) == path_freq.end() ) 
            //     path_freq[pid] = 0;
            // path_freq[pid] += num;
            if ( mt == path_freq.end() ) 
                path_freq.insert(std::pair<PathId, size_t>(pid,num));
            else mt->second += num;
            //path_freq[pid] += min;
        }
    }
}

void KmerPathIndex::findPaths( PathFreqMap &path_freqs, KmerArray &kmers ) 
{		
    for ( KmerArray::iterator it = kmers.begin(); it != kmers.end(); ++it ) {
        KmerId kmer = *it;
        KmerPathIndexMap::iterator jt = index.find(kmer);
        if ( jt == index.end() ) continue;
		
        PathFreqList::iterator kt;
        for ( kt = jt->second.begin(); kt != jt->second.end(); ++kt ) {
            PathId pid = kt->first;
            size_t num = kt->second;
            
            PathFreqMap::iterator mt = path_freqs.find(pid);
            if ( mt == path_freqs.end() ) 
                path_freqs.insert(std::pair<PathId, size_t>(pid,num));
            else mt->second += num;
        }
    }
}


void KmerPathIndex::findPaths( PathFreqMap &path_freqs, const KmerArray *kmers ) 
{		
    for ( KmerArray::const_iterator it = kmers->begin(); it != kmers->end(); ++it ) {
        KmerId kmer = *it;
        KmerPathIndexMap::iterator jt = index.find(kmer);
        if ( jt == index.end() ) continue;
		
        PathFreqList::iterator kt;
        for ( kt = jt->second.begin(); kt != jt->second.end(); ++kt ) {
            PathId pid = kt->first;
            size_t num = kt->second;
            
            PathFreqMap::iterator mt = path_freqs.find(pid);
            if ( mt == path_freqs.end() ) 
                path_freqs.insert(std::pair<PathId, size_t>(pid,num));
            else mt->second += num;
        }
    }
}

