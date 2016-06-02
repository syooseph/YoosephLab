#ifndef PAIR_KEY_MAP_H
#define PAIR_KEY_MAP_H

#include <map>

template<class T>
struct KeyInfo
{
	T key1;
	T key2;
	
	KeyInfo( T k1, T k2 ) 
	{
		key1 = k1; key2 = k2;
	}

    bool operator<(const KeyInfo& other) const
    { 
		return key1 < other.key1;
	}
};

template<class K, class V>
struct Map
{
    typedef std::map<K, V> Type;
};


#endif
