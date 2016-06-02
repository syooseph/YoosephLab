//==============================================================================
// Thu 2010-11-04 05:16:55 PM
//
//==============================================================================

#ifndef UTIL_H
#define UTIL_H

#include <map>
#include <tr1/unordered_map>
#include <algorithm>
#include <iterator>
#include <vector>
#include <set>
#include <boost/unordered_map.hpp>

//using namespace std;

namespace util {
	template<typename K, typename V>
		std::map<K, V> convert_map(const boost::unordered_map<K, V> &in) 
		{
			std::map<K, V> out;
			copy(in.begin(), in.end(), std::inserter(out, out.begin()));
			return out;
		}

	template<typename K, typename V>
		std::pair<const V, K> flip_pair(const std::pair<const K, V>& p)
		{
			return std::pair<const V, K>(p.second, p.first);
		}

	template<typename K, typename V>
		std::multimap<V, K> invert_map(const std::map<K, V>& in)
		{
			std::multimap<V, K> out;
			std::transform(in.begin(), in.end(), std::inserter(out, out.begin()), flip_pair<K, V>);
			return out;
		}

	template<typename K, typename V>
		std::multimap<V, K> sortByValue(const std::map<K, V> &key_map) 
		{
			return invert_map( key_map );
		}

	template<typename K, typename V>
		std::multimap<V, K> sortByValue(const boost::unordered_map<K, V> &unordered) 
		{
			return invert_map( convert_map(unordered) );
		}


	template<typename K, typename V>
		std::pair<int, const K> get_size(const std::pair<const K, const V> &p)
	{
		return std::pair<int, K>( (p.second).size(), p.first);
	}

	template<typename K, typename V>
		std::multimap<int, K> size_map(const std::map<K, V> &in)
		{
			std::multimap<int, K> out;
			std::transform(in.begin(), in.end(), std::inserter(out, out.begin()), get_size<K, V>);
			return out;
		}

	template<typename K, typename V>
		std::multimap<int, K> size_map(const boost::unordered_map<const K, V>& in)
		{
			return size_map( convert_map(in) );
		}


	template<typename K, typename V>
		std::multimap<int, K> sortBySize(std::map<K, V> &key_map) 
		{
			return size_map(key_map);
		}

	
	template<typename K, typename V>
		std::multimap<int, K> sortBySize(boost::unordered_map<K, V> &key_map) 
		{
			return size_map(convert_map(key_map));
		}
	
	

	template<typename K, typename V>
		std::vector<V> extractValues(std::multimap<K, V> &in)
		{
			std::vector<V> values;
			typename std::multimap<K, V>::const_iterator it;
			for ( it = in.begin(); it != in.end(); ++it )
				values.push_back(it->second);
			return values;
		}
	template<typename K, typename V>
		std::vector<V> extractValues(std::map<K, V> &in)
		{
			std::vector<V> values;
			typename std::multimap<K, V>::const_iterator it;
			for ( it = in.begin(); it != in.end(); ++it )
				values.push_back(it->second);
			return values;
		}

	template<typename K>
		K min(K a, K b) 
		{
			if ( a < b ) return a;
			return b;
		}

	template <typename T> void swap( T &a, T &b ) {
		T t = a;
		a = b;
		b = t;
	}

	template <typename T> void swapPtr( T *a, T *b ) {
		T t = *a;
		*a = *b;
		*b = t;
	}

}


#endif
