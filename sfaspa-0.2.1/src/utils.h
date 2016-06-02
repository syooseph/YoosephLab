/** 
 * \file      utils.h
 * \brief     Some utility functions.
 *            
 * \author    Youngik Yang
 * \version   0.2
 * \date      2010-2011
 * \date      modified: Fri 2013-12-13 02:10:26 PM
 * \bug       None.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef UTIL_H
#define UTIL_H

#include <map>
#include <algorithm>
#include <iterator>
#include <vector>
#include <set>
#include <boost/unordered_map.hpp>

/**
 * \brief Hash type container manipulation
 */
namespace util 
{
	/** 
	 * Convert hash to ordered map
	 */
	template<typename K, typename V>
	std::map<K, V> convert_map(const std::tr1::unordered_map<K, V> &in) 
	{
		std::map<K, V> out;
		copy(in.begin(), in.end(), std::inserter(out, out.begin()));
		return out;
	}

	/**
	 * Swap the order of a pair object
	 */
	template<typename K, typename V>
	std::pair<const V, K> flip_pair(const std::pair<const K, V>& p)
	{
		return std::pair<const V, K>(p.second, p.first);
	}

	/** 
	 * Swap key and value of hash
	 */
	template<typename K, typename V>
	std::multimap<V, K> invert_map(const std::map<K, V>& in)
	{
		std::multimap<V, K> out;
		std::transform(in.begin(), in.end(), std::inserter(out, out.begin()), flip_pair<K, V>);
		return out;
	}

	/**
	 * Sort hash by value
	 * Key and value will be swapped
	 */
	template<typename K, typename V>
	std::multimap<V, K> sortByValue(const std::map<K, V> &key_map) 
	{
		return invert_map( key_map );
	}

	/**
	 * Sort hash by value
	 * Key and value will be swapped
	 */
	template<typename K, typename V>
	std::multimap<V, K> sortByValue(const std::tr1::unordered_map<K, V> &unordered) 
	{
		return invert_map( convert_map(unordered) );
	}

	/**
	 * Get the size of values of the given key
	 */
	template<typename K, typename V>
	std::pair<int, const K> get_size(const std::pair<const K, const V> &p)
	{
		return std::pair<int, K>( (p.second).size(), p.first);
	}

	/**
	 * Make hash of key and the size of values
	 */
	template<typename K, typename V>
	std::multimap<int, K> size_map(const std::map<K, V> &in)
	{
		std::multimap<int, K> out;
		std::transform(in.begin(), in.end(), std::inserter(out, out.begin()), get_size<K, V>);
		return out;
	}

	/**
	 * Make hash of key and the size of values
	 */
	template<typename K, typename V>
	std::multimap<int, K> size_map(const std::tr1::unordered_map<const K, V>& in)
	{
		return size_map( convert_map(in) );
	}

	/**
	 * Sort hash by the size of values
	 * Sizes become keys and keys become value
	 */
	template<typename K, typename V>
	std::multimap<int, K> sortBySize(std::map<K, V> &key_map) 
	{
		return size_map(key_map);
	}

	/**
	 * Sort hash by the size of values
	 * Sizes become keys and keys become value
	 */
	
	template<typename K, typename V>
	std::multimap<int, K> sortBySize(std::tr1::unordered_map<K, V> &key_map) 
	{
		return size_map(convert_map(key_map));
	}
	
	
	/**
	 * Extract values from hash
	 */
	template<typename K, typename V>
	std::vector<V> extractValues(std::multimap<K, V> &in)
	{
		std::vector<V> values;
		typename std::multimap<K, V>::const_iterator it;
		for ( it = in.begin(); it != in.end(); ++it )
			values.push_back(it->second);
		return values;
	}


	/**
	 * Extract values from hash
	 */
	template<typename K, typename V>
	std::vector<V> extractValues(std::map<K, V> &in)
	{
		std::vector<V> values;
		typename std::multimap<K, V>::const_iterator it;
		for ( it = in.begin(); it != in.end(); ++it )
			values.push_back(it->second);
		return values;
	}

	/**
	 * Get minimum of two values
	 */
	template<typename K>
	K min(K a, K b) 
	{
		if ( a < b ) return a;
		return b;
	}

	/**
	 * Swap two values
	 */
	template <typename T> 
	void swap( T &a, T &b ) 
	{
		T t = a;
		a = b;
		b = t;
	}

	/** 
	 * Swap two pointers
	 */
	template <typename T> 
	void swapPtr( T *a, T *b ) 
	{
		T t = *a;
		*a = *b;
		*b = t;
	}
}


#endif
