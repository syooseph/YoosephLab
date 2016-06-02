/**
 * \file container.h 
 * STL container util.
 * This header includes conversion from one container to the other, 
 * extracting values and keys from map objects.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \bug       Not known.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef CONTAINER_H
#define CONTAINER_H

#include <vector>
#include <queue>
#include <list>
#include <set>

///**
// * Conversion from STL set to STL vector.
// */
//template<typename T> std::vector<T> SetToVector(std::set<T> &set)
//{
//    return std::vector<T>(set.begin(), set.end());
//}
//
///**
// * Conversion from STL list to STL set.
// */
//template<typename T> std::set<T> ListToSet(std::list<T> &l)
//{
//    return std::set<T>(l.begin(), l.end());
//}
//
///**
// * Conversion from STL list to STL vector.
// */
//template<typename T> std::vector<T> ListToVector(std::list<T> &l)
//{
//    return std::vector<T>(l.begin(), l.end());
//}
//
///**
// * Conversion from STL vector to STL list.
// */
//
//template<typename T> std::list<T> VectorToList( std::vector<T> &v )
//{
//    return std::list<T>(v.begin(), v.end());
//}
//
///**
// * Conversion from STL vector to STL set.
// */
//
//template<typename T> std::set<T> VectorToSet( std::vector<T> &v )
//{
//    return std::set<T>(v.begin(), v.end());
//}
//
///**
// * Conversion from STL vector to STL queue.
// */
//template<typename T> std::queue<T> VectorToQueue( std::vector<T> &v )
//{
//    std::queue<T> q;
//    typename std::vector<T>::iterator it;
//    for ( it = v.begin(); it != v.end(); ++it )
//        q.push(*it);
//    return q;
//}
//
///**
// * Conversion from STL queue to STL vector.
// */
//template<typename T> std::vector<T> QueueToVector( std::queue<T> &q )
//{
//    std::vector<T> v;
//    while ( !q.empty() ) {
//        v.push_back(q.front()); q.pop();
//    }
//    return v;
//}

/**
 * Conversion from STL queue to STL list.
 */
template<typename T> std::list<T> QueueToList( std::queue<T> &q )
{
    std::list<T> l;
    while ( !q.empty() ) {
        l.push_back(q.front()); q.pop();
    }
    return l;
}

/**
 * Reclaim memory of queue object.
 */
template<typename T> void clear(std::queue<T> &q)
{
    while ( !q.empty()) q.pop();
}

//
///**
// * Conversion from STL list to STL queue.
// */
//template<typename T> std::queue<T> ListToQueue( std::list<T> &l )
//{
//	std::queue<T> q;
//	typename std::list<T>::iterator it;
//    for ( it = l.begin(); it != l.end(); ++it )
//        q.push(*it);
//
//    return q;
//}
//
///**
// * Conversion from STL set to STL list.
// */
//template<typename T> std::list<T> SetToList( std::set<T> &s )
//{
//	std::list<T> l;
//	typename std::set<T>::iterator it;
//	for ( it = s.begin(); it != s.end(); ++it )
//		l.push_back(*it);
//	return l;
//}

/**
 * Extract keys from map.
 */
template<typename T, typename V> std::vector<T> GetKeys( std::tr1::unordered_map<T,V> &map )
{
	std::list<T> l;
	typename std::tr1::unordered_map<T,V>::iterator it;
	for ( it = map.begin(); it != map.end(); ++it ) 
		l.push_back(it->first);
	return std::vector<T>(l.begin(), l.end());
}

/**
 * Extract values from map.
 */
template<typename T, typename V> std::vector<T> GetVals( std::tr1::unordered_map<T,V> &map )
{
	std::list<T> l;
	typename std::tr1::unordered_map<T,V>::iterator it;
	for ( it = map.begin(); it != map.end(); ++it ) 
		l.push_back(it->second);
	return std::vector<T>(l.begin(), l.end());
}

#endif
