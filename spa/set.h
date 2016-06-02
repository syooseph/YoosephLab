#ifndef __SET_H__
#define __SET_H__

#include <vector>
#include <set>
//#include <algorithm>

//using namespace std;

namespace Set
{
	template<typename T> std::vector<T> 
	inline Intersect( std::vector<T> &v1,
					  std::vector<T> &v2,
					  bool sorted) 
    {
		if ( !sorted ) {
			std::sort(v1.begin(), v1.end());
			std::sort(v2.begin(), v2.end());
		}
		
		std::vector<T> c(v1.size());
		typename std::vector<T>::iterator it;
		it = std::set_intersection( v1.begin(), v1.end(),
									v2.begin(), v2.end(),
									c.begin() );
		
		return std::vector<T>(c.begin(), c.begin()+int(it-c.begin()));
	}
	
	template<typename T> std::vector<T> 
	inline Union( std::vector<T> &v1,
				  std::vector<T> &v2,
				  bool sorted) 
	{
		if ( !sorted ) {
			std::sort(v1.begin(), v1.end());
			std::sort(v2.begin(), v2.end());
		}
		
		std::vector<T> u(v1.size()+v2.size());
		typename std::vector<T>::iterator it;
		it = std::set_union( v1.begin(), v1.end(),
							 v2.begin(), v2.end(),
							 u.begin() );
		
		return std::vector<T>(u.begin(), u.begin()+int(it-u.begin()));
	}
	
	template<typename T> std::vector<T> 
	inline Difference( std::vector<T> &v1,
					   std::vector<T> &v2,
					   bool sorted) 
    {
		if ( !sorted ) {
			std::sort(v1.begin(), v1.end());
			std::sort(v2.begin(), v2.end());
		}
		
		std::vector<T> d(v1.size());
		typename std::vector<T>::iterator it;
		it = std::set_difference( v1.begin(), v1.end(),
							 v2.begin(), v2.end(),
							 d.begin() );
		
		return std::vector<T>(d.begin(), d.begin()+int(it-d.begin()));
	}
	
}
#endif
