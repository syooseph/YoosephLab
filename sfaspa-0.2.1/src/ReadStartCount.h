#ifndef __READSTARTCOUNT_H
#define __READSTARTCOUNT_H

#include "gsa.h"

class ReadStartCount
{
 private:
	unsigned  size;
	unsigned* cumulatives;



 public:
	ReadStartCount();
	ReadStartCount( const GSA *gsa, const unsigned &size );
	void init( const GSA *gsa, const unsigned &size);
	~ReadStartCount();

	unsigned getCount(int a, int b);
	unsigned get(int a);
};


#endif
