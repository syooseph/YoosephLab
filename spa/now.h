#ifndef __NOW_H__
#define __NOW_H__

#include <cstdio>
#include <ctime>

inline void now(char buffer[], size_t N) 
{
    time_t rawtime;
    struct tm * timeinfo;

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    strftime (buffer,N,"%Y.%m.%d.%H.%M.%S",timeinfo);
}

inline void stamp(char buffer[], size_t N) 
{
    time_t rawtime;
    struct tm * timeinfo;

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    strftime (buffer,N,"%Y%m%d%H%M%S",timeinfo);
}

inline std::string timeStamp() 
{
	char buffer[256];
	sprintf(buffer, "%f", mytime());
	
	size_t N = 80;
	char curr[N];
	stamp(curr, N);
	return std::string(curr);
}

#endif
