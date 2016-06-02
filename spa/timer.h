#ifndef __TIMER_H__
#define __TIMER_H__

#include <time.h>

inline double mytime (void)
{
    int flag;
    clockid_t cid = CLOCK_REALTIME; // CLOCK_MONOTONE might be better
    timespec tp;
    double timing;
	
    flag = clock_gettime(cid, &tp);
    if (flag == 0) timing = tp.tv_sec + 1.0e-9*tp.tv_nsec;
    else           timing = -17.0;         // If timer failed, return non-valid time
	
    return(timing);
}

inline void printLocalTime()
{
    time_t tim=time(NULL);
    tm *now=localtime(&tim);
    printf("[Time: %02d:%02d:%02d]\n", now->tm_hour, now->tm_min, now->tm_sec);
    printf("[Date: %d/%02d/%02d]\n", now->tm_year+1900, now->tm_mon+1, now->tm_mday);
}

#endif
