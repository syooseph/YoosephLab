#ifndef __LOGGER_H__
#define __LOGGER_H__

#include <string>
#include <cstdio>
#include <cstdlib>
#include "timer.h"

inline void log( std::string message, double ts )
{
    message = "[%11.2f] " + message + "\n";
    fprintf( stderr, message.c_str(), mytime()-ts );
}


#endif
