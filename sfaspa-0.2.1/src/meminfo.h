// See also http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process

#include <sys/sysinfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

size_t getTotalRAM();

int parseLine(char* line);

int getMemUsage();

void printMemoryUsage();
