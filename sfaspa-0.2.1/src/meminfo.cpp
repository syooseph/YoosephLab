#include "meminfo.h"

size_t getTotalRAM()
{
	struct sysinfo info;
	sysinfo( &info );
	return (size_t)info.totalram * (size_t)info.mem_unit;
}

int parseLine(char* line)
{
	int i = strlen(line);
	while (*line < '0' || *line > '9') line++;
	line[i-3] = '\0';
	i = atoi(line);
	return i;
}
    

//Note: this value is in KB!
int getMemUsage()
{ 
	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];
    
	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmSize:", 7) == 0){
			result = parseLine(line);
			break;
		}
	}
	fclose(file);
	return (size_t)result;
}

void printMemoryUsage()
{
    printf("[Memory usage:% d KB]\n", getMemUsage());
}
