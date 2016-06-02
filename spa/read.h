#ifndef __READ_H__
#define __READ_H__

typedef unsigned ReadId;
typedef unsigned CoverageType;

class Read 
{
 public:
	CoverageType size;
	ReadId *rid;
	
	Read() { size = 0; }
	Read(CoverageType s) { rid = new ReadId[s]; size = 0; }
	~Read() { delete[] rid; }
	void add(ReadId &r) {rid[size++] = r; } 
};

#endif
