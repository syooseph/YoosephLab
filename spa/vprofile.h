#ifndef __PROFILE_VECTOR_H__
#define __PROFILE_VECTOR_H__

#include <vector>
#include <list>
#include "profile.h"

struct pEntry
{
	unsigned pos; // position in consensus
	unsigned col; // amino acid
	unsigned num; // count

	pEntry() {}
	pEntry(const pEntry &source) 
	{
		pos = source.pos; col = source.col; num = source.num;
	}
	pEntry& operator= (const pEntry &source) 
	{
		pos = source.pos; col = source.col; num = source.num;
		return *this;
	}
	pEntry(unsigned r, char c, unsigned n) 
	{
		pos = r; col = c; num = n;
	}
	void dump(std::ostream &out) 
	{
		out.write((char*)&pos, sizeof(unsigned));
		out.write((char*)&col, sizeof(unsigned));
		out.write((char*)&num, sizeof(unsigned));
	}
	void load(std::istream &in) 
	{
		in.read((char*)&pos, sizeof(unsigned));
		in.read((char*)&col, sizeof(unsigned));
		in.read((char*)&num, sizeof(unsigned));
	}

};

class ProfileVector
{
 private:
	unsigned  size;
	pEntry *entry;
	//std::vector<pEntry> entry;
 public:
	ProfileVector();
	ProfileVector(Profile &profile);
	ProfileVector(const ProfileVector &source);
	ProfileVector& operator = (const ProfileVector &source);
	~ProfileVector();
	//void convert(Profile &profile);
	Profile convert();
	void build(Profile &profile);
	void dump(std::ostream &out);
	void load(std::istream &in);
	void clear();
	unsigned getSize() { return size; }
 private:
	void __copy(const ProfileVector &source);
};




#endif
