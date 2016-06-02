#include "extender.h"

Extender::Extender()
{
}

Extender::Extender(std::string s) 
{
    sequence = s;
}

Extender::~Extender()
{
}

void Extender::init()
{
	sequence = "";
	members = PathIdSet();
	positions = StartPosMap();
    aligns = PathAlignMap();
    lstop = rstop = -1;
    ltrim = rtrim = -1;
}

void Extender::dump( std::fstream &out ) 
{
    size_t len = sequence.size();
    const char *cstr = sequence.c_str();
    out.write((char*)&len, sizeof(size_t));
    out.write(cstr, len);
	
    size_t ctmem = members.size();
    out.write((char*)&ctmem, sizeof(size_t));
    for ( PathIdSet::iterator it = members.begin(); it != members.end(); ++it ) {
        PathId pid = *it;
        out.write((char*)&(pid), sizeof(PathId));
    }

    size_t ctpos = positions.size();
    out.write((char*)&ctpos, sizeof(size_t));
    for ( StartPosMap::iterator it = positions.begin(); it != positions.end(); ++it ) {
        PathId pid = it->first;
        int    pos = it->second;
        out.write((char*)&pid, sizeof(PathId));
        out.write((char*)&pos, sizeof(int));
    }
	
    size_t ctaln = aligns.size();
    out.write((char*)&ctaln, sizeof(size_t));
    for ( PathAlignMap::iterator it = aligns.begin(); it != aligns.end(); ++it ) {
        PathId pid = it->first;
        out.write((char*)&pid, sizeof(PathId));
        (*it).second.dump(out);
    }

    out.write((char*)&lstop, sizeof(int));
    out.write((char*)&rstop, sizeof(int));

    out.write((char*)&ltrim, sizeof(int));
    out.write((char*)&rtrim, sizeof(int));
}

void Extender::load( std::fstream &in ) 
{
    size_t len;
    in.read((char*)&len, sizeof(size_t));
    assert(len>0);
    char *seq = new char[len+1];
    in.read(seq, len);
    seq[len]='\0';
    sequence = std::string(seq);
    delete[] seq;

    size_t ctmem;
    in.read((char*)&ctmem, sizeof(size_t));
    for ( size_t i = 0; i < ctmem; i++ ) {
        PathId pid;
        in.read((char*)&pid, sizeof(PathId));
        members.insert(pid);
    }

    size_t ctpos;
    in.read((char*)&ctpos, sizeof(size_t));
    for ( size_t i = 0; i < ctpos; i++ ) {
        PathId pid; in.read((char*)&pid, sizeof(PathId));
        int    pos; in.read((char*)&pos, sizeof(int));
        positions.insert(std::pair<PathId, int>(pid, pos));
    }
	
    size_t ctaln;
    in.read((char*)&ctaln, sizeof(size_t));
    for ( size_t i = 0; i < ctaln; i++ ) {
        PathId pid; in.read((char*)&pid, sizeof(PathId));
        AlignSummary sum; sum.load(in);
        aligns[pid] = sum;
    }

    in.read((char*)&lstop, sizeof(int));
    in.read((char*)&rstop, sizeof(int));

    in.read((char*)&ltrim, sizeof(int));
    in.read((char*)&rtrim, sizeof(int));
}
