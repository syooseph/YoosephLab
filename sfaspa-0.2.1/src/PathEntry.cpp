#include "PathEntry.h"

PathEntry::PathEntry() 
{

}

PathEntry::PathEntry( std::string &s, int ls, int rs, int lt, int rt )
{
    seq = s;
    lstop = ls; rstop = rs;
    ltrim = lt; rtrim = rt;
}

void PathEntry::dump( std::fstream &out )
{
    size_t nchar = seq.size();
    out.write((char*)&nchar, sizeof(size_t));
    out.write(seq.c_str(), nchar);
    
    out.write((char*)&lstop, sizeof(int));
    out.write((char*)&rstop, sizeof(int));
    out.write((char*)&ltrim, sizeof(int));
    out.write((char*)&rtrim, sizeof(int));
}

void PathEntry::load( std::fstream &in )
{
    size_t len;
    in.read((char*)&len, sizeof(size_t));
    assert(len>0);
    char *str = new char[len+1];
    in.read(str, len);
    str[len]='\0';
    seq = std::string(str);
    delete[] str;

    in.read((char*)&lstop, sizeof(int));
    in.read((char*)&rstop, sizeof(int));
    in.read((char*)&ltrim, sizeof(int));
    in.read((char*)&rtrim, sizeof(int));
}

