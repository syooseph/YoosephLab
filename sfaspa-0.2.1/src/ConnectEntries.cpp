#include "ConnectEntries.h"

JoinEntry::JoinEntry( PathId p, int n, int s )
{
    pid = p; nkmers = n; seqlen = s;
}

//////////////////////////////////////////////////////////////////////
// ReadEntry
//////////////////////////////////////////////////////////////////////
ReadEntry::ReadEntry()
{

}

ReadEntry::ReadEntry( ReadId i, unsigned r, unsigned p ) 
{
    read = i; rpos = r; ppos = p;
}

void ReadEntry::dump( std::fstream &out )
{
    out.write( (char*)&read, sizeof(ReadId) );
    out.write( (char*)&rpos, sizeof(unsigned) );
    out.write( (char*)&ppos, sizeof(unsigned) );
}

void ReadEntry::load( std::fstream &in )
{
    in.read( (char*)&read, sizeof(ReadId) );
    in.read( (char*)&rpos, sizeof(unsigned) );
    in.read( (char*)&ppos, sizeof(unsigned) );
}

