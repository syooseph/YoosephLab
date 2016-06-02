#include "rplace.h"

void ReadPlacement::print( std::ostream &out )
{
    out << "rid:" << rid << "\t"
        << "rp:" << read_pos << "\t"
        << "pp:" << path_pos << "\t"
        << "len:" << length << "\t"
        << "#ins:" << ilist.size() << "\t"
        << "#del:" << dlist.size() << "\n";
}

void ReadPlacement::dump( std::fstream &out )
{
    out.write( (char*)&rid, sizeof(ReadId) );
    out.write( (char*)&read_pos, sizeof(int) );
    out.write( (char*)&path_pos, sizeof(int) );
    out.write( (char*)&length, sizeof(int) );
    size_t count = ilist.size();
    out.write( (char*)&count, sizeof(size_t) );
    for ( std::list<Mismatch>::iterator it = ilist.begin(); it != ilist.end(); ++it )
        it->dump(out);
    count = dlist.size();
    out.write( (char*)&count, sizeof(size_t) );
    for ( std::list<Mismatch>::iterator it = dlist.begin(); it != dlist.end(); ++it )
        it->dump(out);
}

void ReadPlacement::load( std::fstream &in )
{
    in.read( (char*)&rid, sizeof(ReadId) );
    in.read( (char*)&read_pos, sizeof(int) );
    in.read( (char*)&path_pos, sizeof(int) );
    in.read( (char*)&length, sizeof(int) );
    size_t count;
    in.read( (char*)&count, sizeof(size_t) );
    for ( size_t i = 0; i < count; i++ ) {
        Mismatch mis; mis.load(in);
        ilist.push_back(mis);
    }
    in.read( (char*)&count, sizeof(size_t) );
    for ( size_t i = 0; i < count; i++ ) {
        Mismatch mis; mis.load(in);
        dlist.push_back(mis);
    }
}

