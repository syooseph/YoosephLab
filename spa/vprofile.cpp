#include "vprofile.h"

ProfileVector::ProfileVector()
{
    size = 0;
    entry = NULL;
}

ProfileVector::ProfileVector(Profile &profile)
{
    build(profile);
}

ProfileVector::ProfileVector(const ProfileVector &source)
{
    __copy(source);
}

ProfileVector& ProfileVector::operator= (const ProfileVector &source)
{
    __copy(source);
    return *this;
}

void ProfileVector::__copy(const ProfileVector &source)
{
    size = source.size;
    if ( size ) entry = new pEntry[size];
    for ( size_t i = 0; i < size; i++ ) 
        entry[i] = source.entry[i];
}

ProfileVector::~ProfileVector()
{
    clear();
}

void ProfileVector::build(Profile &profile)
{
    std::list<pEntry> tmplist;
    for ( size_t i = 0; i < profile.ncol; i++ )
        for ( size_t j = 0; j < (size_t)naas; j++ ) 
            if ( profile.matrix[i][j] > 0 ) 
                tmplist.push_back( pEntry(i,j,profile.matrix[i][j]) );

//     //entry = std::vector<pEntry>(tmplist.begin(), tmplist.end());
//     entry = std::vector<pEntry>(tmplist.size());
//     int i = 0;
//     for ( std::list<pEntry>::iterator it = tmplist.begin(); it != tmplist.end(); ++it ) {
//         entry[i] = *it; 
//         i++;
//     }
    
    size = tmplist.size();
    if ( size == 0 ) return;

    entry = new pEntry[size];
    
    size_t i = 0;
    for ( std::list<pEntry>::iterator it = tmplist.begin(); it != tmplist.end(); ++it ) {
        entry[i] = *it; 
        i++;
    }
}

//void ProfileVector::convert(Profile &profile)
Profile ProfileVector::convert()
{
    Profile profile;
//     size_t size = entry.size();
    if ( size == 0 ) return profile;


    int npos = entry[size-1].pos + 1;
    // I don't think this is necessary
    //     if ( npos == 0 ) return;

    //std::cout << "Npos:" << npos << "\n";

//     if (npos <=0 || npos > size) 
//         std::cout << "[Warning] profile size is smaller than the position of last column (size:" << size << "\tpos:" << npos << ")\n";

    //assert(npos>0 && npos <= size);
    assert(npos>0);// && npos <= size);
    profile.init(npos);
    for ( unsigned i = 0; i < size; i++ ) 
        profile.matrix[ entry[i].pos ][ entry[i].col ] = entry[i].num;
    return profile;
}

void ProfileVector::clear()
{
    if ( size ) {
        delete[] entry;
        size = 0;
    }
}

	
void ProfileVector::dump(std::ostream &out)
{
    out.write((char*)&size, sizeof(unsigned));
    for ( unsigned i = 0; i < size; i++ ) 
        entry[i].dump(out);
}

void ProfileVector::load(std::istream &in) 
{
    in.read((char*)&size, sizeof(unsigned));
    entry = new pEntry[size];
    for ( unsigned i = 0; i < size; i++ ) 
        entry[i].load(in);
}

