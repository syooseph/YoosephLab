#include "alignsummary.h"

AlignSummary::AlignSummary() 
{
    init();
}

AlignSummary::AlignSummary(const AlignSummary &source) 
{
    score    = source.score;
    length   = source.length;
    match    = source.match;
    mismatch = source.mismatch;
    positive = source.positive;
    posrate  = source.posrate;
    range    = source.range;
    outer    = source.outer;
    lgap     = source.lgap;
    egap     = source.egap;
    s1se     = source.s1se;
    s2se     = source.s2se;
    ilist    = source.ilist;
    dlist    = source.dlist;
}

AlignSummary& AlignSummary::operator= (const AlignSummary &source)
{
    if (this == &source) return *this;
    
    score    = source.score;
    length   = source.length;
    match    = source.match;
    mismatch = source.mismatch;
    positive = source.positive;
    posrate  = source.posrate;
    range    = source.range;
    outer    = source.outer;
    lgap     = source.lgap;
    egap     = source.egap;
    s1se     = source.s1se;
    s2se     = source.s2se;
    ilist    = source.ilist;
    dlist    = source.dlist;
    
    return *this;
}

void AlignSummary::init()
{
    score = length = match = mismatch = positive = 0;
    posrate = 0.0;
    lgap = egap = s1se = s2se = range = outer = std::pair<int,int>(0,0);
    ilist.clear(); dlist.clear();
}
	   
void AlignSummary::print(std::ostream &out)
{
    out << "\nAlignment Summmary\n";
    out << "length:" << length << "\n";
    out << "match:" << match << "\n";
    out << "mismatch:" << mismatch << "\n";
    out << "positive:" << positive << "\n";
    out << "score:" << score << "\n";
    out << "%positive:" << posrate << "\n";
    out << "insertion:" << ilist.size() << "\n";
    out << "deletion:" << dlist.size() << "\n";
    out << "sbjct range:" << s1se.first << "\t" << s1se.second << "\n";
    out << "query range:" << s2se.first << "\t" << s2se.second << "\n";
    out << "leading gap:" << lgap.first << "\t" << lgap.second << "\n";
    out << "trailing gap:" << egap.first << "\t" << egap.second << "\n";
    out << "outer range:" << outer.first << "\t" << outer.second << "\n";
    out << "inner range:" << range.first << "\t" << range.second << "\n";
    out << "\n";
}

void AlignSummary::printGap(std::ostream &out)
{
    AlignPosList::iterator it;    
    for ( it = ilist.begin(); it != ilist.end(); ++it ) {
        out << "Insertion:\t";
        it->print(out);
    }
    for ( it = dlist.begin(); it != dlist.end(); ++it ) {
        out << "Deletion:\t";
        it->print(out);
    }
        
}

void AlignSummary::dump( std::fstream &out ) 
{
    out.write( (char*)&score, sizeof(int) );
    out.write( (char*)&length, sizeof(int) );
    out.write( (char*)&match, sizeof(int) );
    out.write( (char*)&mismatch, sizeof(int) );
    out.write( (char*)&positive, sizeof(int) );
    out.write( (char*)&posrate, sizeof(double) );		
    out.write( (char*)&range.first, sizeof(int) );		
    out.write( (char*)&range.second, sizeof(int) );			
    out.write( (char*)&outer.first, sizeof(int) );		
    out.write( (char*)&outer.second, sizeof(int) );		
    out.write( (char*)&lgap.first, sizeof(int) );		
    out.write( (char*)&lgap.second, sizeof(int) );		
    out.write( (char*)&egap.first, sizeof(int) );		
    out.write( (char*)&egap.second, sizeof(int) );		
    out.write( (char*)&s1se.first, sizeof(int) );		
    out.write( (char*)&s1se.second, sizeof(int) );		
    out.write( (char*)&s2se.first, sizeof(int) );		
    out.write( (char*)&s2se.second, sizeof(int) );		
    size_t npos = ilist.size();
    out.write( (char*)&npos, sizeof(size_t) );		
    for ( AlignPosList::iterator it = ilist.begin(); it != ilist.end(); ++it )
        it->dump(out);
    npos = dlist.size();
    out.write( (char*)&npos, sizeof(size_t) );		
    for ( AlignPosList::iterator it = dlist.begin(); it != dlist.end(); ++it )
        it->dump(out);
}

void AlignSummary::load( std::fstream &in ) 
{
    in.read( (char*)&score, sizeof(int) );
    in.read( (char*)&length, sizeof(int) );
    in.read( (char*)&match, sizeof(int) );
    in.read( (char*)&mismatch, sizeof(int) );
    in.read( (char*)&positive, sizeof(int) );
    in.read( (char*)&posrate, sizeof(double) );		
    in.read( (char*)&range.first, sizeof(int) );		
    in.read( (char*)&range.second, sizeof(int) );			
    in.read( (char*)&outer.first, sizeof(int) );		
    in.read( (char*)&outer.second, sizeof(int) );		
    in.read( (char*)&lgap.first, sizeof(int) );		
    in.read( (char*)&lgap.second, sizeof(int) );		
    in.read( (char*)&egap.first, sizeof(int) );		
    in.read( (char*)&egap.second, sizeof(int) );		
    in.read( (char*)&s1se.first, sizeof(int) );		
    in.read( (char*)&s1se.second, sizeof(int) );		
    in.read( (char*)&s2se.first, sizeof(int) );		
    in.read( (char*)&s2se.second, sizeof(int) );		
    size_t npos;
    in.read( (char*)&npos, sizeof(size_t) );		
    for ( size_t i = 0; i < npos; i++ ) {
        AlignIndex aln;
        aln.load(in);
        ilist.push_back(aln);
    }
    in.read( (char*)&npos, sizeof(size_t) );		
    for ( size_t i = 0; i < npos; i++ ) {
        AlignIndex aln;
        aln.load(in);
        dlist.push_back(aln);
    }
}

//====================================================================
// Shift indel positions of query/sbjct sequence.
// ref_seq: (1:sbjct and 0:query)
//====================================================================
void AlignSummary::shift( int offset, bool ref_seq )
{
    for ( AlignPosList::iterator it = ilist.begin(); it != ilist.end(); ++it )
        ref_seq ? it->ref_pos += offset : it->qry_pos += offset;

    for ( AlignPosList::iterator it = dlist.begin(); it != dlist.end(); ++it )
        ref_seq ? it->ref_pos += offset : it->qry_pos += offset;
}


void AlignSummary::self( size_t n )
{
    score    = 1;
    length   = n;
    match    = n;
    mismatch = 0;
    positive = n;
    posrate  = 1;
    range    = std::pair<int,int>(0,n-1);
    outer    = std::pair<int,int>(0,n-1);
}
