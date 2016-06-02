#include "lalignment.h"


LocalAlignPair::LocalAlignPair() 
{
    // no nothing
}

LocalAlignPair::LocalAlignPair(TString &seq1, TString &seq2)
{
    align(seq1, seq2);
    setSummary();
}

// seq1:reference
// seq2:query
LocalAlignPair::LocalAlignPair(std::string &seq1, std::string &seq2)
{
    align(seq1, seq2);
    setSummary();
}

LocalAlignPair::LocalAlignPair(std::string &sbjct, std::string &query, int gex, int gop)
{
    setGapPenalty(gex, gop);
    align(sbjct, query, anchor);
    setSummary();
}


void LocalAlignPair::align(TString &seq1, TString &seq2) 
{
    using namespace seqan;
    Blosum62 scoring(gapext, gapopen);

    TStringSet sequences;
    appendValue(sequences, seq1);
    appendValue(sequences, seq2);

    TAlignGraph g(sequences);
    localAlignment(g, scoring, SmithWaterman());
    alignG = g;
}
	
void LocalAlignPair::align(std::string &seq1, std::string &seq2) 
{
    TString tstr1 = TString(seq1);
    TString tstr2 = TString(seq2);
    align(tstr1, tstr2);
}

void LocalAlignPair::setSummary()
{
    using namespace seqan;

    std::string matrix;
    convertAlignment( alignG, matrix );

    std::string seq1 = matrix.substr(0, matrix.size()/2);
    std::string seq2 = matrix.substr(matrix.size()/2, matrix.size()/2);

    setSequenceRange(seq1, seq2);
    setLeadingGaps();
    setEndGaps();

    alignmentOuterBoundary();
    alignmentInnerRange(seq1, seq2);


    setGapPositions();
    setAlignmentStats();
}

void LocalAlignPair::alignmentInnerRange()
{
    intPair ob = summary.outer;

    unsigned s = ob.first;
    while ( countGaps(row(aln,0), s) > 0 || countGaps(row(aln,1), s) > 0 ) s++;
    
    unsigned e = ob.second;
    while ( countGaps(row(aln,0), e) > 0 || countGaps(row(aln,1), e) > 0 ) 
        e--;

    summary.range = intPair(s, e);
}

void LocalAlignPair::alignmentOuterBoundary()
{
//     int s1 = beginPosition(row(aln, 0));
//     int s2 = beginPosition(row(aln, 1));
//     int e1 = endPosition(row(aln, 0))-1;
//     int e2 = endPosition(row(aln, 1))-1;

    int s1 = summary.s1se.first;
    int s2 = summary.s2se.first;
    int e1 = summary.s1se.second;
    int e2 = summary.s2se.second;
    
    intPair boundary;
    
    ( s1 < s2 ) ? ( boundary.first = s1 ) : ( boundary.first = s2 );
    ( e1 < e2 ) ? ( boundary.second = e1 ) : ( boundary.second = e2 );

    // possible seqan bug??
    //( e1 < e2 ) ? ( boundary.second = e2 ) : ( boundary.second = e1 );
    // this puts 'A' in ending gap positions
    
    summary.outer = boundary;
}


void LocalAlignPair::setSequenceRange( std::string &seq1, std::string &seq2 )
{
    int s = 0;
    int e = (int)seq1.size()-1;
    while ( seq1[s] == '-' && seq2[s] == '-' ) s++;
    while ( seq1[e] == '-' && seq2[e] == '-' ) e--;

    int g1 = 0;
    int g2 = 0;
    for ( int i = 0; i <= s; i++ ) {
        if ( seq1[i] == '-' ) g1++;
        if ( seq2[i] == '-' ) g2++;
    }
    int s1 = s-g1;
    int s2 = s-g2;

    int e1 = e;
    int e2 = e;
    for ( int i = e; i > s; i-- ) {
        if ( seq1[i] == '-' ) e1--;
        if ( seq2[i] == '-' ) e2--;
    }
    
    summary.s1se = intPair(s1,e1);
    summary.s2se = intPair(s2,e2);
}

// debugging alignment outrange purpose
void LocalAlignPair::printGaps()
{
//     int s1 = beginPosition(row(aln, 0));
//     int s2 = beginPosition(row(aln, 1));
//     int e1 = endPosition(row(aln, 0))-1;
//     int e2 = endPosition(row(aln, 1))-1;

    int s1 = summary.s1se.first;
    int s2 = summary.s2se.first;
    int e1 = summary.s1se.second;
    int e2 = summary.s2se.second;
    std::cout << "gaps from:" << s1 << " to " << e1 << "\n";
    for ( int i = s1; i <= e1; i++ )
        std::cout << countGaps(row(aln,0), i) << " ";
    std::cout << "\n";
    for ( int i = s2; i <= e2; i++ )
        std::cout << countGaps(row(aln,1), i) << " ";
    std::cout << "\n";


}

void LocalAlignPair::setLeadingGaps( )
{
    int p1 = summary.s1se.first;
    int p2 = summary.s2se.first;

    if ( p1 == p2 ) 
        summary.lgap = intPair(0,0);
    else if ( p1 > p2 )
        summary.lgap = intPair( p1-p2, 0);
    else
        summary.lgap = intPair( 0, p2-p1 );
}

void LocalAlignPair::setEndGaps()
{
    int p1 = summary.s1se.second;
    int p2 = summary.s2se.second;

    if ( p1 == p2 ) 
        summary.egap = intPair(0,0);
    else if ( p1 > p2 )
        summary.egap = intPair( 0, p1-p2);
    else
        summary.egap = intPair( p2-p1, 0 );
}

void LocalAlignPair::setAlignmentStats()
{
    for ( int i = summary.range.first; i <= summary.range.second; i++ ) {
//         if ( countGaps(row(aln,0),i) > 0 ) summary.del.push_back(i); //summary.gap++;
//         else if ( countGaps(row(aln,1),i) > 0 ) summary.ins.push_back(i); //summary.gap++;
        if ( countGaps(row(aln,0),i) > 0 ) summary.ins.push_back(i); //summary.gap++;
        else if ( countGaps(row(aln,1),i) > 0 ) summary.del.push_back(i); //summary.gap++;
        else {
            char c1 = (char)getValue(row(aln, 0), i);
            char c2 = (char)getValue(row(aln, 1), i);
            if ( c1 == c2 ) { summary.match++; summary.positive++; }
            else { 
                summary.mismatch++;
                if ( scoring::getScore(c1, c2, BLOSUM62) > 0 ) summary.positive++;
            }
        }
    }
    summary.gap = summary.ins.size() + summary.del.size();
    summary.length = summary.range.second - summary.range.first + 1;//e-s+1;
    summary.posrate = (double)summary.positive/summary.length;
}

void LocalAlignPair::setGapPositions()
{
    int cins = 0;
    int cdel = 0;
    for ( int i = summary.range.first; i <= summary.range.second; i++ ) {
        if ( countGaps(row(aln,0),i) > 0 ) {
//             summary.dlist.push_back( AlignIndex(i, i-cdel-summary.lgap.first, i-cins-summary.lgap.second) );
//             cdel++;
            //summary.ilist.push_back( AlignIndex(i, i-cins-summary.lgap.first, i-cdel-summary.lgap.second) );
            summary.ilist.push_back( AlignIndex(i, i-cdel-summary.lgap.second, i-cins-summary.lgap.first) );
            cins++;
        }
        else if ( countGaps(row(aln,1),i) > 0 ) {
//             summary.ilist.push_back( AlignIndex(i, i-cins-summary.lgap.second, i-cdel-summary.lgap.first) );
//             cins++;
            summary.dlist.push_back( AlignIndex(i, i-cdel-summary.lgap.second, i-cins-summary.lgap.first) );
            //summary.dlist.push_back( AlignIndex(i, i-cins-summary.lgap.first, i-cdel-summary.lgap.second) );
            cdel++;
        }
    }
}

// void LocalAlignPair::setAllGapPositions()
// {
//     int cins = 0;
//     int cdel = 0;
//     for ( int i = summary.outer.first; i <= summary.outer.second; i++ ) {
//         if ( countGaps(row(aln,0),i) > 0 ) {
//             std::cout << "row:0\t" << i << "\t" << i-cdel-summary.lgap.first << "\t" << i-cins-summary.lgap.second << "\n";
//             summary.alldels.push_back( AlignIndex(i, i-cdel-summary.lgap.first, i-cins-summary.lgap.second) );
//             //summary.alldels.push_back( AlignIndex(i, i-cdel, i-cins) );
//             cdel++;
//         }
//         else if ( countGaps(row(aln,1),i) > 0 ) {
//             std::cout << "row:1\t" << i << "\t" << i-cdel-summary.lgap.second << "\t" << i-cins-summary.lgap.first << "\n";
//             //summary.allinss.push_back( AlignIndex(i, i-cins, i-cdel) );
//             summary.allinss.push_back( AlignIndex(i, i-cins-summary.lgap.second, i-cdel-summary.lgap.first) );
//             cins++;
//         }
//     }
// }

void LocalAlignPair::printAlignment(std::ostream &out)
{
    out << aln;
}
