#include "galignment.h"


GlobalAlignPair::GlobalAlignPair() 
{
    // no nothing
}

GlobalAlignPair::GlobalAlignPair(TString &seq1, TString &seq2)
{
    align(seq1, seq2);
    setSummary();
}

// seq1:reference
// seq2:query
GlobalAlignPair::GlobalAlignPair(std::string &seq1, std::string &seq2)
{
    align(seq1, seq2);
    setSummary();
}


/*
http://trac.seqan.de/wiki/Tutorial/Alignments/AssignmentPairwiseGlobalAlignment2
Now we use an AlignConfig to configure our alignment to be semi-global (the second sequence being contained in the first sequence). The signature is AlignConfig<TTop, TLeft, TRight, TBottom>. TTop is true meaning that the first row of the DP matrix is initialized with zeros (gaps before the start of the second sequence are free), TLeft and TRight to false meaning that all gaps in the first sequence receive a penalty, and TBottom to true, leaving gaps at the end of the second sequence unpunished. 
*/
//TAlign 
void GlobalAlignPair::align(TString &seq1, TString &seq2) 
{
    using namespace seqan;
    // semi-global setup
    // 1.true: no penalty for leading gaps in second sequence 
    // 2. false: penalty for leading gaps in first sequence
    // 3. false: penalty for trailing gaps in first sequence
    // 3. true: no penalty for trailing gaps in second sequence
    AlignConfig<true,false,false,true> ac;    
    Blosum62 score_type(-1, -11);
		
    appendValue(rows(aln), seq1);
    appendValue(rows(aln), seq2);
    summary.score = globalAlignment(aln, score_type, ac, Gotoh());
    //return aln_pair;
    //std::cout << "Aligns Seq1[" << clippedBeginPosition(row(aln, 0)) << ":" << (clippedEndPosition(row(aln, 0))-1) << "]";
    //std::cout << " and Seq2[" << clippedBeginPosition(row(aln, 1)) << ":" <<  (clippedEndPosition(row(aln, 1))-1) << "]" << ::std::endl << ::std::endl;
}
	
//TAlign 
void GlobalAlignPair::align(std::string &seq1, std::string &seq2) 
{
    TString tstr1 = TString(seq1);
    TString tstr2 = TString(seq2);
    //return assembly::align(tstr1, tstr2);
    align(tstr1, tstr2);
}

void GlobalAlignPair::setSummary()
{
    using namespace seqan;

//     TString str1 = row(aln, 0);
//     TString str2 = row(aln, 1);
    
    setSequenceRange();
    setLeadingGaps();
    setEndGaps();

//     std::cout << "range\n";
//     std::cout << summary.s1se.first << "\t" << summary.s1se.second << "\n";
//     std::cout << summary.s2se.first << "\t" << summary.s2se.second << "\n";

    alignmentOuterBoundary();
    alignmentInnerRange();


    setGapPositions();
    //setAllGapPositions();
    setAlignmentStats();
}

void GlobalAlignPair::alignmentInnerRange()
{
    intPair ob = summary.outer;

    unsigned s = ob.first;
    while ( countGaps(row(aln,0), s) > 0 || countGaps(row(aln,1), s) > 0 ) s++;
    
    unsigned e = ob.second;
    while ( countGaps(row(aln,0), e) > 0 || countGaps(row(aln,1), e) > 0 ) 
        e--;

    summary.range = intPair(s, e);
}

void GlobalAlignPair::alignmentOuterBoundary()
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
    //( e1 < e2 ) ? ( boundary.second = e1 ) : ( boundary.second = e2 );

    // possible seqan bug??
    ( e1 < e2 ) ? ( boundary.second = e2 ) : ( boundary.second = e1 );
    // this puts 'A' in ending gap positions
    
    summary.outer = boundary;
}


void GlobalAlignPair::setSequenceRange()
{
    int s1 = beginPosition(row(aln, 0));
    int s2 = beginPosition(row(aln, 1));
    int e1 = endPosition(row(aln, 0))-1;
    int e2 = endPosition(row(aln, 1))-1;
    
    //std::cout << "s1:" << s1 << ",e1:" << e1 << "\ts2:" << s2 << ",e2:" << e2 << "\n";
    summary.s1se = intPair(s1,e1);
    summary.s2se = intPair(s2,e2);
}

// debugging alignment outrange purpose
void GlobalAlignPair::printGaps()
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

void GlobalAlignPair::setLeadingGaps()
{
    int p1 = summary.s1se.first;
    int p2 = summary.s2se.first;

    if ( p1 == p2 ) 
        summary.lgap = intPair(0,0);
    else if ( p1 > p2 )
        summary.lgap = intPair( p1-p2, 0);
    else
        summary.lgap = intPair( 0, p2-p1 );

    //std::cout << "leading gaps:" << summary.lgap.first << "\t" << summary.lgap.second << "\n";
}

void GlobalAlignPair::setEndGaps()
{
    int p1 = summary.s1se.second;
    int p2 = summary.s2se.second;

    if ( p1 == p2 ) 
        summary.egap = intPair(0,0);
    else if ( p1 > p2 )
        summary.egap = intPair( 0, p1-p2);
    else
        summary.egap = intPair( p2-p1, 0 );
    //std::cout << "end gaps:" << summary.egap.first << "\t" << summary.egap.second << "\n";
}

void GlobalAlignPair::setAlignmentStats()
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

void GlobalAlignPair::setGapPositions()
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

// void GlobalAlignPair::setAllGapPositions()
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

void GlobalAlignPair::printAlignment(std::ostream &out)
{
    out << aln;
}
