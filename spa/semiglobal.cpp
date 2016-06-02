#include "semiglobal.h"


SemiGlobalAlign::SemiGlobalAlign() 
{
    // do nothing
}

SemiGlobalAlign::SemiGlobalAlign(TString &sbjct, TString &query, int anchor)
{
    align(sbjct, query, anchor);
    setSummary();
    //if ( !adjustHead() ) summary.posrate = summary.score = 0;
    refine();
}


SemiGlobalAlign::SemiGlobalAlign(std::string &sbjct, std::string &query, int anchor)
{
    align(sbjct, query, anchor);
    setSummary();
    //if ( !adjustHead() ) summary.posrate = summary.score = 0;
    refine();
}

SemiGlobalAlign::SemiGlobalAlign(std::string &sbjct, std::string &query, int anchor, int gex, int gop)
{
    setGapPenalty(gex, gop);
    align(sbjct, query, anchor);
    setSummary();
    //if ( !adjustHead() ) summary.posrate = summary.score = 0;
    refine();
}


/*
http://trac.seqan.de/wiki/Tutorial/Alignments/AssignmentPairwiseGlobalAlignment2
Now we use an AlignConfig to configure our alignment to be semi-global (the second sequence being contained in the first sequence). The signature is AlignConfig<TTop, TLeft, TRight, TBottom>. TTop is true meaning that the first row of the DP matrix is initialized with zeros (gaps before the start of the second sequence are free), TLeft and TRight to false meaning that all gaps in the first sequence receive a penalty, and TBottom to true, leaving gaps at the end of the second sequence unpunished. 
*/
//TAlign 
void SemiGlobalAlign::align(TString &sbjct, TString &query, int anchor) 
{
    using namespace seqan;
    
    Blosum62 score_type(gapext, gapopen);
    appendValue(rows(aln), sbjct);
    appendValue(rows(aln), query);
    
    if ( anchor == ANCHOR_CENTER ) {
        // semi-global setup
        // 1.true: no penalty for leading gaps in second sequence 
        // 2. false: penalty for leading gaps in first sequence
        // 3. false: penalty for trailing gaps in first sequence
        // 3. true: no penalty for trailing gaps in second sequence
        AlignConfig<true,false,false,true> ac;    
        summary.score = globalAlignment(aln, score_type, ac, Gotoh());
    } 
    else if ( anchor == ANCHOR_RIGHT ) {
        AlignConfig<true,false,true,false> ac;    
        summary.score = globalAlignment(aln, score_type, ac, Gotoh());
        
    } 
    else if ( anchor == ANCHOR_LEFT ) {
        AlignConfig<false,true,false,true> ac;    
        summary.score = globalAlignment(aln, score_type, ac, Gotoh());
    }
}
	
void SemiGlobalAlign::align(std::string &sbjct, std::string &query, int anchor) 
{
    TString tstr1 = TString(sbjct);
    TString tstr2 = TString(query);
    //return assembly::align(tstr1, tstr2);
    align(tstr1, tstr2, anchor);
}

void SemiGlobalAlign::setSummary()
{
    using namespace seqan;

    setSequenceRange();
    setLeadingGaps();
    setEndGaps();

    alignmentOuterBoundary();
    alignmentInnerRange();


    setGapPositions();
    //setAllGapPositions();
    setAlignmentStats();
}

void SemiGlobalAlign::alignmentInnerRange()
{
    intPair ob = summary.outer;

    unsigned s = ob.first;
    while ( countGaps(row(aln,0), s) > 0 || countGaps(row(aln,1), s) > 0 ) s++;
    
    unsigned e = ob.second;
    while ( countGaps(row(aln,0), e) > 0 || countGaps(row(aln,1), e) > 0 ) 
        e--;

    summary.range = intPair(s, e);
}

void SemiGlobalAlign::alignmentOuterBoundary()
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


void SemiGlobalAlign::setSequenceRange()
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
void SemiGlobalAlign::printGaps()
{
//     int s1 = beginPosition(row(aln, 0));
//     int s2 = beginPosition(row(aln, 1));
//     int e1 = endPosition(row(aln, 0))-1;
//     int e2 = endPosition(row(aln, 1))-1;

//     int s1 = summary.s1se.first;
//     int s2 = summary.s2se.first;
//     int e1 = summary.s1se.second;
//     int e2 = summary.s2se.second;
//     std::cout << "gaps from:" << s1 << " to " << e1 << "\n";
//     for ( int i = s1; i <= e1; i++ ) {
//         int ngap = countGaps(row(aln,0), i);
//         if ( ngap > 0 ) std::cout << i << ":" << ngap << " ";
//     }
//     std::cout << "\n";

//     for ( int i = s2; i <= e2; i++ ) {
//         int ngap = countGaps(row(aln,1), i);
//         if ( ngap > 0 ) std::cout << i << ":" << ngap << " ";
//     }
//     std::cout << "\n";

    std::list<AlignIndex>::iterator it;
    std::cout << "Insertions:" << "\n";
    for ( it = summary.ilist.begin(); it != summary.ilist.end(); ++it )
        std::cout << it->aln_pos << ":" << it->seq_pos << ":" << it->ref_pos << " ";
    std::cout << "\n";

    std::cout << "Deletions:" << "\n";
    for ( it = summary.dlist.begin(); it != summary.dlist.end(); ++it )
        std::cout << it->aln_pos << ":" << it->seq_pos << ":" << it->ref_pos << " ";
    std::cout << "\n";
}

void SemiGlobalAlign::setLeadingGaps()
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

void SemiGlobalAlign::setEndGaps()
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

void SemiGlobalAlign::setAlignmentStats()
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

void SemiGlobalAlign::setGapPositions()
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

// void SemiGlobalAlign::setAllGapPositions()
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


void SemiGlobalAlign::refine()
{
    int nshift = refineHead();
    refineTail(nshift);
    //printGaps();
}

void SemiGlobalAlign::printAlignment(std::ostream &out)
{
    out << aln;
}

/*
0     .    :    .    :    .    :    .    :    .    :
        KP--MSKFLDRFRYFKQKGETFADGHGQLLNTNRDWEDGYRQRWQHDKIV
            ||| ||||||||||||||||||||    |||||| |||||| ||||
        --HVMSKLLDRFRYFKQKGETFADGHGQVMHSNRDWEDSYRQRWQFDKIV
*/

void SemiGlobalAlign::__trimIndels( int s, int e )
{
    for (AlignPosList::iterator it = summary.dlist.begin(); it != summary.dlist.end();  ) {
        if ( it->aln_pos < s || it->aln_pos >= e ) ++it;
        else summary.dlist.erase(it++);
    }
    for (AlignPosList::iterator it = summary.ilist.begin(); it != summary.ilist.end();  ) {
        if ( it->aln_pos < s || it->aln_pos >= e ) ++it;
        else summary.ilist.erase(it++);
    }
}

void SemiGlobalAlign::__adjustIndels( int nshift )
{
    // adjust indels
    for (AlignPosList::iterator it = summary.dlist.begin(); it != summary.dlist.end(); ++it ) {
        it->aln_pos -= nshift; it->ref_pos -= nshift; it->seq_pos -= nshift;
    }
    for (AlignPosList::iterator it = summary.ilist.begin(); it != summary.ilist.end(); ++it ) {
        it->aln_pos -= nshift; it->ref_pos -= nshift; it->seq_pos -= nshift;
    }
}

//-DEFGH
//CD--GH
//>>>>>>>
//DEFGH
//-CDGH
//=========
//AB---FGHI
//--CDEFGHI
//>>>>>>>>>
//  -ABFGHI
//  CDEFGHI
//=================================
//     .    :   
// A---EF----KLMN
// -BCD--GHIJKLMN
//>>>>>>>>>>>>>>>
//     .    :
// ----AEFKLMN
// BCDGHIJKLMN
int SemiGlobalAlign::refineHead()
{
    if ( summary.lgap.first == 0 && summary.lgap.second == 0 ) return 0;

    int gap1 = 0;
    int gap2 = 0;
    
    for ( int i = 0; i < summary.range.first; i++ ) {
        if ( countGaps(row(aln,0),i) > 0 ) gap1++;
        if ( countGaps(row(aln,1),i) > 0 ) gap2++;
    }

    if ( gap1 == 0 || gap2 == 0 ) return 0;

    /* no of new leading gaps */
    int ngaps = gap1 - gap2;

    /* no of bases to be reduced */
    int shrink = gap1 > gap2 ? gap2 : gap1;

    // adjust indels
    __trimIndels( 0, summary.range.first );
    __adjustIndels(shrink);

    // adjust leading gaps
    summary.lgap = ngaps > 0 ? IntPair(ngaps, 0) : IntPair(0, -1*ngaps);

    //summary.range = IntPair( summary.range.first-shrink, summary.range.second-shrink );
    summary.s1se  = ngaps > 0 ? IntPair(ngaps, summary.s1se.second-shrink) : IntPair(0, summary.s1se.second-shrink);
    summary.s2se  = ngaps > 0 ? IntPair(0, summary.s2se.second-shrink) : IntPair(-1*ngaps, summary.s2se.second-shrink);

    // Alignment range
    summary.range.first -= (2*shrink);
    summary.range.second -= shrink;
    summary.outer.second -= shrink;

    //std::cout << "Head refinement:" << shrink << "\n";
    return shrink;
}

// shift: previous shift from refining head
void SemiGlobalAlign::refineTail(int shift)
{
    if ( summary.egap.first == 0 && summary.egap.second == 0 ) return;

    int gap1 = 0;
    int gap2 = 0;
    
    int start = summary.range.second+1+shift;
    int end = summary.outer.second+shift;
    for ( int i = start; i <= end; i++ ) {
        if ( countGaps(row(aln,0),i) > 0 ) gap1++;
        if ( countGaps(row(aln,1),i) > 0 ) gap2++;
    }

    if ( gap1 == 0 || gap2 == 0 ) return;

    /* no of new leading gaps */
    int ngaps = gap1 - gap2;

    /* no of bases to be reduced */
    int shrink = gap1 > gap2 ? gap2 : gap1;

    // adjust indels
    __trimIndels( start-shift, end-shift );
    

    
    //std::cout << "old1:" << summary.s1se.second << "\told2:" << summary.s2se.second << "\n";
    summary.s1se.second -= gap1 > gap2 ? abs(gap1-summary.egap.first) : abs(ngaps-summary.egap.first);
summary.s2se.second -= gap1 > gap2 ? abs(gap2-summary.egap.second) : abs(ngaps-summary.egap.second);    
//std::cout << "new1:" << summary.s1se.second << "\tnew2:" << summary.s2se.second << "\n";

summary.range.second += shrink; // range increases
summary.outer.second -= shrink; // outer decreases

    // adjust ending gaps
    summary.egap = ngaps > 0 ? IntPair(ngaps, 0) : IntPair(0, -1*ngaps);

    //std::cout << "Tail refinement:" << shrink << "\n";

}

// bool SemiGlobalAlign::adjustHead()
// {
//     if ( summary.lgap.first == 0 && summary.lgap.second == 0 ) return true;
    
//     int lg1 = 0;
//     int lg2 = 0;
//     for ( int i = 0; i < summary.range.first; i++ ) {
//         if ( countGaps(row(aln,0),i) > 0 ) lg1++;
//         if ( countGaps(row(aln,1),i) > 0 ) lg2++;
//     }

//     if ( lg1 == 0 || lg2 == 0 ) return true;

//     std::cout << "lg1:" << lg1 << "\tlg2:" << lg2 << "\n";
//     std::cout << "Adjusting leading gaps\n";
//     std::cout << "Before:\n";
//     summary.print(std::cout);

//     int min = lg1; 
//     if ( lg2 < lg1 ) min = lg2;
    
//     int diff = abs(lg1-lg2);

//     std::cout << "min:" << min << "\tdiff:" << diff << "\n";

//     //AB---FGHI
//     //--CDEFGHI
//     //=========
//     //  -ABFGHI
//     //  CDEFGHI
//     //=========
//     //---DEFG
//     //ABC--FG
//     //=======
//     //-DEFG
//     //ABCFG
//     if ( summary.s1se.first == 0 && summary.s2se.first == 0 ) {
//         std::cout << "Warning: unable to correct alignment end\n";
//         return false;
//     }
    
//     if ( lg1 >= lg2 ) {
//         summary.s1se.first = diff;
//         summary.s2se.first = 0; 
        
//         summary.lgap.first  = diff;
//         summary.lgap.second = 0;
//     }
//     else {
//         summary.s2se.first = diff;
//         summary.s1se.first = 0;
        
//         summary.lgap.first  = 0;
//         summary.lgap.second = diff;
//     }
    
//     summary.s1se.second -= min;
//     summary.s2se.second -= min;
    
//     summary.range.first  = diff;
//     summary.range.second -= min;
    
//     summary.outer.second -= min;
    
//     // length??
    
//     for (AlignPosList::iterator it = summary.dlist.begin(); it != summary.dlist.end();  ) {
//         it->aln_pos -= min; it->ref_pos -= min; it->seq_pos -= min;
//         if ( it->aln_pos <= summary.range.first ) it = summary.dlist.erase(it++);
//         else ++it;
//     }
    
//     for (AlignPosList::iterator it = summary.ilist.begin(); it != summary.ilist.end();  ) {
//         it->aln_pos -= min; it->ref_pos -= min; it->seq_pos -= min;
//         if ( it->aln_pos <= summary.range.first ) it = summary.ilist.erase(it++);
//         else ++it;
//     }

//     std::cout << "After:\n";
//     summary.print(std::cout);

//     return true;
// }
