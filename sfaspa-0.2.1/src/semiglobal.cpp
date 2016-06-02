#include "semiglobal.h"


SemiGlobalAlign::SemiGlobalAlign() 
{
    verbose = false;
    banded = false;
}

SemiGlobalAlign::SemiGlobalAlign(TString &sbjct, TString &query, int anchor)
{
    verbose = false;
    banded = false;
    align(sbjct, query, anchor);
    setSummary();
    refine();
}

SemiGlobalAlign::SemiGlobalAlign(std::string &sbjct, std::string &query, int anchor)
{
    verbose = false;    
    banded = false;
    align(sbjct, query, anchor);
    setSummary();
    refine();
}

SemiGlobalAlign::SemiGlobalAlign(std::string &sbjct, std::string &query, int anchor, int gex, int gop)
{
    verbose = false;
    banded = false;
    setGapPenalty(gex, gop);
    align(sbjct, query, anchor);
    setSummary();
    refine();
}

SemiGlobalAlign::SemiGlobalAlign(std::string &sbjct, std::string &query, int anchor, int gex, int gop, int ldig, int udig)
{
    verbose = false;
    banded = true;
    setGapPenalty(gex, gop);
    lower = ldig;
    upper = udig;
    align(sbjct, query, anchor);
    setSummary();
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

    // semi-global setup
    // 1. true: no penalty for leading gaps in second sequence 
    // 2. false: penalty for leading gaps in first sequence
    // 3. false: penalty for trailing gaps in first sequence
    // 3. true: no penalty for trailing gaps in second sequence

    Blosum62 score_type(gapext, gapopen);

    resize(rows(aln), 2);
    assignSource(row(aln,0),sbjct);
    assignSource(row(aln,1),query);

    /* old seqan v1.3 */
    // appendValue(rows(aln), sbjct);
    // appendValue(rows(aln), query);
    if ( anchor == ANCHOR_CENTER ) {
        // semi-global setup
        // Align query to center of sbjct
        //======================================
        // -------------xxooooooxx---            sbjct
        //              xxooooooxx               query 
        //======================================
        // 1. true: no penalty for leading gaps in second sequence 
        // 2. false: penalty for leading gaps in first sequence
        // 3. false: penalty for trailing gaps in first sequence
        // 3. true: no penalty for trailing gaps in second sequence
        AlignConfig<true,false,false,true> ac;    
        summary.score = ! banded ? 
            globalAlignment(aln, score_type, ac, Gotoh()):
            globalAlignment(aln, score_type, ac, lower, upper, Gotoh());
    } 
    else if ( anchor == ANCHOR_RIGHT ) {
        // Align query to suffix of sbjct
        //======================================
        // ---------------oooooo                 sbjct
        //                oooooo---------------- query
        //======================================
        AlignConfig<false,true,false,true> ac;    
        summary.score = ! banded ?
            globalAlignment(aln, score_type, ac, Gotoh()):
            globalAlignment(aln, score_type, ac, lower, upper, Gotoh());
        
    } 
    else if ( anchor == ANCHOR_LEFT ) {
        // Align query to prefix of sbjct
        //======================================
        //                oooooo---------------- sbjct
        // ---------------oooooo                 query
        //======================================
        AlignConfig<true,false,true,false> ac;    
        summary.score = ! banded ?
            globalAlignment(aln, score_type, ac, Gotoh()) :
            globalAlignment(aln, score_type, ac, lower, upper, Gotoh());
    }
}
	
void SemiGlobalAlign::align(std::string &sbjct, std::string &query, int anchor) 
{
    TString tstr1 = TString(sbjct);
    TString tstr2 = TString(query);

    align(tstr1, tstr2, anchor);
}


void SemiGlobalAlign::setSummary()
{
    using namespace seqan;

    int s = clippedBeginPosition(row(aln,0));
    int e = clippedEndPosition(row(aln,0));
    summary.outer = intPair(s,e-1);
    if ( verbose ) std::cout << "Outer:" << s << " " << e-1 << "\n";
    if ( verbose ) {
        std::cout << "Clipped positions:" << s << " " << e << "\n";
        s = clippedBeginPosition(row(aln,1));
        e = clippedEndPosition(row(aln,1));
        std::cout << "Clipped positions:" << s << " " << e << "\n";

        std::cout << "Sbjct begin/end positions:" << beginPosition(row(aln, 0)) << " " << endPosition(row(aln, 0)) << "\n";
        std::cout << "Query begin/end positions:" << beginPosition(row(aln, 1)) << " " << endPosition(row(aln, 1)) << "\n";
        std::cout << "\n";
    }
    setSequenceRange();
    setLeadingGaps();
    setEndGaps();

    alignmentInnerRange();


    setGapPositions();
    setAlignmentStats();
}

void SemiGlobalAlign::alignmentInnerRange()
{
    using namespace seqan;
    
    Gaps<TString> g1(row(aln,0)), g2(row(aln,1));
    Iterator< Gaps<TString> >::Type its = begin(g1), ite = end(g1);
    Iterator< Gaps<TString> >::Type jts = begin(g2), jte = end(g2);

    int s = summary.outer.first;
    while( countGaps(its) > 0 || countGaps(jts) > 0 ) {
        if ( its == ite || jts == jte ) break;
        its++; jts++; s++;
    }
    
    int e = summary.outer.second;
    ite--; jte--;
    while( countGaps(ite) > 0 || countGaps(jte) > 0 ) {
        if ( atBegin(ite) || atBegin(jte) ) break;
        ite--; jte--; e--;
    }
    
    summary.range = intPair(s,e);

    if ( verbose ) std::cout << "range:" << s << " " << e << "\n";
}

void SemiGlobalAlign::setSequenceRange()
{
    using namespace seqan;

    summary.s1se = intPair( beginPosition(row(aln, 0)), endPosition(row(aln, 0)) );
    summary.s2se = intPair( beginPosition(row(aln, 1)), endPosition(row(aln, 1)) );

    // Gaps<TString> g1(row(aln,0)), g2(row(aln,1));
    // Iterator< Gaps<TString> >::Type its = begin(g1), ite = end(g1)-1;
    // Iterator< Gaps<TString> >::Type jts = begin(g2), jte = end(g2)-1;
    // int s1 = summary.outer.first, e1 = summary.outer.second;
    // int s2 = summary.outer.first, e2 = summary.outer.second;
    
    // while( countGaps(its) > 0 ) {
    //     if ( its == end(g1) ) break;
    //     s1++; its++;
    // }
    // while( countGaps(ite) > 0 ) {
    //     if ( ite == begin(g1) ) break;
    //     e1--; ite--;
    // }
    // while( countGaps(jts) > 0 ) {
    //     if ( jts == end(g2) ) break;
    //     s2++; jts++;
    // }
    // while( countGaps(jte) > 0 ) {
    //     if ( jte == begin(g2) ) break;
    //     e2--; jte--;
    // }

    // if ( verbose ) std::cout << "s1:" << s1 << ",e1:" << e1 << "\ts2:" << s2 << ",e2:" << e2 << "\n";

    // summary.s1se = intPair(s1,e1);
    // summary.s2se = intPair(s2,e2);
}

void SemiGlobalAlign::printGaps(std::ostream &out)
{
    for ( AlignPosList::iterator it = summary.ilist.begin(); it != summary.ilist.end(); ++it )
        out << "Insertion:\ta:" << it->aln_pos << "\tr:" << it->ref_pos << "\tq:" << it->qry_pos << "\n";
    for ( AlignPosList::iterator it = summary.dlist.begin(); it != summary.dlist.end(); ++it )
        out << "Deletion :\ta:" << it->aln_pos << "\tr:" << it->ref_pos << "\tq:" << it->qry_pos << "\n";
}

void SemiGlobalAlign::setLeadingGaps()
{
    using namespace seqan;

    Gaps<TString> g1(row(aln,0)), g2(row(aln,1));
    Iterator< Gaps<TString> >::Type its = begin(g1), ite = end(g1);
    Iterator< Gaps<TString> >::Type jts = begin(g2), jte = end(g2);

    int lgap = 0;
    while( its!=ite ) {
        if ( countGaps(its) > 0 && countGaps(jts) == 0 ) lgap++;
        else break;
        its++; jts++;
    }
    summary.lgap.first = lgap;

    lgap = 0;
    its = begin(g1), jts = begin(g2);
    while( jts!=jte ) {
        if ( countGaps(jts++) > 0 && countGaps(its++) == 0 ) lgap++;
        else break;
    }
    summary.lgap.second = lgap;
}

void SemiGlobalAlign::setEndGaps()
{
    using namespace seqan;
    Gaps<TString> g1(row(aln,0)), g2(row(aln,1));
    Iterator< Gaps<TString> >::Type its = begin(g1), ite = end(g1)-1;
    Iterator< Gaps<TString> >::Type jts = begin(g2), jte = end(g2)-1;

    int egap = 0;
    while( its!=ite ) {
        if ( countGaps(ite--) && !countGaps(jte--) ) egap++;
        else break;
    }
    summary.egap.first = egap;

    egap = 0;
    ite = end(g1)-1, jte = end(g2)-1;
    while( jts!=jte ) {
        if ( countGaps(jte--) && !countGaps(ite--) ) egap++;
        else break;
    }
    summary.egap.second = egap;
}

void SemiGlobalAlign::setAlignmentStats()
{
    using namespace seqan;
    Gaps<TString> g1(row(aln,0)), g2(row(aln,1));
    Iterator< Gaps<TString> >::Type its = begin(g1), ite = end(g1);
    Iterator< Gaps<TString> >::Type jts = begin(g2), jte = end(g2);
    its += summary.range.first;
    jts += summary.range.first;

    for ( int i = summary.range.first; i <= summary.range.second; i++ ) {
        if ( its == ite || jts == jte ) break;

        if ( countGaps(its) == 0 && countGaps(jts) == 0 ) {
            TString s1 = getValue(its);
            TString s2 = getValue(jts);
            char c1 = s1[0];
            char c2 = s2[0];
            if ( c1 == c2 ) { summary.match++; summary.positive++; }
            else { 
                summary.mismatch++;
                if ( scoring::getScore(c1, c2, BLOSUM62) > 0 ) summary.positive++;
            }
        }
        its++; jts++;
    }
    summary.length = summary.range.second - summary.range.first + 1;//e-s+1;
    summary.posrate = (double)summary.positive/summary.length;
}

void SemiGlobalAlign::setGapPositions()
{
    using namespace seqan;
    Gaps<TString> g1(row(aln,0)), g2(row(aln,1));
    Iterator< Gaps<TString> >::Type its = begin(g1), ite = end(g1);
    Iterator< Gaps<TString> >::Type jts = begin(g2), jte = end(g2);

    int spos = 0; // sbjct_pos
    int qpos = 0; // query_pos

    for ( int i = summary.outer.first; i<= summary.outer.second; i++ ) {
        if ( its == ite || jts == jte ) break;
        int sgap = countGaps(its);
        int qgap = countGaps(jts);
        
        if ( sgap > 0  ) summary.ilist.push_back( AlignIndex(i, qpos, spos) );
        if ( qgap > 0  ) summary.dlist.push_back( AlignIndex(i, qpos, spos) );

        if ( sgap == 0 ) spos++;
        if ( qgap == 0 ) qpos++;

        its++; jts++;
    }
    if ( verbose ) {
        for ( AlignPosList::iterator it = summary.ilist.begin(); it != summary.ilist.end(); ++it )
            std::cout << "Insertion:\ta:" << it->aln_pos << "\tr:" << it->ref_pos << "\tq:" << it->qry_pos << "\n";
        for ( AlignPosList::iterator it = summary.dlist.begin(); it != summary.dlist.end(); ++it )
            std::cout << "Deletion :\ta:" << it->aln_pos << "\tr:" << it->ref_pos << "\tq:" << it->qry_pos << "\n";
    }
        
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

void SemiGlobalAlign::refine()
{
    int hshift = refineHead();
    int tshift = refineTail();
    trimIndels(hshift, tshift);
}

void SemiGlobalAlign::trimIndels(int hoff, int toff)
{
    for (AlignPosList::iterator it = summary.dlist.begin(); it != summary.dlist.end();  ) {
        if ( it->aln_pos <= summary.range.first + hoff ) summary.dlist.erase(it++);
        else break;
    }
    for (AlignPosList::iterator it = summary.ilist.begin(); it != summary.ilist.end();  ) {
        if ( it->aln_pos <= summary.range.first + hoff ) summary.ilist.erase(it++);
        else break;
    }

    for (AlignPosList::iterator it = summary.dlist.begin(); it != summary.dlist.end();  ) {
        if ( it->aln_pos >= summary.range.second - toff ) summary.dlist.erase(it++);
        else ++it;
    }
    for (AlignPosList::iterator it = summary.ilist.begin(); it != summary.ilist.end();  ) {
        if ( it->aln_pos >= summary.range.second - toff ) summary.ilist.erase(it++);
        else ++it;
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
    using namespace seqan;

    if ( summary.lgap.first == 0 && summary.lgap.second == 0 ) return 0;

    Gaps<TString> g1(row(aln,0)), g2(row(aln,1));
    Iterator< Gaps<TString> >::Type its = begin(g1), ite = end(g1);
    Iterator< Gaps<TString> >::Type jts = begin(g2), jte = end(g2);

    int gap1 = 0;
    int gap2 = 0;
    
    for ( int i = 0; i < summary.range.first; i++ ) {
        if ( its == ite || jts == jte ) break;
        if ( countGaps(its++) > 0 ) gap1++;
        if ( countGaps(jts++) > 0 ) gap2++;
    }

    if ( verbose ) std::cout << "gap1:" << gap1 << "\tgap2:" << gap2 << "\n";

    if ( gap1 == 0 || gap2 == 0 ) return 0;

    /* no of new leading gaps */
    int ngaps = gap1 - gap2;

    /* no of bases to be reduced */
    int shrink = gap1 > gap2 ? gap2 : gap1;


    // adjust leading gaps
    summary.lgap = ngaps > 0 ? intPair(ngaps,0) : intPair(0,-1*ngaps);
    //summary.lgap = ngaps > 0 ? intPair(0, ngaps) : intPair(-1*ngaps, 0) ;

    // if ( verbose ) {
    //     std::cout << "Head refinement:" << shrink << "\n";
    //     std::cout << "old1:" << summary.s1se.first << "\t" << summary.s1se.second << "\n";
    //     std::cout << "old2:" << summary.s2se.first << "\t" << summary.s2se.second << "\n";
    // }
    // summary.s1se  = ngaps > 0 ? intPair(ngaps, summary.s1se.second-shrink) : intPair(0, summary.s1se.second-shrink);
    // summary.s2se  = ngaps > 0 ? intPair(0, summary.s2se.second-shrink) : intPair(-1*ngaps, summary.s2se.second-shrink);
    // if ( verbose ) {
    //     std::cout << "new1:" << summary.s1se.first << "\t" << summary.s1se.second << "\n";
    //     std::cout << "new2:" << summary.s2se.first << "\t" << summary.s2se.second << "\n";
    // }


    // Alignment range
    if ( verbose ) std::cout << "Old Range:" << summary.range.first << "\t" << summary.range.second << "\n";
    summary.range.first -= (2*shrink);
    summary.range.second -= shrink;
    summary.outer.second -= shrink;

    if ( verbose ) {
        std::cout << "Range:" << summary.range.first << "\t" << summary.range.second << "\n";
        std::cout << "Outer:" << summary.outer.first << "\t" << summary.outer.second << "\n";
    }

    return shrink;
}

int SemiGlobalAlign::refineTail()
{
    using namespace seqan;

    if ( summary.egap.first == 0 && summary.egap.second == 0 ) return 0;

    Gaps<TString> g1(row(aln,0)), g2(row(aln,1));
    Iterator< Gaps<TString> >::Type its = begin(g1), ite = end(g1)-1;
    Iterator< Gaps<TString> >::Type jts = begin(g2), jte = end(g2)-1;

    int gap1 = 0;
    int gap2 = 0;

    while(true) {
        if ( ite == its || jte == jts ) break;
        int g1 = countGaps(ite--);
        int g2 = countGaps(jte--);
        if ( g1 > 0 ) gap1++;
        if ( g2 > 0 ) gap2++;
        if ( g1+g2==0) break;
    }
    
    
    if ( gap1 == 0 || gap2 == 0 ) return 0;

    /* no of new leading gaps */
    int ngaps = gap1 - gap2;

    /* no of bases to be reduced */
    int shrink = gap1 > gap2 ? gap2 : gap1;


    // if ( verbose ) {
    //     std::cout << "old1:" << summary.s1se.first << "\t" << summary.s1se.second << "\n";
    //     std::cout << "old2:" << summary.s2se.first << "\t" << summary.s2se.second << "\n";
    // }
    // summary.s1se.second -= gap1 > gap2 ? abs(gap1-summary.egap.first) : abs(ngaps-summary.egap.first);
    // summary.s2se.second -= gap1 > gap2 ? abs(gap2-summary.egap.second) : abs(ngaps-summary.egap.second);    
    // if ( verbose ) {
    //     std::cout << "new1:" << summary.s1se.first << "\t" << summary.s1se.second << "\n";
    //     std::cout << "new2:" << summary.s2se.first << "\t" << summary.s2se.second << "\n";
    // }
    
    summary.range.second += shrink; // range increases
    summary.outer.second -= shrink; // outer decreases
    
    // adjust ending gaps
    //summary.egap = ngaps > 0 ? intPair(0, ngaps) : intPair(-1*ngaps, 0) ;
    summary.egap = ngaps > 0 ? intPair(ngaps,0) : intPair(0,-1*ngaps);
    
    if ( verbose ) 
        std::cout << "Tail refinement:" << shrink << "\n";

    return shrink;
}

