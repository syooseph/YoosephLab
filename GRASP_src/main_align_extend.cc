#include "scoring_function.h"
#include "seq_align.h"
#include "seq_align_extend.h"
#include <iostream>
#include <string>

using namespace std;

int main()  {
  string seqA = "ACGACTATCATCC";
	string seqB = "ACAGTACGATT";
	string seed_seqA = "";
	string seed_seqB = "";
	string extend_seqA = "CGACGTAGCTAC";
	string extend_seqB = "CGACTGACT";
	//ScoringFunction<int>* simple_score = new ScoringFunction<int>;
	//cout << simple_score->mismatch << endl;
	
	
	/*
	ScoringFunction<int> simple_score;
	SeqAlign<int> simple_align(seqA, seqB, &simple_score, GLOBAL, 3, 5, false);
	simple_align.Align();
	simple_align.PrintAlignment();
	
	
	
	int i, j;
	for(i = 0; i <= seqA.length(); ++ i)	{
		for(j = 0; j <= seqB.length(); ++ j)	{
			cout << simple_align.AccessDPScoreMatrix(i, j, FULL) << "\t";
		}
		cout << endl;
	}
	
	
	int k = simple_align.GetAlignmentScore();
	cout << k << endl;
	*/
	
	
	vector<int> init_score_row(1); 
	init_score_row = {0};
	vector<int> init_score_column(1);
	init_score_column = {0};
	
	EdgeScoreBundle<int> in_score;
	in_score.row = init_score_row;
	in_score.row_del = init_score_row;
	in_score.row_ins = init_score_row;
	in_score.column = init_score_column;
	in_score.column_del = init_score_column;
	in_score.column_ins = init_score_column;
	in_score.begin_row = 0;
	in_score.begin_column = 0;
	
	ScoringFunction<int> simple_score_extend;
	SeqAlignExtend<int> simple_align_extend(
		seed_seqA, seed_seqB, extend_seqA, extend_seqB,
		in_score, &simple_score_extend, GLOBAL
	);
	simple_align_extend.Align();
	
	int i, j;
	for(i = 0; i < seed_seqA.length() + extend_seqA.length() + 1; ++ i)	{
		for(j = 0; j < seed_seqB.length() + extend_seqB.length() + 1; ++ j)	{
			cout << simple_align_extend.AccessDPScoreMatrix(i, j, FULL) << "\t";
		}
		cout << endl;
	}
  return 0;
}
