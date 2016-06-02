#!/usr/bin/perl -w
use strict;

my @contents;
foreach(</usr/local/depot/projects/SPA/czhong/Grasp/Results/SGO_*.bst>){
    push @contents, $_;
}

my $i;
my $j;
for($i = 0; $i < scalar(@contents) - 1; ++ $i){
    for($j = $i + 1; $j < scalar(@contents); ++ $j){
	system "perl count_overlap_Blast.pl $contents[$i] $contents[$j]";
    }
}
