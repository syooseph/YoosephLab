#!/usr/bin/perl -w
use strict;

my $file = shift;

open my $IN, "<$file" or die "Cannot open file: $!\n";
while(<$IN>)	{
	chomp;
	print "$_\n";
	system "cat /home/cczhong/Works/GuidedAssemble/Data/Cellulase_Close/Grasp/2OSW_A.fa.N14.grasp.summary | grep '$_' |  wc -l";
	system "cat /home/cczhong/Works/GuidedAssemble/Data/Cellulase_Close/Short/2OSW_A.fa.N14.blast.summary | grep '$_:' |  wc -l";
}
close $IN;
