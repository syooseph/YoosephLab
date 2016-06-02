#!/usr/bin/perl -w
use strict;

my $file = shift;	# this is expected to be a kegg gene-pfam file

my $results_dir = "/usr/local/depot/projects/SPA/czhong/Grasp/Results";
my $hmmer_exe = "/usr/local/packages/hmmer3/bin/hmmsearch";
my $hmm_dir = "/usr/local/depot/projects/SPA/czhong/Data/Pfam_HMM";

open my $IN, "<$file" or die "Cannot open file: $!\n";
while(<$IN>)	{
	chomp;
	/\:(.*)\s+.*\:(.*)/;
	my $gene = $1;
	my $pfam = $2;
	if(-e "$results_dir/$gene.asm.faa")	{
		print "running hmmsearch on $_ ...\n";
		system "$hmmer_exe $hmm_dir/$pfam.hmm $results_dir/$gene.asm.faa >$results_dir/$gene-$pfam.hmm.search";
	}
}
close $IN;
