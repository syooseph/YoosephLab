#!/usr/bin/perl -w
use strict;

my $query_file = shift;	# the query sequence
my $target_file = shift;	# the target multi-fasta sequence
my $outfile = shift;	# the file for blast output

if(-e $outfile)	{
	die "#The output file exists.";
}

open my $IN, "<$target_file" or die "Cannot open file: $!\n";
while(<$IN>)	{
	chomp;
	if(/^>/)	{
		system "echo \"$_\" >.temp";
		my $seq = <$IN>;
		chomp $seq;
		system "echo \"$seq\" >>.temp";
	}
	system "/usr/local/packages/ncbi-blast+-2.2.28/bin/blastp -query $query_file -subject .temp -outfmt 7 >>$outfile";
}
close $IN;
