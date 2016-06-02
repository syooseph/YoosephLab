#!/usr/bin/perl -w
use strict;

my $mfasta = shift;
my $ID_to_filter = shift;

open my $IN, "<$mfasta" or die "Cannot open file: $!\n";

while(<$IN>)	{
	chomp;
	if(/^>/ && !(/$ID_to_filter/))	{
		print "$_\n";
		my $sequence = <$IN>;
		print $sequence;
	}
}
close $IN;
