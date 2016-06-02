#!/usr/bin/perl -w
use strict;

my $mfasta_file = shift;
open my $IN, "<$mfasta_file" or die "Cannot open file: $!\n";

my $min_id;
my $min_len = 999999;
my $max_id;
my $max_len = -1;

while(<$IN>)	{
	chomp;
	if(/^>/)	{
		my $header = $_;
		my $sequence = <$IN>;
		chomp $sequence;
		if(length($sequence) < $min_len)	{
			$min_len = length($sequence);
			$min_id = $header;
		}	elsif(length($sequence) > $max_len)	{
			$max_len = length($sequence);
			$max_id = $header;
		}
	}
}
close $IN;

print "$min_id	$min_len	$max_id	$max_len\n";
