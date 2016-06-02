#!/usr/bin/perl -w
use strict;

my $file = shift;

my @contents;
open my $IN, "<$file" or die "Cannot open file: $!\n";

while(<$IN>)	{
	chomp;
	push @contents, $_;
}

close $IN;

my $i;
my $sequence = "";
my $header = "";
for($i = 0; $i < scalar(@contents); ++ $i)	{
	if($contents[$i] =~ /^>/)	{
		if($header ne "" && $sequence ne "")	{
			print "$header\n$sequence\n";
		}
		$header = "";
		$sequence = "";
		my @decom = split /\s+/, $contents[$i];
		$header = $decom[0];
	}	else	{
		$sequence .= $contents[$i];
	}
}

print "$header\n$sequence\n";
