#!/usr/bin/perl -w
use strict;
use warnings;

my $file = shift;

open my $IN, "<$file" or die "Cannot open file: $!\n";

my $temp = <$IN>;
$temp = <$IN>;

while(<$IN>)	{
	chomp;
	my @decom = split /\s+/, $_;
	my $i;
	print "	{";
	for($i = 1; $i < scalar(@decom); ++ $i)	{
		if($decom[$i] < 0)	{
			print "$decom[$i], ";
		}	else	{
			print " $decom[$i], ";
		}
	}
	print "\b\b},\n";
}
close $IN;
