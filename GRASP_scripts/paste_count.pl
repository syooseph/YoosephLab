#!/usr/bin/perl -w
use strict;

my $file = shift;

my $index = -1;
my @specificities;
open my $IN, "<$file" or die "Cannot open file: $!\n";
my @contents;
while(<$IN>)	{
	chomp;
	push @contents, $_;
}
close $IN;

my $i;
for($i = 0; $i < scalar(@contents); ++ $i)	{
	if($contents[$i] =~ /^\#/)	{
		next;
	}
	if($contents[$i] =~ /^LBPG/)	{
		++ $index;
		push @{$specificities[$index]}, $contents[$i];
	}	else	{
		my @decom = split /\s+/, $contents[$i];
		my $spe = $decom[0];
		push @{$specificities[$index]}, $spe;
	}
}

my @overall;
push @overall, "overall";
my $j;
for($i = 0; $i < scalar(@specificities); ++ $i)	{
	#print "$specificities[$i][0]	";
	for($j = scalar(@{$specificities[$i]}) - 1; $j >= 1; -- $j)	{
		$overall[$j] += $specificities[$i][$j];
	#	print "$specificities[$i][$j]	";
	}
	#print "\n";
}

print "$overall[0]	";
for($i = 1; $i < scalar(@overall); ++ $i)	{
	my $p = $overall[$i];
	print "$p	";
}
print  "\n";
