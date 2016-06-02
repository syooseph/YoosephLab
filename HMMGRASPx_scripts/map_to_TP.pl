#!/usr/bin/perl -w
use strict;

my $map_file = shift;

open my $IN, "<$map_file" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  my @decom2 = split /\|\|/, $decom[2];
  print "$decom2[2]	$decom2[1]	>$decom[1]\n";
}
close $IN;
