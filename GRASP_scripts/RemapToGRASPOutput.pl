#!/usr/bin/perl -w
use strict;

my $remap_file = shift;

open my $IN, "<$remap_file" or die "Cannot open file: $!\n";

print "\# Fields: query_sequence_ID, read_ID, best_Evalue_aligned\n";
while(<$IN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  print "0	$decom[0]	$decom[3]\n";
}
