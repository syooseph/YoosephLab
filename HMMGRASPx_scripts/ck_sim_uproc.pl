#!/usr/bin/perl -w
use strict;

my $raw_file = shift;		# the file contains the raw reads (input) in FASTA format
my $uproc_map = shift;		# the uproc output in HMM-GRASPx map format

my $total = 0;
open my $IN, "<$raw_file" or die "Cannot open file: $!\n";
while(<$IN>)  {
  ++ $total;
}
close $IN;

$total = $total / 2;

open my $MIN, "<$uproc_map" or die "Cannot open file : $!\n";
my $ntp = 0;
my $nfp = 0;
while(<$MIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  my @decom1 = split /\:\:/, $decom[1];
  my @decom2 = split /\|\|/, $decom[2];
  if($decom1[2] eq $decom2[1])  {
    ++ $ntp;
  }  else  { 
    ++ $nfp;
  }
}
close $MIN;

my $recall = $ntp / $total;
my $precision = $ntp / ($ntp + $nfp);

print "#TP	FP	TOTAL	RECALL	PRECISION\n";
print "$ntp	$nfp	$total	$recall	$precision\n";
