#!/usr/bin/perl -w
use strict;

# filter out reads that map to contigs whose alignment E-value is above the cutoff

my $file = shift;
my $cutoff = shift;

open my $IN, "<$file" or die "Cannot open file: $!\n";
while(<$IN>)  {
  my $line = $_;
  chomp;
  my @decom1 = split /\s+/, $_;
  my @decom2 = split /\|\|/, $decom1[2];
  my @decom3 = split /\:\:/, $decom2[1];
  if($decom3[2] < $cutoff)  {
    print $line;
  }
}
close $IN;
