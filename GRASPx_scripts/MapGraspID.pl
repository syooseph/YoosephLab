#!/usr/bin/perl -w
use strict;

my $fasta_file = shift;
my $grasp_file = shift;
my $cutoff = shift;

my @ID;
open my $IN, "<$fasta_file" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  if(/^>(.*)/)  {
    push @ID, $1;
  }
}
close $IN;

open my $GIN, "<$grasp_file" or die "Cannot open file: $!\n";
my $temp = <$GIN>;
while(<$GIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  if($decom[2] <= $cutoff)  {
    print "$ID[$decom[1]]\t$decom[2]\n";
  }
}
