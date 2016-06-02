#!/usr/bin/perl -w
use strict;

my $fasta_file = shift;
my $id_file = shift;

my @ID;
open my $IN, "<$fasta_file" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  if(/^>(.*)/)  {
    push @ID, $1;
  }
}
close $IN;

open my $GIN, "<$id_file" or die "Cannot open file: $!\n";
my $temp = <$GIN>;
while(<$GIN>)  {
  chomp;
  print "$ID[$_]\n";
}
