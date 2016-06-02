#!/usr/bin/perl -w
use strict;

# expecting multi-fasta file
my $file = shift;

open my $IN, "<$file" or die "Cannot open file: $!\n";
my $n = 0;
while(<$IN>)  {
  if(/^>/)  {
    print ">$n\n";
    ++ $n;
  }  else  {
    print;
  }
}
