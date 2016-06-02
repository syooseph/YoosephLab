#!/usr/bin/perl -w
use strict;

my $file = shift;
open my $IN, "<$file" or die "Cannot open file: $!\n";

while(<$IN>)  {
  chomp;
  if(/^(\>PF\d+)/)  {
    print "$1\n";
  }  else  {
    print "$_\n";
  }
}
close $IN;
