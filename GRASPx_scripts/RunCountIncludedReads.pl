#!/usr/bin/perl -w
use strict;

my $folder = shift;
my $simulated_data = shift;

foreach(<$folder/*.summary>) {
  print "$_ ";
  #system "perl CountFullLengthAlignment.pl $_";
  system "perl CountIncludedReads.pl $_ $simulated_data"
}
