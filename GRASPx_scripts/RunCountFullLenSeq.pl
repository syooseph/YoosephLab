#!/usr/bin/perl -w
use strict;

my $folder = shift;
my $len = shift;

foreach(<$folder/*.summary>) {
  #print "$_ ";
  system "perl CountFullLengthAlignment.pl $_ $len";
  #system "perl CountIncludedReads.pl $_ /home/cczhong/Works/GuidedAssemble/Data/sim.PF02834.random.fa"
}
