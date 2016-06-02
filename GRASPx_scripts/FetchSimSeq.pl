#!/usr/bin/perl -w
use strict;

my $file = "/usr/local/depot/projects/SPA/czhong/Works/GRASPx/Data/marine_sim_glycolysis/sim.fa";
my $range1 = 1884000;
my $range2 = 1888000;

open my $IN, "<$file" or die "Cannot open file: $!\n";

while(<$IN>)  {
  chomp;
  if(/^>/)  {
    my @decom = split /\_/, $_;
    #print "$decom[1]	$decom[2]\n";
    if($decom[2] >= $range1 && $decom[2] <= $range2 && $decom[3] >= $range1 && $decom[3] <= $range2)  {
      print "$_\n";
      my $temp = <$IN>;
      print $temp;
    }
  }
}
close $IN;
