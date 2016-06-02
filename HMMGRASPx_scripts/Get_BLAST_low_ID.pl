#!/usr/bin/perl -w
use strict;

# parsing a BLAST result and print queries that have no HSPs (0-hits or identity < 75%)

my $file = shift;
open my $IN, "<$file" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  if(/\# Query\:\s(\S+)/)  {
    my $name = $1;
    my $tmp = <$IN>;
    while(!($tmp =~ /\# (\d+) hits found/)) {
      $tmp = <$IN>;
    }
    $tmp =~ /\# (\d+) hits found/;
    my $num_hit = $1;
    if($num_hit > 0)  {
      $tmp = <$IN>;
      my @decom = split /\s+/, $tmp;
      if($decom[2] < 80)  {
        print "$name\n";
      }
    }  else  {
      print "$name\n";
    }
  }
}
close $IN;
