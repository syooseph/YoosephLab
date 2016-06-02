#!/usr/bin/perl -w
use strict;

my $folder = shift;

foreach(<$folder/*.summary>) {
  open my $IN, "<$_" or die "Cannot open file: $!";
  my $good = 0;
  while(<$IN>) {
    chomp;
    if(/^\#\s+(\d+)\s+hits\s+found/)  {
      print "$1\n";
      $good = 1;
    }
  }
  if(!$good)  {
    print "ERROR\n";
  }
  close $IN;
}
