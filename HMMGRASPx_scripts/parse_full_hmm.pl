#!/usr/bin/perl -w
use strict;

my $hmm_full_file = shift;

open my $IN, "<$hmm_full_file" or die "Cannot open file: $!\n";

my $name;
my $id;
my $target;
while(<$IN>)  {
  chomp;
  if(/^Query:\s+(\S+)/)  {
    $name = $1;
  } elsif(/^Accession:\s+(PF\d+)\./)  {
    $id = $1;
  } elsif(/^>>\s+(\S+)/)  {
    $target = $1;
  } elsif(/\!/)  {
    my @decom = split /\s+/, $_;
    print "$name	$id	$target	$decom[10]	$decom[11]\n";
  }
}
