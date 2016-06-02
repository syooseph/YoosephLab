#!/usr/bin/perl -w
use strict;

my $included_file = shift;
my $read_blast_file = shift;

my %inc_hash;
open my $IIN, "<$included_file" or die "Cannot open included reads file: $!\n";
while(<$IIN>) {
  chomp;
  $inc_hash{$_} = 1;
}
close $IIN;

open my $RIN, "<$read_blast_file" or die "Cannot open blast file : $!\n";
while(<$RIN>) {
  chomp;
  my $line = $_;
  if(/^#/)  {
    next;
  }
  my @decom = split /\s+/, $line;
  if(!(exists $inc_hash{$decom[1]}))  {
    print "$line\n";
  }
}
close $RIN;
