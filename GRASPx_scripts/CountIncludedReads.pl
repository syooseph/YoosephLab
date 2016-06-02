#!/usr/bin/perl -w
use strict;

my $blast_file = shift;
my $read_file = shift;

my %aligned_location;
open my $IN, "<$blast_file" or die "Cannot open blast result file: $!\n";
while(<$IN>) {
  chomp;
  if(/^\#/)  {
    next;
  }
  my @decom = split /\s+/, $_;
  $aligned_location{$decom[1]} = [$decom[8] - 1, $decom[9] - 1];
}
close $IN;

my $num_reads = 0;
open my $RIN, "<$read_file" or die "Cannot open read file: $!\n";
while(<$RIN>) {
  chomp;
  if(/^>(.*)\s+(\d+)\:(\d+)/)  {
    my $id = $1;
    my $begin = $2;
    my $end = $3;
    if(exists $aligned_location{$id} && 
        $aligned_location{$id}[0] <= $begin && 
        $aligned_location{$id}[1] >= $end)  {
      ++ $num_reads;
    }
  }
}
close $RIN;
print "$num_reads\n";
