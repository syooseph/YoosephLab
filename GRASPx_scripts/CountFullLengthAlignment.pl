#!/usr/bin/perl -w
use strict;

my $blast_file = shift;
my $full_len = shift;

my %aligned_location;
open my $IN, "<$blast_file" or die "Cannot open blast result file: $!\n";
my $num_reads = 0;
while(<$IN>) {
  chomp;
  if(/^\#/)  {
    next;
  }
  my @decom = split /\s+/, $_;
  my $span = $decom[9] - $decom[8] + 1;
  if($span == $full_len)  {
    ++ $num_reads;
  }  
}
close $IN;
print "$num_reads\n";
