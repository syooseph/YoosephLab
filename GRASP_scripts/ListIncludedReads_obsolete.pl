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
  if(exists $aligned_location{$decom[1]})  {
    push @{$aligned_location{$decom[1]}}, $decom[8] - 1;
    push @{$aligned_location{$decom[1]}}, $decom[9] - 1;
  } else  {
    $aligned_location{$decom[1]} = [$decom[8] - 1, $decom[9] - 1];
  }
}
close $IN;

#foreach(keys %aligned_location) {
#  print "$_\n";
#  print "@{$aligned_location{$_}}\n";
#}
#die;

my $num_reads = 0;
open my $RIN, "<$read_file" or die "Cannot open read file: $!\n";
while(<$RIN>) {
  chomp;
  my $full_header = $_;
  if(/^>(.*)\:\:(\d+)-(\d+)/)  {
    my $id = $1;
    my $begin = $2;
    my $end = $3;
    if(exists $aligned_location{$id}) {
      my $i;
      for($i = 0; $i < scalar(@{$aligned_location{$id}}); $i += 2) {
        my $recruit_begin = ${aligned_location{$id}}[$i];
        my $recruit_end = ${aligned_location{$id}}[$i + 1];
        if(!($recruit_begin > $end || $recruit_end < $begin))  {
          $full_header =~ /^>(.*)/;
          print "$1\n";
        }
      }   
    }
  }
}

#
#
#(($aligned_location{$id}[0] <= $begin && $aligned_location{$id}[1] >= $begin) || 
#         ($aligned_location{$id}[0] <= $end && $aligned_location{$id}[1] >= $end)))  {
#      $full_header =~ /^>(.*)/;
#      print "$1\n";

close $RIN;
#print "$num_reads\n";
