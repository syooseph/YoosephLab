#!/usr/bin/perl -w
use strict;

my $blast_file = shift;	 # the full-sequence search results for the specific query sequence
my $full_folder = shift; # the folder contains the blast full-sequence results
my $read_file = shift;   # the set of simulated read files

# reading information from the result that takes the sequence as query
my %aligned_location;
my $of_interest_ID;
open my $IN, "<$blast_file" or die "Cannot open blast result file: $!\n";
while(<$IN>) {
  chomp;
  if(/^\#/)  {
    next;
  }
  my @decom = split /\s+/, $_;
  if($decom[0] eq $decom[1] || $decom[10] > 1e-3)  {
    next;
  }
  $of_interest_ID = $decom[0];
  if(exists $aligned_location{$decom[1]})  {
    push @{$aligned_location{$decom[1]}}, $decom[8] - 1;
    push @{$aligned_location{$decom[1]}}, $decom[9] - 1;
    push @{$aligned_location{$decom[1]}}, $decom[11];
  } else  {
    $aligned_location{$decom[1]} = [$decom[8] - 1, $decom[9] - 1, $decom[11]];
  }
}
close $IN;

foreach(<$full_folder/*.blast.summary>) {
  open my $AIN, "<$_" or die "Cannot open blast additional file: $!\n";
  while(<$AIN>) {
    chomp;
    if(/^\#/)  {
      next;
    }
    my @decom = split /\s+/, $_;
    if($decom[0] eq $decom[1] || $decom[10] > 1e-3)  {
      next;
    }
    if($decom[1] eq $of_interest_ID)  {
      if(exists $aligned_location{$decom[0]})  {
        push @{$aligned_location{$decom[0]}}, $decom[6] - 1;
        push @{$aligned_location{$decom[0]}}, $decom[7] - 1;
        push @{$aligned_location{$decom[0]}}, $decom[11];
      } else  {
        $aligned_location{$decom[0]} = [$decom[6] - 1, $decom[7] - 1, $decom[11]];
      }
    }
  }
  close $AIN;
}

#foreach(keys %aligned_location) {
#  print "$_\n";
#  print "@{$aligned_location{$_}}\n";
#}
#die;

my $num_reads = 0;
open my $RIN, "<$read_file" or die "Cannot open read file: $!\n";
my %included_reads;
while(<$RIN>) {
  chomp;
  my $full_header = $_;
  if(/^>(.*)\:\:(\d+)-(\d+)/)  {
    my $id = $1;
    my $begin = $2;
    my $end = $3;
    if(exists $aligned_location{$id}) {
      my $i;
      for($i = 0; $i < scalar(@{$aligned_location{$id}}); $i += 3) {
        my $recruit_begin = ${aligned_location{$id}}[$i];
        my $recruit_end = ${aligned_location{$id}}[$i + 1];
        my $bit_score = ${aligned_location{$id}}[$i + 2];
        if(!($recruit_begin > $end || $recruit_end < $begin))  {
          $full_header =~ /^>(.*)/;
          #if(!(exists $included_reads{$1}) || $included_reads{$1} < $bit_score)  {
          $included_reads{$1} = $num_reads;
          #}
        }
      }   
    }
    ++ $num_reads;
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

foreach(keys %included_reads) {
  print "$_ $included_reads{$_}\n";
}
