#!/usr/bin/perl -w
use strict;
#my $query_file = shift;	 # the query sequence in fasta format
my $blast_file = shift;	 # the full-sequence search results for the specific query sequence
my $read_file = shift;   # the set of simulated read files
my $map_file = shift;	# the mapping file tells where each reads is aligned to the full-sequences (read_remap.csv)

# reading sequence form the query file to resolve the length
#open my $QIN, "<$query_file" or die "Cannot open query file: $!\n";
#my $header = <$QIN>;
#my $sequence = <$QIN>;
#my $query_length = length($sequence) - 1;
#close $QIN;

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
  #print "$query_length	$decom[10]	$decom[9]\n";
  #if($decom[0] eq $decom[1] || abs($decom[10] - $decom[9]) < 0.8 * 3 * $query_length)  {
  if($decom[0] eq $decom[1])  {
    #print "region skipped\n";
    next;
  }
  $of_interest_ID = $decom[0];
  if(exists $aligned_location{$decom[1]})  {
    push @{$aligned_location{$decom[1]}}, $decom[9] - 1;
    push @{$aligned_location{$decom[1]}}, $decom[10] - 1;
    push @{$aligned_location{$decom[1]}}, $decom[11];
  } else  {
    $aligned_location{$decom[1]} = [$decom[9] - 1, $decom[10] - 1, $decom[11]];
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
my %included_reads;
my %map_reads;
my $index = 0;
while(<$RIN>) {
  chomp;
  if(/^>(.*)/)  {
    $map_reads{$1} = $index;
    ++ $index;
  }
}
close $RIN;

open my $MIN, "<$map_file" or die "Cannot open file: $!\n";
while(<$MIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  #print "$_\n";
  if(exists $aligned_location{$decom[1]}) {
    #print "exists!\n";
    my $i;
    my $begin = $decom[3];
    my $end = $decom[4];
    for($i = 0; $i < scalar(@{$aligned_location{$decom[1]}}); $i += 3) {
      my $recruit_begin = ${aligned_location{$decom[1]}}[$i];
      my $recruit_end = ${aligned_location{$decom[1]}}[$i + 1];
      my $evalue = ${aligned_location{$decom[1]}}[$i + 2];
      #print "$recruit_begin	$recruit_end	$begin	$end\n";
      #if($recruit_begin <= $recruit_end && !($recruit_begin - $end > 100 || $begin - $recruit_end > 100)  ||
      #   $recruit_begin >= $recruit_end && !($recruit_end - $begin > 100 || $end - $recruit_begin > 100)
      #)  {
      if(($recruit_begin <= $recruit_end && $recruit_begin <= $begin && $recruit_end >= $end)  ||
         ($recruit_begin >= $recruit_end && $recruit_end <= $end && $recruit_begin >= $begin)
      )  {
        if($evalue < 1e-10)  {
          #print "recruited\n";
          $included_reads{$decom[0]} = $map_reads{$decom[0]};
        }
      }
    }
  }
}
close $MIN;

=pod
foreach(keys %map_reads)  {
  my @decom = split /\_/;
  my $begin;
  my $end;
  my $full_ID = $_;
  $full_ID =~ /^(.*\|)\_(\d+)\_(\d+)\_.*\/(\d+)\_(\d+)\_(\d+)\_(.)$/;
  my $id = $1;
  my $frag_begin = $2;
  my $frag_end = $3;
  my $pair_loc = $4;
  my $pair_begin = $5;
  my $pair_end = $6;
  my $strand = $7;
  #print "$id	$frag_begin	$frag_end	$pair_loc	$pair_begin	$pair_end	$strand\n";
  if($pair_loc == 1)  {
    $begin = $frag_begin + $pair_begin - 1;
    $end = $frag_begin + $pair_end - 1;
  }  else  {
    $begin = $frag_end - $pair_end + 1;
    $end = $frag_end - $pair_begin + 1;
    if($strand eq '+')  {
      $strand = '-';
    }  else  {
      $strand = '+';
    }
  }

  if(exists $aligned_location{$id}) {
    my $i;
    for($i = 0; $i < scalar(@{$aligned_location{$id}}); $i += 3) {
      my $recruit_begin = ${aligned_location{$id}}[$i];
      my $recruit_end = ${aligned_location{$id}}[$i + 1];
      my $evalue = ${aligned_location{$id}}[$i + 2];
      if($strand eq '+' && $recruit_begin <= $recruit_end && !($recruit_begin > $end || $recruit_end < $begin)  ||
         $strand eq '-' && $recruit_begin >= $recruit_end && !($recruit_begin < $begin || $recruit_end > $end)
      )  {
        if($evalue < 1e-1)  {
          $included_reads{$full_ID} = $map_reads{$full_ID};
        }
      }
    }   
  }
  ++ $num_reads;
}
=cut
#
#
#(($aligned_location{$id}[0] <= $begin && $aligned_location{$id}[1] >= $begin) || 
#         ($aligned_location{$id}[0] <= $end && $aligned_location{$id}[1] >= $end)))  {
#      $full_header =~ /^>(.*)/;
#      print "$1\n";

#close $RIN;
#print "$num_reads\n";

foreach(keys %included_reads) {
  print "$_ $included_reads{$_}\n";
}
