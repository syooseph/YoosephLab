#!/usr/bin/perl -w
use strict;

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
  $decom[0] =~ /(^PF\d+)\./;
  $of_interest_ID = $1;
}
close $IN;

my $num_reads = 0;
open my $RIN, "<$read_file" or die "Cannot open read file: $!\n";
my %included_reads;
while(<$RIN>) {
  chomp;
  my $full_header = $_;
  if($full_header =~ /^>/)  {
    if($full_header =~ /$of_interest_ID/)  {
      $full_header =~ s/^>//g;
      print "$full_header	$num_reads\n";
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
