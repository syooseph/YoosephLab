#!/usr/bin/perl -w
use strict;

my $result_file = shift;	# expecting tabulated blast result file
my $db_file = shift;	# the multi-fasta file which was used as the db for the search

# read file and construct hash
my %header_hash;
my %seq_hash;
my $n = 0;
open my $DIN, "<$db_file" or die "Cannot open file: $!\n";
while(<$DIN>)  {
  chomp;
  if(/^>(.*)/)  {
    $header_hash{$1} = $n;
    ++ $n;
    my $seq = <$DIN>;
    chomp $seq;
    $seq_hash{$1} = $seq;
  }
}
close $DIN;

$n = 0;
open my $RIN, "<$result_file" or die "Cannot open file: $!\n";
while(<$RIN>)  {
  chomp;
  if(/^#/)  {
    next;
  }
  my @decom = split /\s+/, $_;
  if(scalar(@decom) >= 11 && exists $header_hash{$decom[1]})  {
    print ">assembled_sequence:$n	$decom[10]	$header_hash{$decom[1]},0;\n";
    print "$seq_hash{$decom[1]}\n";
    ++ $n;
  }
}
close $RIN;
