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
  if(/^\>\>/)  {
    $_ =~ /^>>(\S+)/;
    my $name = $1;
    #print "$name\n";
    #print "$_\n";
    $_ = <$RIN>;
    chomp;
    #print "$_\n";
    my @decom = split /\s+/, $_;
    if(exists $header_hash{$name})  {
      print ">assembled_sequence:$n	$decom[12]	$header_hash{$name},0;\n";
      print "$seq_hash{$name}\n";
      ++ $n;
    }
  }
}
close $RIN;
