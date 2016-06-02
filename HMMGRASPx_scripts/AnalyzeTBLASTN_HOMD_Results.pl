#!/usr/bin/perl -w
use strict;

# the script is used to analyze the tblastn results of searching GRASPxp contigs against the Human Oral Microbiome Database (HOMD)
# the purpose is to summarize the number of reads that are mapped to each of the strain/species

my $tblastn_results = shift;		# the TBLASTN results for searching GRASPxp contigs against HOMD
my $recruit_list = shift;		# the recruitmented reads computed by GRASPxp
my $homd_header = shift;		# the headers of the HOMD sequences (headers in the fasta with '>')

# get the first mapped entries for each contig
my %top_entry;
open my $IN, "<$tblastn_results" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  if(/\#\s+(\d+) hits found/ and $1 > 0)  {
    my $next = <$IN>;
    chomp $next;
    my @decom = split /\s+/, $next;
    $top_entry{$decom[0]} = $next;
  }
}
close $IN;

# count the reads mapped to each contig
my %count_hash;
open my $RIN, "<$recruit_list" or die "Cannot open file: $!\n";
while(<$RIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  if(exists $top_entry{$decom[2]})  {
    ++ $count_hash{$decom[2]};
  }
}
close $RIN;

# fetch the annotation from homd header
my %homd_annot;
open my $HIN, "<$homd_header" or die "Cannot open file: $!\n";
while(<$HIN>)  {
  chomp;
  my $line = $_;
  my @decom = split /\s+/, $line;
  $decom[0] =~ /^>(.*)/;
  $homd_annot{$1} = $line;
}
close $HIN;

# output the results
foreach(keys %top_entry)  {
  my @decom = split /\s+/, $top_entry{$_};
  if(exists $count_hash{$_} and $count_hash{$_} > 0)  {
    print "$_	$count_hash{$_}	$homd_annot{$decom[1]}\n";
  }
}
