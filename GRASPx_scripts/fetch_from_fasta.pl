#!/usr/bin/perl -w
use strict;

my $fasta_file = shift;  # the parsed GENBANK file in multi-fasta format
my $id_file = shift;  # a list of IDs that are to be retrieved

open my $IN, "<$id_file" or die "Cannot open file: $!\n";
my @to_pick;
while(<$IN>){
  chomp;
  my @decom = split /\s+/, $_;  
  push @to_pick, $decom[0];
}
close $IN;


my @contents;
open my $FIN, "<$fasta_file" or die "Cannot open file: $!\n";
while(<$FIN>){
  chomp;
  push @contents, $_;
}
my $i = 0;
my $j = 0;
for($i = 0; $i < scalar(@contents); ++ $i)  {
  for($j = 0; $j < scalar(@to_pick); ++ $j){
    #print "$contents[$i]	$to_pick[$j]\n";
    if($contents[$i] =~ /^>/ && $contents[$i] =~ /$to_pick[$j]\s+/){
      my @decom = split /\s+/, $contents[$i];
      print "$decom[0]\n";
      my $seq = "";
      ++ $i;
      while(!($contents[$i] =~ /^>/) && $i < scalar(@contents))  {
        $seq .= $contents[$i];
        ++ $i;
      }
      -- $i;
      print "$seq\n";
    }
  }
}
