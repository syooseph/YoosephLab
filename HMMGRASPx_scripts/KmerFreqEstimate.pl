#!/usr/bin/perl -w
use strict;

my $file = shift;	# estimating FASTA file
my $klen = shift;	# the k-mer length
my $cutoff = shift;	# the minimum frequency

my %kmer_hash;
open my $IN, "<$file" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  if(/^\>/)  {
    next;
  }
  my $seq = $_;
  for(my $i = 0; $i < length($seq) - $klen; ++ $i)  {
    my $sseq = substr($seq, $i, $klen);
    ++ $kmer_hash{$sseq};
  }
}
close $IN;

my %freq_hash;
my $above = 0;
foreach(keys %kmer_hash)  {
  if($kmer_hash{$_} <= $cutoff)  {
    ++ $freq_hash{$kmer_hash{$_}};
  }  else  {
    ++ $above;
  }
}

foreach(sort keys %freq_hash)  {
  print "$_	$freq_hash{$_}\n";
}
print "above:	$above\n";
