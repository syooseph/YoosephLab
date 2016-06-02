#!/usr/bin/perl -w
use strict;

my $ko_map_file = shift;
my $psi_file = shift;
my $cutoff = shift;

sub ParseKEGGKOFile($)  {
  my $file = shift;
  my %g2k_hash; 
  open my $IN, "<$file" or die "Cannot open file: $!\n";
  while(<$IN>)  {
    chomp;
    my @decom = split /\s+/, $_;
    $g2k_hash{$decom[0]} = $decom[1];
  }
  close $IN;
  return %g2k_hash;
}

my %g2k_hash = ParseKEGGKOFile($ko_map_file);
my %map_hash_evalue;
my %map_hash_dest;

# go through the psi-blast results
open my $IN, "<$psi_file" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  if(/^#/ || /^\//)  {
    next;
  }
  my @decom = split /\s+/, $_;
  if(scalar(@decom) >= 11 && $decom[10] <= $cutoff && exists $g2k_hash{$decom[0]})  {
    if(!(exists $map_hash_dest{$decom[1]}) || $map_hash_evalue{$decom[1]} > $decom[10])  {
      $map_hash_evalue{$decom[1]} = $decom[10];
      $g2k_hash{$decom[0]} =~ /(\d+)/;
      $map_hash_dest{$decom[1]} = "K" . $1;
    }
    #$g2k_hash{$decom[0]} =~ /(\d+)/;
    #print "0	K$1\n";
  }
}
close $IN;

foreach(keys %map_hash_dest)  {
  print "$_	$map_hash_dest{$_}\n";
}
