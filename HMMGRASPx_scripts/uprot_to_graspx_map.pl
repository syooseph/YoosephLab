#!/usr/bin/perl -w
use strict;

my $profile = shift;	# the multi-HMMER3 format input that contains families to be classified
my $file = shift;	# the uproc-prot output file

open my $PIN, "<$profile" or die "Cannot open file: $!\n";
my $name;
my $ID;
my %id_hash;
while(<$PIN>)  {
  chomp;
  if(/^NAME\s+(\S+)/)  {
    $name = $1;
  }  elsif(/^ACC\s+(PF\d+)/)  {
    $ID = $1;
    $id_hash{$ID} = $name;
  }
}
close $PIN;

open my $IN, "<$file"  or die "Cannot open file: $!\n";

while(<$IN>)  {
  chomp;
  my @decom = split /\,/, $_;
  #my @decom2 = split /\:\:/, $decom[3];
  print "XX	$decom[1]	CHECKED||$id_hash{$decom[3]}\n";
}
close $IN;
