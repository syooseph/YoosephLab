#!/usr/bin/perl -w
use strict;

my $profile = shift;	# the multi-HMMER3 format input that contains families to be classified
my $folder = shift;	# the salt output file

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

foreach(<$folder/*/*.out>)  {
  open my $IN, "<$_"  or die "Cannot open file: $!\n";
  my $fam;
  while(<$IN>)  {
    chomp;
    if(/^>(.*)/)  {
      $fam = $id_hash{$1};
    }  else  {
      if(/ref/) {
        print "XX	$_	CHECKED||$fam\n";
      }
    }
  }
  close $IN;
}
