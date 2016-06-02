#!/usr/bin/perl -w
use strict;

my $read_file = shift;
my $ref_file = shift;
my $id_file = shift;

my @ID;
open my $IN, "<$read_file" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  if(/^>(.*)/)  {
    push @ID, $1;
  }
}
close $IN;

my %ref;
open my $RIN, "<$ref_file" or die "Cannot open file: $!\n";
while(<$RIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  $ref{$decom[0]} = 1;
}
close $RIN;

my @fp;
my @tp;

open my $IIN, "<$id_file" or die "Cannot open file: $!\n";
my $temp = <$IIN>;
while(<$IIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  my $keyid = $ID[$decom[1]];
  if(exists $ref{$keyid})  {
    push @tp, $keyid;
    delete $ref{$keyid};
  }  else  {
    push @fp, $keyid;
  }
}
close $IIN;

print "***************TP****************\n";
foreach(sort @tp)  {
  print "$_\n";
}

print "***************FP****************\n";
foreach(sort @fp)  {
  print "$_\n";
}


print "***************FN****************\n";
foreach(sort keys %ref)  {
  print "$_\n";
}
