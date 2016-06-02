#!/usr/bin/perl -w
use strict;

my $file = shift;	# expecting the multi-pfam hmm file

# printing the Pfam_ID, length, and name information

open my $IN, "<$file" or die "Cannot open file: $!\n";
my $ID;
my $len;
my $name;
my $desc;
while(<$IN>)  {
  chomp;
  if(/^NAME\s+(\S+)/)  {
    $name = $1;
  }  elsif(/^ACC\s+(RF\d+)/)  {
    $ID = $1;
  }  elsif(/^DESC\s+(.*)/)  {
    $desc = $1;
    #print "$ID	$len	$name\n";
  }  elsif(/^LENG\s+(\d+)/)  {
    $len = $1;
    print "$ID	$len	$name	$desc\n";
  }
}
close $IN;
