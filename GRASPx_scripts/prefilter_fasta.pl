#!/usr/bin/perl -w
use strict;

# the script is used to clean up the header of the fasta file
my $file = shift;
open my $IN, "<$file" or die "Cannot open file: $!\n";
my @contents;

while(<$IN>)  {
  chomp;
  push @contents, $_;
}
close $IN;

for(my $i = 0; $i < scalar(@contents); ++ $i)  {
  if($contents[$i] =~ /^>/)  {
    my @decom = split /\s+/, $contents[$i];
    print "$decom[0]\n";
    ++ $i;
    my $seq = "";
    while($i < scalar(@contents) && !($contents[$i] =~ /^>/))  {
      $seq = $seq . $contents[$i];
      ++ $i;
    }
    print "$seq\n";
    -- $i;
  }
}
