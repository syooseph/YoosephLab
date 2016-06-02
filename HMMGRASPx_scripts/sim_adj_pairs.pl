#!/usr/bin/perl -w
use strict;

my $file = shift;
my @entries;
open my $IN, "<$file" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  push @entries, $_;
}
close $IN;

for(my $i = 0; $i < 1000; ++ $i)  {
  my @perm = @entries;
  for(my $j = 0; $j < 500; ++ $j)  {
    my $a = int(rand(scalar(@perm)));
    my $b = int(rand(scalar(@perm)));
    my $tmp = $perm[$a];
    $perm[$a] = $perm[$b];
    $perm[$b] = $tmp;
  }
  #print "@perm\n";
  my $adj = 0;
  for(my $k = 0; $k < scalar(@perm) - 1; ++ $k)  {
    if($perm[$k] eq $perm[$k + 1])  {
      ++ $adj;
    }
  }
  print "$adj\n";
}
