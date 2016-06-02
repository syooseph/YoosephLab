#!/usr/bin/perl -w
use strict;

my $file1 = shift;
my $file2 = shift;
my $file3 = shift;

open my $IN1, "<$file1" or die "Cannot open file: $!\n";
open my $IN2, "<$file2" or die "Cannot open file: $!\n";
open my $IN3, "<$file3" or die "Cannot open file: $!\n";

while(<$IN1>)  {
  chomp;
  my $num1 = $_;
  my $num2 = <$IN2>;
  my $num3 = <$IN3>;
  chomp $num2;
  chomp $num3;
  if($num1 != $num2 + $num3)  {
    print "Error! $num1 $num2 $num3\n";
  }
}
close $IN1;
close $IN2;
close $IN3;
