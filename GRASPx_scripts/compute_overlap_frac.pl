#!/usr/bin/perl -w
use strict;

my $file = shift;

open my $IN, "<$file" or die "Cannot open file: $!\n";
while(<$IN>){
    chomp;
    my @decom = split /\s+/, $_;
    my $frac = $decom[0] * 2 / ($decom[1] + $decom[2]);
    print "$frac\t$decom[3]\t$decom[4]\n"
}
close $IN;
