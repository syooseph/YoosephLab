#!/usr/bin/perl -w
use strict;

my @contents;
foreach(</usr/local/depot/projects/SPA/czhong/Grasp/Results/SGO_*.reads>){
    push @contents, $_;
}

my $i;
my %id_hash;
for($i = 0; $i < scalar(@contents); ++ $i){
    open my $IN, "<$contents[$i]" or die "Cannot open file: $!\n";
    while(<$IN>){
	chomp;
	if(/^\#/){
	    next;
	}
	my @decom = split /\s+/, $_;
	$id_hash{$decom[1]} = 1;
    }
    close $IN;
}

my $a = scalar keys %id_hash;
print "$a\n";
