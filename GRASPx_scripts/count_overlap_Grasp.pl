#!/usr/bin/perl -w
use strict;

my $file_1 = shift;
my $file_2 = shift;

my %seq_IDs;

open my $IN1, "<$file_1" or die "Cannot open file: $!\n";
my $num_entries_1 = 0;
my $num_entries_2 = 0;
my $num_overlap = 0;
while(<$IN1>){
    chomp;
    my @decom = split /\s+/, $_;
    $seq_IDs{$decom[1]} = 1;
    ++ $num_entries_1;
}
close $IN1;

open my $IN2, "<$file_2" or die "Cannot open file: $1\n";
while(<$IN2>){
    chomp;
    my @decom = split /\s+/, $_;
    if(exists $seq_IDs{$decom[1]}){
	++ $num_overlap;
    }
    ++ $num_entries_2;
}

print "$num_overlap\t$num_entries_1\t$num_entries_2\t$file_1\t$file_2\n";
