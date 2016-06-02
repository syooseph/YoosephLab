#!/usr/bin/perl -w
use strict;

my $id_list = shift;  # should be a list of IDs that are expected to present in the header of the fasta format
my $target_file = shift;  # should be a multi-fasta file that is expected to contain the sequences

open my $IN, "<$id_list" or die "Cannot open ID file: $!\n";
my %lookup_ID;
while(<$IN>)  {
    chomp;
    my @decom = split /\s+/, $_;
    $lookup_ID{$decom[0]} = 1;
}
close $IN;

open my $SEQIN, "<$target_file" or die "Cannot open target sequence file: $!\n";
my $currentID = 'foo';
my $seq = '';
while(<$SEQIN>) {
    chomp;
    if(/^>/)  {
	if(exists $lookup_ID{$currentID})  {
	    print ">$currentID\n";
	    print "$seq\n";
	}
	if(/\|(.*)\|/)  {
	    $currentID = $1;
	    $seq = '';
	} else  {
	    $currentID = 'foo';
	    $seq = '';
	}
    } else  {
	$seq .= $_;
    }
}
close $SEQIN;

