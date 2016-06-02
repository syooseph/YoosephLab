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
my @contents;
while(<$SEQIN>)	{
	chomp;
	push @contents, $_;
}
close $SEQIN;

my $i;
my $currentID;
my $seq;
for($i = 0; $i < scalar(@contents); ++ $i)	{
    if($contents[$i] =~ /^>/)  {
	my @decom = split /\s+/, $contents[$i];
	$decom[0] =~ s/^>//g;
	$currentID = $decom[0];
	if(exists $lookup_ID{$currentID})  {
		++ $i;
		$seq = "";
		while(!($contents[$i] =~ /^>/))	{
			$seq .= $contents[$i];
			++ $i;
		}
		-- $i;
		print ">$currentID\n";
		print "$seq\n";
	}
    }
}
close $SEQIN;

