#!/usr/bin/perl -w
use strict;

my $asm_file = shift;
my $bst_file = shift;
my $asm_cutoff = shift;

my $blast_significant_cutoff = 1e-5;

# load all significant blast hits

my %good_paths;
open my $BIN, "<$bst_file" or die "Cannot open file: $!\n";
while(<$BIN>)	{
	chomp;
	if(/(\d+)\s+hits\s+found/)	{
		if($1 > 0)	{
			my $best_hit = <$BIN>;
			my @decom = split /\s+/, $best_hit;
			$good_paths{$decom[1]} = $decom[10];
		}
	}
}
close $BIN;

my %all_reads;
my %good_reads;
open my $AIN, "<$asm_file" or die "Cannot open file: $!\n";
while(<$AIN>)	{
	chomp;
	if(/^>/)	{
		chomp;
		my @decom = split /\s+/, $_;
		if($decom[1] > $asm_cutoff)	{
			next;
		}
		$decom[0] =~ s/^>//g;
		# record all read if the path is better than cutoff
		if(!(defined $decom[2]))	{
			next;
		}
		my @reads = split /\;/, $decom[2];
		foreach(@reads)	{
			if(/(\d+)\,/)	{
				my $rid = $1;
				if(exists $good_paths{$decom[0]} && $good_paths{$decom[0]} < $blast_significant_cutoff)	{
					$good_reads{$rid} = 1;
				}
				$all_reads{$rid} = 1;
			}
		}
	}
}
close $AIN;

my $num_good = scalar keys %good_reads;
my $num_all = scalar keys %all_reads;
print "$num_good	$num_all\n";
