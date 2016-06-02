#!/usr/bin/perl -w
use strict;

my $fasta = shift;	# the fasta file that needs to be subsampled
my $list = shift;	# the list of subsampling instruction; the first field is the regex
			# the second field is the subsamping rate (prob of keeping this read)

my %id_hash;			
open my $LIN, "<$list" or die "Cannot open file: $!\n";
while(<$LIN>)	{
	chomp;
	my @decom = split /\s+/, $_;
	$id_hash{$decom[0]} = $decom[1];
}
close $LIN;

open my $IN, "<$fasta" or die "Cannot open file: $!\n";
while(<$IN>)	{
	chomp;
	if(/^>/)	{
		my $line = $_;
		my $seq= <$IN>;
		my $ran_num = rand(1);
		my $rate = 0;
		for(keys %id_hash)	{
			my $k = $_;
			if($line =~ /$k/)	{
				$rate = $id_hash{$k};
			}
		}
		if($ran_num <= $rate)	{
			print "$line\n";
			print "$seq";
		}
	}
}
close $IN;
