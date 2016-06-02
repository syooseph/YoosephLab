#!/usr/bin/perl -w
use strict;

my $file = shift;
my $folder = shift;

open my $IN, "<$file" or die "Cannot open file: $!\n";
while(<$IN>)	{
	chomp;
	if(/^>/)	{
		my $line = $_;
		my @decom = split /\s+/, $_;
		my $id = $decom[0];
		$id =~ s/\>//g; $id =~ s/\|/\_/g;
		open my $OUT, ">$folder/$id.fasta" or die "Cannot create file: $!\n";
		print $OUT ">$line\n";
		my $seq = <$IN>;
		print $OUT "$seq";
	}
}
