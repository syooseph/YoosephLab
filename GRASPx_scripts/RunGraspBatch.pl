#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';

my $fa_file = shift; # the multi-fasta file that needs to be run
my $target_file = shift; # the file contains the reads
$target_file = abs_path($target_file);

open my $IN, "<$fa_file" or die "Cannot open file: $!\n";

chdir "/usr/local/depot/projects/SPA/czhong/Grasp/" or die "Cannot change directory: $!\n";

while(<$IN>){
    if(/^>/){
	my @name = split /\s+/, $_;
	$name[0] =~ s/^>//g;
	# write the sequence into a separate file
	open my $OUT, ">./WorkSpace/$name[0]\.fa" or die "Cannot create file: $!\n";
	print $OUT ">$name[0]\n";
	my $seq = <$IN>;
	chomp $seq;
	print $OUT "$seq\n";
	close $OUT;
	# run grasp on the query
	print "running $name[0]...\n";
	system "/usr/bin/time -v ./Search ./WorkSpace/$name[0]\.fa $target_file --write_certificate 1 --num_threads 4";
    }
}
