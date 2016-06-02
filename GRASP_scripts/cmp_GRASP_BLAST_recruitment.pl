#!/usr/bin/perl -w
use strict;

my $bcutoff = shift;
my $gcutoff = shift;

chdir "/usr/local/depot/projects/SPA/czhong/Grasp/Results" or die "Cannot change directory: $!\n";

foreach(<SGO_*.bst>){
    if(!($_ =~ /SGO\_\d+\.bst/))    {
	next;
    }
    #print "$_\n";
    my $bfile = $_;
    $bfile =~ /(.*)\.bst/;
    my $query_name = $1;
    my $gfile = $1 . "\.fa\.0\.reads";
    if(! -e "$gfile")	{
        next;
    }
    #print "$bfile\t$gfile\n";
    # loop over the blast file
    my $num_blast_hit = 0;
    open my $BIN, "<$bfile" or die "Cannot open file: $!\n";
    while(<$BIN>){
	chomp;
	if(/^\#/){
	    next;
	}
	my @decom = split /\s+/, $_;
	if($decom[10] <= $bcutoff){
	    ++ $num_blast_hit;
	}
    }
    close $BIN;
    # loop over the grasp file
    my $num_grasp_hit = 0;
    open my $GIN, "<$gfile" or die "Cannot open file: $!\n";
    while(<$GIN>){
	chomp;
	if(/^\#/){
	    next;
	}
	my @decom = split /\s+/, $_;
	if($decom[2] <= $gcutoff){
	    ++ $num_grasp_hit;
	}
    }
    close $GIN;
    print "$query_name\t$num_grasp_hit\t$num_blast_hit\n";
}
