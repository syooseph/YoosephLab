#!/usr/bin/perl -w
use strict;

my %first_hash;
my %second_hash;

for(</usr/local/depot/projects/SPA/czhong/Grasp/Results/SGO*.bst>)     {
        if(/reBLAST/)   {
                next;
        }
        open my $IN, "<$_" or die "Cannot open file: $!\n";
        while(<$IN>)    {
                chomp;
                if(/^\#/)       {
                        next;
                }
                my @decom = split /\s+/, $_;
                if(!(exists $first_hash{$decom[1]}) || $first_hash{$decom[1]} > $decom[10])     {
                        $first_hash{$decom[1]} = $decom[10];
                }
        }
        close $IN;
}

for(</usr/local/depot/projects/SPA/czhong/Grasp/Results/LBPG*.bst>)     {
        if(/reBLAST/)   {
                next;
        }
        open my $IN, "<$_" or die "Cannot open file: $!\n";
        while(<$IN>)    {
                chomp;
                if(/^\#/)       {
                        next;
                }
                my @decom = split /\s+/, $_;
                if(!(exists $second_hash{$decom[1]}) || $second_hash{$decom[1]} > $decom[10])     {
                        $second_hash{$decom[1]} = $decom[10];
                }
        }
        close $IN;
}

my $i;
for($i = 1; $i >= -10; -- $i)   {
        my $cutoff = "1e" . $i;
        my $num_hit1 = 0;
        my $num_hit2 = 0;
        my $num_overlap = 0;
        foreach(keys %first_hash)       {
                if($first_hash{$_} < $cutoff)   {
                        ++ $num_hit1;
                        if(exists $second_hash{$_} && $second_hash{$_} < $cutoff)       {
                                ++ $num_overlap;
                        }
                }
        }
        foreach(keys %second_hash)       {
                if($second_hash{$_} < $cutoff)   {
                        ++ $num_hit2;
                }
        }
        print "$cutoff  $num_overlap    $num_hit1       $num_hit2\n";
}
