#!/usr/bin/perl -w
use strict;

my %Grasp_hash;
my %Blast_hash;

for(</usr/local/depot/projects/SPA/czhong/Grasp/Results/LBPG*.aln>)	{
	open my $IN, "<$_" or die "Cannot open file: $!\n";
	while(<$IN>)	{
		chomp;
		if(/ALIGNMENT/)	{
			my $stat = <$IN>;
			$stat = <$IN>;
			chomp $stat;
			$stat =~ /Evalue:\s+(.*)/;
			my $evalue = $1;
			#print "$evalue\n";
			my $line = "";
			while(!($line =~ /^included/))	{
				$line = <$IN>;
			} 
			$line = <$IN>;
			chomp $line;
			my @decom = split /\;/, $line;
			foreach(@decom)	{
				#print "$_\n";
				if(/(.*)\,/)	{
					if(!(exists $Grasp_hash{$1}) || $Grasp_hash{$1} > $evalue)	{
						$Grasp_hash{$1} = $evalue;
						#print "$1	$evalue\n";
					}
				}
			}
		}
	}
	close $IN;
}

for(</usr/local/depot/projects/SPA/czhong/Grasp/Results/LBPG*.bst>)	{
	if(/reBLAST/)	{
		next;
	}
	open my $IN, "<$_" or die "Cannot open file: $!\n"; 
	while(<$IN>)	{
		chomp;
		if(/^\#/)	{
			next;
		}
		my @decom = split /\s+/, $_;
		if(!(exists $Blast_hash{$decom[1]}) || $Blast_hash{$decom[1]} > $decom[10])	{
			$Blast_hash{$decom[1]} = $decom[10];
		}
	}
	close $IN;	
}

my $i;
for($i = 1; $i >= -10; -- $i)	{
	my $cutoff = "1e" . $i;
	my $num_ghit = 0;
	my $num_bhit = 0;
	foreach(keys %Grasp_hash)	{
		if($Grasp_hash{$_} < $cutoff)	{
			++ $num_ghit;
		}
	}
	foreach(keys %Blast_hash)	{
		if($Blast_hash{$_} < $cutoff)	{
			++ $num_bhit;
		}
	}
	my $fold_change = $num_ghit / $num_bhit;
	print "$cutoff	$num_ghit	$num_bhit	$fold_change\n";
}
