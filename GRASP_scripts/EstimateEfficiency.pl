#!/usr/bin/perl -w
use strict;

my $file = shift;
my $cutoff = shift;

my $total_reads = 0;
my $identified_reads = 0;
my @saved_time;
my @total_time;

$saved_time[0] = $saved_time[1] = $saved_time[2] = 0;
$total_time[0] = $total_time[1] = $total_time[2] = 0;

open my $IN, "<$file" or die "Cannot open file: $!\n";
while(<$IN>)	{
	chomp;
	my $reads = $_;
	my $time = <$IN>;
	chomp $time;
	if($time =~ /^\[/ && $reads =~ /^\[/)	{
		$reads =~ /^\[(.*)\]/;
		my $reads_phase = $1;
		if($reads_phase >= $cutoff)	{
			$identified_reads += $reads_phase;
		}
		$total_reads += $reads_phase;
		$time =~ /^\[(.*)\]/;
		my @decom = split /\:/, $1;
		if($reads_phase < $cutoff)	{
			$saved_time[0] += $decom[0];
			$saved_time[1] += $decom[1];
			$saved_time[2] += $decom[2];
		}
		$total_time[0] += $decom[0];
		$total_time[1] += $decom[1];
		$total_time[2] += $decom[2];
	}	else	{
		last;
	}
}

print "$total_reads	$identified_reads\n";
print "$total_time[0]:$total_time[1]:$total_time[2]\n";
print "$saved_time[0]:$saved_time[1]:$saved_time[2]\n";
