#!/usr/bin/perl -w
use strict;

foreach(</usr/local/depot/projects/SPA/czhong/Grasp/Results/SGO*.gsp>)	{
	my $file = $_;
	$file =~ /.*\/(.*)\.gsp/;
	my $stem = $1;
	#print "$stem\n";
	open my $IN, "<$file" or die "Cannot open file: $!\n";
	while(<$IN>)	{
		chomp;
		if(/traversed:\s+(\d+)/)	{
			my $touched = $1;
			system "cat /usr/local/depot/projects/SPA/czhong/Grasp/Results/$stem.fa.0.reads | wc -l >temp";
			open my $TIN, "<temp" or die "Cannot open file: $!\n";
			my $num = <$TIN>;
			close $TIN;
			chomp $num;
			if($num > 1000)	{
				print "$stem	$touched	$num\n";
			}
			last;
		}
	}
	close $IN;
}
