#!/usr/bin/perl -w
use strict;

chdir "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/marine_sim_v2/Coreset" or die "Cannot cnahge directory: $!";

foreach(<*.bst>)	{
	my $file_name = $_;
	$_ =~ /(.*)\.bst/;
	my $id = $1;
	print "Running $id\n";
	system "perl /usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Scripts/ListIncludedReads.pl $file_name /usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/marine.core.faa /usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/marine_sim_v2/read_remap.csv >/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/marine_sim_v2/Reference/$id.ref";
}
