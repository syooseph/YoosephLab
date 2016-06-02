#!/usr/bin/perl -w
use strict;

foreach(</usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/marine_sim_glycolysis/Fastm/deh*.fastm>)  {
	my $old_name = $_;
	my $new_name = $old_name;
	$new_name =~ s/deh\_/deh\:/g;
	system "mv $old_name $new_name";
}
