#!/usr/bin/perl -w
use strict;

foreach(</usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/SGO*.fa>)	{
	my $file = $_;
	$file =~ /.*\/(.*)/;
	my $stem = $1;
	system "/usr/local/depot/projects/SPA/czhong/Grasp/c.out $file 30 >/usr/local/depot/projects/SPA/czhong/Grasp/Results/$stem.sel.sum";
}
