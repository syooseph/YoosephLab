#!/usr/bin/env perl

BEGIN { push @INC, `echo -n \$SPA` }

use strict;
use warnings;
use Set::Scalar;
use io;
use Getopt::Long;

my $fgs_out;
GetOptions( "fgs=s" => \$fgs_out );
if ( ! defined $fgs_out ) {
	print "Usage:$0 -f fgs-out\n";
	exit(1);
}

my $fh = io::openInput($fgs_out);
my ( $id, $indel ) = ( undef, 0 );
while (<$fh>) {
	chomp;
	if ( /^>([^\s]+)/ ) {
		print $id, "\n" if defined $id && $indel;
		$id = $1;
		$indel = 0;
	}
	#1       100     +       1       1.362077        I:54,   D:
	else {
		$indel = 1 if /I:\d+/;
		$indel = 1 if /D:\d+/;
	}
}
print $id, "\n" if defined $id && $indel;
