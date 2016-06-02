#!/usr/bin/env perl

# Created : Jul 28 15:49:27 EDT 2009
# Modified: Wed Jan  2 14:00:59 EST 2013

# Youngik Yang
#

use strict;
use warnings;
use Getopt::Long;

my ($input, $clean, $help);
my ($id, $seq) = (undef, "");
GetOptions( "input=s" => \$input,
			"clean" => \$clean,
			"help" => \$help );
usage() if defined $help;

my $in = openInput( $input );
while (<$in>) {
	#chomp;
	s/^\s*//g; ## remove any leading whitespaces
	s/\s*$//g; ## remove any trailing whitespaces
	#next if /^\s*$/; ## no need
	if ( /^>/ ) {
		&flush if defined $id;
		$id = $';
	} else {
		$seq .= $_;
	}
}
&flush;

sub flush {
	$seq =~ s/[^A-Za-z]+//g if defined $clean;
	print join("\t", $id, length($seq)), "\n";
	$seq = "";
}

## Justin's open module
sub openInput
{
        my ($fileName) = @_;

        return *STDIN unless defined $fileName;

        my ($fd);
        open($fd, $fileName =~ /.gz(ip)?$/ ? "zcat $fileName |" : $fileName =~ /.bz(ip)?2$/ ? "bzcat $fileName |" 
			 : $fileName) || die("Open error: $fileName");
        return $fd;
}

sub usage
{
	print "Usage:$0 [options]\n";
	print "      -i (--input) input-fasta\n";
	print "      -c (--clean) alphabet only (no gap, no stop codon)\n";
	print "      -h (--help) print this message\n";
	exit;
}
