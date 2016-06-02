#!/usr/bin/env perl

################################################################################
#
#
#
################################################################################

use strict;
use warnings;

package io;



##==============================================================================
## Justin's open module
##==============================================================================
sub openInput
{
	my ($fileName) = @_;
	return *STDIN unless defined $fileName;

	my ($fd);
	open($fd, 
		 $fileName =~ /.gz(ip)?$/  ? "zcat $fileName |"  : 
		 $fileName =~ /.bz(ip)?2$/ ? "bzcat $fileName |" : 
		 $fileName) || die("Open error: $fileName");
	return $fd;
}

sub openOutput
{
    my ($fileName) = @_;
    return *STDOUT unless defined $fileName;
    
    my ($fd);
    open($fd, 
		 $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : 
		 $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : 
		 ">$fileName") || die("Open error: $fileName");
    return $fd;
}
##==============================================================================

1;
