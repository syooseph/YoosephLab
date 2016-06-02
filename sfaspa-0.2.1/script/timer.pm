#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(floor fmod);
use Time::HiRes qw(gettimeofday tv_interval);
#use Time::localtime;

package timer;

sub printLocalTime
{
	my ( $sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst ) = localtime(time);
	print STDERR sprintf("[%d/%02d/%02d %02d:%02d:%02d]\n", 
						 $year+1900, $mon+1, $mday,
						 $hour, $min, $sec);
}

sub printElapsed
{
	my ( $ts, $task, $sum ) = @_;
	my $te = time();
	my $elapsed = $te-$ts;
	if ( $sum ) {
		print STDERR sprintf ("[Total %02.0f:%02.0f:%05.2f]\t%s\n", 
							  POSIX::floor($elapsed/3600.0),
							  POSIX::floor(POSIX::fmod($elapsed,3600.0)/60.0), 
							  POSIX::fmod($elapsed,60.0),
							  $task);
	} else {
		print STDERR sprintf ("[%02.0f:%02.0f:%05.2f]\t%s\n", 
							  POSIX::floor($elapsed/3600.0),
							  POSIX::floor(POSIX::fmod($elapsed,3600.0)/60.0), 
							  POSIX::fmod($elapsed,60.0),
							  $task);
	}
}


1;
