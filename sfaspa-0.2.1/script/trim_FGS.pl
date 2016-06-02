#!/usr/bin/env perl

BEGIN { push @INC, `echo -n \$SPA_SCRIPT` }

use strict;
use warnings;
use Set::Scalar;
use io;
use Getopt::Long;

my ( $fgs_seq, $ids_file, $new_seq ) = ( undef, undef, undef, undef );
GetOptions( "fgs=s" => \$fgs_seq,
			"ids=s" => \$ids_file,
			"new=s" => \$new_seq );

if ( !defined $fgs_seq || !defined $ids_file ) {
	print "[Error] FGS seq file and ID files must be given\n"; 
	usage();
}


my $ifh = io::openInput($ids_file);
my @Ids = <$ifh>; chomp @Ids;
close $ifh;

my $Target = Set::Scalar->new(@Ids);

my ( $seq, $id );
my $sfh = io::openInput($fgs_seq);
my $ofh = io::openOutput($new_seq);
while ( <$sfh> ) {
	chomp;
	if ( /^>([^\s]+)/ ) {
		if ( defined $id ) {
			$id =~ /_\d+_\d+_[+-]$/;
			my $sid = $`;
			print $ofh join("\n", ">$id", $seq), "\n" if defined $seq && ! $Target->has($sid);
		}
		$id = $1; $seq = "";
	} else {
		$seq .= $_;
	}
}
if ( defined $id ) {
	$id =~ /_\d+_\d+_[+-]$/;
	my $sid = $`;
	print $ofh join("\n", ">$id", $seq), "\n" if defined $seq && ! $Target->has($sid);
}

close($sfh);

sub usage
{
	print "Usage:$0 -f fgs-seq -i ids-file [-n new-seq]\n"; exit(1);
}

