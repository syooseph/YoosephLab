#!/usr/bin/perl -w
use strict;

my $file = "/usr/local/projdata/0599/projects/SPA/czhong/GRASPxp_current/mix3fams.fa";
open my $IN, "<$file" or die "Cannot open file: $!\n";

while(<$IN>)  {
  chomp;
  if(/^>(\S+)/)  {
    my $line = $1;
    my @decom = split /\:\:/, $line;
    print "$decom[0]	$decom[1]	>$line\n";
    print "XX	$line	contig||$decom[1]\n";
  }
}
close $IN;
