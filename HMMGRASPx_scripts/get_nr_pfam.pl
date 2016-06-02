#!/usr/bin/perl -w
use strict;

my $folder = "/usr/local/projdata/0599/projects/SPA/czhong/Works/GRASPxp/Data/KEGG_info";

my %pfam_hash;
foreach(<$folder/*_pfam.list>)  {
  open my $IN, "<$_" or die "Cannot open file: $!\n";
  while(<$IN>)  {
    chomp;
    my @decom = split /\s+/, $_;
    $pfam_hash{$decom[1]} = 1;
  }
  close $IN;
}

foreach(sort keys %pfam_hash)  {
  print "$_\n";
}
