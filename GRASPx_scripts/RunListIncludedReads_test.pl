#!/usr/bin/perl -w
use strict;

my $full_blast_results = shift;
my $reads_set = shift;
my $map_file = shift;
my $out_folder = shift;

foreach(<$full_blast_results/*.bst>) {
  my $blast_file = $_;
  $blast_file =~ /(.*)\/(.*)\.bst/;
  my $file_stem = $2;
  print "Working on $file_stem...\n";
  system "perl /usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Scripts/ListIncludedReads_test.pl $blast_file $reads_set $map_file >$out_folder/$file_stem.ref";
}
