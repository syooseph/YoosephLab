#!/usr/bin/perl -w
use strict;

my $file_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/Mix3Fams/Psi";

foreach(<$file_dir/*.psibst>)  {
  $_ =~ /.*\:\:(.*)/;
  system "mv $_ $file_dir/$1";
}
