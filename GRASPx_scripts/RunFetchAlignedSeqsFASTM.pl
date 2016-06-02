#!/usr/bin/perl -w
use strict;

my $result_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/HMP_saliva/FASTM/raw";
my $out_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/HMP_saliva/FASTM/derived";
my $script_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Scripts";
my $db_file = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/saliva.faa";
foreach(<$result_dir/*.fastm>)  {
  my $file = $_;
  $file =~ /.*\/(.*)\.fastm/;
  my $file_stem = $1;
  system "perl $script_dir/FetchFASTMAlignedSeqs.pl $file $db_file >$out_dir/$file_stem.asm.faa";
}
