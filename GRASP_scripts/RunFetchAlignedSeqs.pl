#!/usr/bin/perl -w
use strict;

my $result_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/HMP_saliva/BLASTP/raw";
my $out_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/HMP_saliva/BLASTP/derived";
my $script_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Scripts";
my $db_file = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/saliva.faa";
foreach(<$result_dir/*.bst>)  {
  my $file = $_;
  $file =~ /.*\/(.*)\.bst/;
  my $file_stem = $1;
  system "perl $script_dir/FetchBLASTAlignedSeqs.pl $file $db_file >$out_dir/$file_stem.asm.faa";
}
