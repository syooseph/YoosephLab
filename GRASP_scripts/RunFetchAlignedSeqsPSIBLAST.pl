#!/usr/bin/perl -w
use strict;

my $result_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/HMP_saliva/PSIBLAST/raw";
my $out_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/HMP_saliva/PSIBLAST/derived";
my $script_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Scripts";
my $db_file = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/saliva.faa";
foreach(<$result_dir/*.psibst>)  {
  my $file = $_;
  $file =~ /.*\/(.*)\.psibst/;
  my $file_stem = $1;
  print "Working on $file_stem...\n";
  system "perl $script_dir/FetchPSIBLASTAlignedSeqs.pl $file $db_file >$out_dir/$file_stem.asm.faa";
}
