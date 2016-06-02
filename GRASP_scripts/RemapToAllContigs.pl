#!/usr/bin/perl -w
use strict;

my $script_folder = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Scripts";
my $certificate_folder = "/usr/local/depot/projects/SPA/czhong/Grasp_stable/Amphora2_marine/";
my $output_folder = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/marine_sim_amphora2/Grasp_remap";

my $remap_exe = "/usr/local/depot/projects/SPA/czhong/GRASPx_current/grasp-map";
my $db_file = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/marine.core.faa";
my $workspace = "/usr/local/depot/projects/SPA/czhong/GRASPx_current/WorkSpace";

foreach(<$certificate_folder/*.fa.aln>)  {
  my $aln_file = $_;
  $aln_file =~ /.*\/(.*)\.fa\.aln/;
  my $stem = $1;
  print "Working on $stem...\n";
  # convert aln file to contigs
  system "perl $script_folder/ConvertCertificateToMFASTA.pl $aln_file >$output_folder/$stem.asm.faa";
  # map the individual reads to the contigs
  system "$remap_exe $db_file $output_folder/$stem.asm.faa $output_folder/$stem.remap --work_space=$workspace --num_errors 3 --portion_mapped 0.6";
  # convert the remap results to GRASP output format
  system "perl $script_folder/RemapToGRASPOutput.pl $output_folder/$stem.remap >$output_folder/$stem.fa.reads";
}
