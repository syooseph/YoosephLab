#!/usr/bin/perl -w
use strict;

my $hmm_result_folder = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/core_marine_glycolysis_annot_hmm";
my $ref_file_folder = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/Marine_Core_Genomes/Annotation";

# build hash for target genomes
my %genome_id;
foreach(<$hmm_result_folder/*.hmm.search>)  {
  $_ =~ /(T\d+\-NC\_\d+)\.hmm\.search/;
  $genome_id{$1} = 1;
}

foreach(keys %genome_id)  {
  my $gid = $_;
  # load the annotation file for this genome
  my %annot_hash;
  open my $RIN, "<$ref_file_folder/$gid.kff" or die "Cannot open file: $!\n";
  while(<$RIN>)  {
    chomp;
    my @decom = split /\s+/, $_;
    $annot_hash{$decom[0]} = $decom[4];
  }
  close $RIN;
  # parse the hmm result file
  foreach(<$hmm_result_folder/*$gid*.hmm.search>)  {
    $_ =~ /.*\/(.*)\.hmm\.search/;
    my $search_stem = $1;
    open my $HIN, "<$_" or die "Cannot open file: $!\n";
    my $n = 0;
    while(<$HIN>)  {
      ++ $n;
      if($n == 15)  {
        chomp;
        my @decom = split /\s+/, $_;
        if(scalar(@decom) >= 11 && $decom[1] < 10e-10)  {
           print "$search_stem	$decom[9]	$annot_hash{$decom[9]}\n";
        }
        last;
      }
    }
    close $HIN;
  }
}
