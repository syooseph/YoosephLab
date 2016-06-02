#!/usr/bin/perl -w
use strict;

my $spe_tag = shift;
my $folder = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/HMP_saliva";

my @performances;
open my $GIN, "<$folder/GRASP/performance/$spe_tag.GRASP.perf.txt" or die  "Cannot open file: $!\n";
my $n = 0;
while(<$GIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  push @{$performances[$n]}, $decom[1];
  push @{$performances[$n]}, $decom[2];
  ++ $n;
}
close $GIN;

open my $FIN, "<$folder/FASTM/performance/$spe_tag.FASTM.perf.txt" or die  "Cannot open file: $!\n";
$n = 0;
while(<$FIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  push @{$performances[$n]}, $decom[1];
  push @{$performances[$n]}, $decom[2];
  ++ $n;
}
close $FIN;

open my $PIN, "<$folder/PSIBLAST/performance/$spe_tag.PSIBLAST.perf.txt" or die  "Cannot open file: $!\n";
$n = 0;
while(<$PIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  push @{$performances[$n]}, $decom[1];
  push @{$performances[$n]}, $decom[2];
  ++ $n;
}
close $PIN;

open my $BIN, "<$folder/BLASTP/performance/$spe_tag.BLASTP.perf.txt" or die  "Cannot open file: $!\n";
$n = 0;
while(<$BIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  push @{$performances[$n]}, $decom[1];
  push @{$performances[$n]}, $decom[2];
  ++ $n;
}
close $BIN;

foreach(@performances)  {
  print "@{$_}\n";
}
