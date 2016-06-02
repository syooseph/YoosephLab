#!/usr/bin/perl -w
use strict;

my $file = shift;

my @accumulate_count;
my $i;
for($i = 0; $i < 9; ++ $i) {
  $accumulate_count[$i] = 0;
}

open my $IN, "<$file" or die "Cannot open file: $!\n";
while(<$IN>) {
  chomp;
  my @decom = split /\s+/, $_;
  for($i = 0; $i < 9; ++ $i) {
    $accumulate_count[$i] += $decom[$i + 1];
  }
}
close $IN;

my $blast_full_sen = $accumulate_count[1] / $accumulate_count[0];
my $blast_full_spe = $accumulate_count[1] / $accumulate_count[2];
my $blast_par_cons_sen = $accumulate_count[3] / $accumulate_count[0];
my $blast_par_cons_spe = $accumulate_count[3] / $accumulate_count[4];
my $blast_par_sen = $accumulate_count[5] / $accumulate_count[0];
my $blast_par_spe = $accumulate_count[5] / $accumulate_count[6];
my $grasp_max_sen = $accumulate_count[7] / $accumulate_count[0];
my $grasp_max_spe = $accumulate_count[7] / $accumulate_count[8];
#my $grasp_min_sen = $accumulate_count[9] / $accumulate_count[0];
#my $grasp_min_spe = $accumulate_count[9] / $accumulate_count[10];

print "$blast_full_sen  $blast_full_spe $blast_par_cons_sen $blast_par_cons_spe $blast_par_sen  $blast_par_spe  $grasp_max_sen  $grasp_max_spe\n";
