#!/usr/bin/perl -w
use strict;

my $file = shift;

my $blast_full_sen = 0;
my $blast_full_spe = 0;
my $blast_par_cons_sen = 0;
my $blast_par_cons_spe = 0;
my $blast_par_sen = 0;
my $blast_par_spe = 0;
my $grasp_sen = 0;
my $grasp_spe = 0;
my $num_count = 0;

open my $IN, "<$file" or die "Cannot open file: $!\n";
while(<$IN>) {
  chomp;
  my @decom = split /\s+/, $_;
  if($decom[1] <= 0)  {
    $blast_full_sen += 1;
    $blast_par_cons_sen += 1;
    $blast_par_sen += 1;
    $grasp_sen += 1;
    #$blast_full_spe += 1;
    #$blast_par_cons_spe += 1;
    #$blast_par_sen += 1;
    #$graps_sen += 1;
    if($decom[3] <= 0)  {
      $blast_full_spe += 1;
    }
    if($decom[5] <= 0)  {
      $blast_par_cons_spe += 1;
    }
    if($decom[7] <= 0)  {
      $blast_par_spe += 1;
    }
    if($decom[9] <= 0)  {
      $grasp_spe += 1;
    }
    ++ $num_count;
    next;
  }
  
  if($decom[3] > 0)  {
    $blast_full_sen += $decom[2] / $decom[1];
    $blast_full_spe += $decom[2] / $decom[3];
  }  else  {
    $blast_full_sen += $decom[2] / $decom[1];
    $blast_full_spe += 1;
  }
  if($decom[5] > 0)  {
    $blast_par_cons_sen += $decom[4] / $decom[1];
    $blast_par_cons_spe += $decom[4] / $decom[5];
  }  else  {
    $blast_par_cons_sen += $decom[4] / $decom[1];
    $blast_par_cons_spe += 1;
  }
  if($decom[7] > 0)  {
    $blast_par_sen += $decom[6] / $decom[1];
    $blast_par_spe += $decom[6] / $decom[7];
  }  else  {
    $blast_par_sen += $decom[6] / $decom[1];
    $blast_par_spe += 1;
  }
  if($decom[9] > 0)  {
    $grasp_sen += $decom[8] / $decom[1];
    $grasp_spe += $decom[8] / $decom[9];
  }  else  {
    $grasp_sen += $decom[8] / $decom[1];
    $grasp_spe += 1;
  }
  ++ $num_count;
}
close $IN;

$blast_full_sen /= $num_count;
$blast_full_spe /= $num_count;
$blast_par_cons_sen /= $num_count;
$blast_par_cons_spe /= $num_count;
$blast_par_sen /= $num_count;
$blast_par_spe /= $num_count;
$grasp_sen /= $num_count;
$grasp_spe /= $num_count;


#my $blast_full_sen = $accumulate_count[1] / $accumulate_count[0];
#my $blast_full_spe = $accumulate_count[1] / $accumulate_count[2];
#my $blast_par_cons_sen = $accumulate_count[3] / $accumulate_count[0];
#my $blast_par_cons_spe = $accumulate_count[3] / $accumulate_count[4];
#my $blast_par_sen = $accumulate_count[5] / $accumulate_count[0];
#my $blast_par_spe = $accumulate_count[5] / $accumulate_count[6];
#my $grasp_max_sen = $accumulate_count[7] / $accumulate_count[0];
#my $grasp_max_spe = $accumulate_count[7] / $accumulate_count[8];
#my $grasp_min_sen = $accumulate_count[9] / $accumulate_count[0];
#my $grasp_min_spe = $accumulate_count[9] / $accumulate_count[10];

print "$blast_full_sen  $blast_full_spe $blast_par_cons_sen $blast_par_cons_spe $blast_par_sen  $blast_par_spe  $grasp_sen  $grasp_spe\n";
