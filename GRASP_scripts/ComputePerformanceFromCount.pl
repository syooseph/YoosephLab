#!/usr/bin/perl -w
use strict;

my $file = shift;

my $blast_par_sen = 0;
my $blast_par_spe = 0;
my $psiblast_par_sen = 0;
my $psiblast_par_spe = 0;
my $fastm_sen = 0;
my $fastm_spe = 0;
my $grasp_sen = 0;
my $grasp_spe = 0;
my $num_count = 0;

open my $IN, "<$file" or die "Cannot open file: $!\n";
while(<$IN>) {
  chomp;
  my @decom = split /\s+/, $_;
  
  if($decom[3] > 0)  {
    if($decom[1] > 0)  {
      $blast_par_sen += $decom[2] / $decom[1];
      $blast_par_spe += $decom[2] / $decom[3];
    } else  {
      $blast_par_sen += 1;
    }
  }  else  {
    if($decom[1] > 0)  {
      $blast_par_spe += 1;
    } else  {
      $blast_par_sen += 1;
      $blast_par_spe += 1;
    }
  }
  
  if($decom[5] > 0)  {
    if($decom[1] > 0)  {
      $psiblast_par_sen += $decom[4] / $decom[1];
      $psiblast_par_spe += $decom[4] / $decom[5];
    } else  {
      $psiblast_par_sen += 1;
    }
  }  else  {
    if($decom[1] > 0)  {
      $psiblast_par_spe += 1;
    } else  {
      $psiblast_par_sen += 1;
      $psiblast_par_spe += 1;
    }
  }

  if($decom[7] > 0)  {
    if($decom[1] > 0)  {
      $fastm_sen += $decom[6] / $decom[1];
      $fastm_spe += $decom[6] / $decom[7];
    } else  {
      $fastm_sen += 1;
    }
  }  else  {
    if($decom[1] > 0)  {
      $fastm_spe += 1;
    } else  {
      $fastm_sen += 1;
      $fastm_spe += 1;
    }
  }

  if($decom[9] > 0)  {
    if($decom[1] > 0)  {
      $grasp_sen += $decom[8] / $decom[1];
      $grasp_spe += $decom[8] / $decom[9];
    } else  {
      $grasp_sen += 1;
    }
  }  else  {
    if($decom[1] > 0)  {
      $grasp_spe += 1;
    } else  {
      $grasp_sen += 1;
      $grasp_spe += 1;
    }
  }
  ++ $num_count;
}
close $IN;

$blast_par_sen /= $num_count;
$blast_par_spe /= $num_count;
$psiblast_par_sen /= $num_count;
$psiblast_par_spe /= $num_count;
$fastm_sen /= $num_count;
$fastm_spe /= $num_count;
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

print "$blast_par_sen  $blast_par_spe  $psiblast_par_sen  $psiblast_par_spe  $fastm_sen  $fastm_spe  $grasp_sen  $grasp_spe\n";
