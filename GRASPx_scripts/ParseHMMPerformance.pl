#!/usr/bin/perl -w
use strict;

# example: perl ParseHMMPerformance.pl $MY_SPA_DIR/Workds/GuidedAssemble/Results/GRASP/derived/SGO $MY_SPA_DIR/Workds/GuidedAssemble/Results/GRASP/hmm_results/SGO
my $prefix_asm = shift;  # the prefix for the assembled sequences, expected to end with ".asm.faa"
my $prefix_hmm = shift;  # the prefix for the hmm search results, expected to end with ".hmm.search"

my $hmm_ecutoff = 1e-3;

my %grasp_evalue_hash;
my %hmm_evalue_hash;
my %read_recruitment_hash;

foreach(<$prefix_asm*.asm.faa>)  {
  # parse information regarding GRASP e-values and recruited reads
  my $fa_file = $_;
  $fa_file =~ /(.*)\.asm\.faa/;
  my $gid = $1;
  #print "$fa_file  $gid\n";
  open my $FIN, "<$fa_file" or die "Cannot open file: $!\n";
  while(<$FIN>)  {
    chomp;
    if(/^>(.*)/)  {
      my @decom = split /\s+/, $1;
      my $key = $gid . '-' . $decom[0];
      if(!exists $grasp_evalue_hash{$key} || $grasp_evalue_hash{$key} > $decom[1])  {
        $grasp_evalue_hash{$key} = $decom[1];
        $read_recruitment_hash{$key} = $decom[2];
      }
    }
  }
  close $FIN;
  # check for HMM alignment scores
  foreach(<$prefix_hmm*.hmm.search>)  {
    #print "$_\n";
    open my $HIN, "<$_" or die "Cannot open file: $!\n";
    my $record_tag = 0;
    while(<$HIN>)  {
      chomp;
      my $line = $_;
      $line =~ s/^\s+//g;
      if(!$record_tag && $line =~ /^\d+/)  {
        $record_tag = 1;
      }
      if($record_tag && $line =~ /^\d+/)  {
        my @decom = split /\s+/, $line;
        my $key = $gid . '-' . $decom[8];
        if(!exists $hmm_evalue_hash{$key} || $hmm_evalue_hash{$key} > $decom[0])  {
          $hmm_evalue_hash{$key} = $decom[0];
        }
      }  elsif($record_tag && !($line =~ /^\d+/))  {
        last;
      }
    }
    close $HIN;
  }
}
  
my $i;
for($i = 1; $i >= -10; -- $i)  {
  my $cutoff = "1e" . $i;
  my %all_hash;
  my %good_hash;
  foreach(keys %grasp_evalue_hash)  {
    my $id = $_;
    if($grasp_evalue_hash{$id} > $cutoff)  {
      next;
    }
    my @decom;
    if(defined $read_recruitment_hash{$id})  {
      @decom = split /\;/, $read_recruitment_hash{$id};
    }
    foreach(@decom)  {
      if(/(.*)\,/)  {
        $all_hash{$1} = 1;
        if(exists $hmm_evalue_hash{$id} && $hmm_evalue_hash{$id} <= $hmm_ecutoff)  {
          $good_hash{$1} = 1;
        }
      }
    }
  }
  my $num_good = scalar keys %good_hash;
  my $num_all = scalar keys %all_hash;
  my $spe = $num_good / $num_all;
  print "$cutoff  $num_good  $spe\n";  
}

