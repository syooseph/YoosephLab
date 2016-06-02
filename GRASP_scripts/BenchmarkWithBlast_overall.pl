#!/usr/bin/perl -w
use strict;

my $inc_reads_folder = shift;
my $grasp_folder = shift;
my $blast_folder = shift;
my $bit_score_cutoff = shift;

my $total_TP = 0;
my $total_full_blast_TP = 0;
my $total_partial_blast_TP = 0;
my $total_grasp_TP = 0;
my $total_full_blast_output = 0;
my $total_partial_blast_output = 0;
my $total_grasp_output = 0;

foreach(<$inc_reads_folder/*.inc.reads>) {
  my $inc_read_file = $_;
  $inc_read_file =~ /(.*)\/(.*)\.inc\.reads/;
  my $file_stem = $2;
  $file_stem =~ /(.*)\.fa/;
  my $seq_ID = $1;
  if(-e "$grasp_folder/$file_stem.grasp.summary" and -e "$blast_folder/$file_stem.blast.summary")  {
    # handle the included reads
    my $total_inc = 0;
    my %inc_hash;
    open my $IIN, "<$inc_read_file" or die "Cannot open included reads file: $!\n";
    while(<$IIN>) {
      chomp;
      my @decom = split /\s+/, $_;
      if($decom[1] < $bit_score_cutoff)  {
        #print "$file_stem $decom[0] $decom[1] $bit_score_cutoff\n";
        #next;
      }
      my $full_header = $decom[0];
      $full_header =~ /(.*)\:\:(\d+)-(\d+)/;
      my $read_len = $3 - $2 + 1;
      $inc_hash{$full_header} = $read_len;
      ++ $total_inc;
    }
    close $IIN;
    # handle the blast reads
    my %blast_full_hash;
    my %blast_partial_hash;
    my $total_blast_full = 0;
    my $total_blast_partial = 0;
    my $blast_full_overlap_with_inc = 0;
    my $blast_partial_overlap_with_inc = 0;
    open my $BIN, "<$blast_folder/$file_stem.blast.summary" or die "Cannot open blast summary file: $!\n";
    while(<$BIN>) {
      chomp;
      if(/^#/)  {
        next;
      }
      my @decom = split /\s+/, $_;
      if($decom[1] =~ /$seq_ID/ || $decom[11] < $bit_score_cutoff)  {
        next;
      }
      $decom[1] =~ /(.*)\:\:(\d+)-(\d+)/;
      my $full_len = $3 - $2 + 1;
      if($decom[7] - $decom[6] + 1 == $full_len)  {
        $blast_full_hash{$decom[1]} = 1;  
        if(exists $inc_hash{$decom[1]})  {
          ++ $blast_full_overlap_with_inc;
        }
        ++ $total_blast_full;
      } 
      $blast_partial_hash{$decom[1]} = 1;
      if(exists $inc_hash{$decom[1]})  {
        ++ $blast_partial_overlap_with_inc;
      }
      ++ $total_blast_partial; 
    }
    close $BIN;
    # handle the grasp reads
    my $total_grasp = 0;
    my $overlap_with_inc = 0;
    my $overlap_with_blast_full = 0;
    my $overlap_with_blast_partial = 0;
    open my $GIN, "<$grasp_folder/$file_stem.grasp.summary" or die "Cannot open grasp summary file: $!\n";
    while(<$GIN>) {
      chomp;
      if(/^\W/)  {
        next;
      }
      my @decom = split /\s+/, $_;
      if($decom[0] =~ /$seq_ID/ || $decom[1] < $bit_score_cutoff)  {
        next;
      }
      if(exists $inc_hash{$decom[0]})  {
        ++ $overlap_with_inc;
      }
      if(exists $blast_full_hash{$decom[0]})  {
        ++ $overlap_with_blast_full;
      }
      if(exists $blast_partial_hash{$decom[0]})  {
        ++ $overlap_with_blast_partial;
      }
      ++ $total_grasp;
    }
    close $GIN;
    # print out information
    #print "$file_stem $total_inc $blast_full_overlap_with_inc $total_blast_full  $blast_partial_overlap_with_inc $total_blast_partial $overlap_with_inc $total_grasp $overlap_with_blast_full  $overlap_with_blast_partial\n";
    $total_TP += $total_inc;
    $total_full_blast_TP += $blast_full_overlap_with_inc;
    $total_partial_blast_TP += $blast_partial_overlap_with_inc;
    $total_full_blast_output += $total_blast_full;
    $total_partial_blast_output += $total_blast_partial;
    $total_grasp_TP += $overlap_with_inc;
    $total_grasp_output += $total_grasp;
  }
}

my $blast_full_sen = $total_full_blast_TP / $total_TP;
my $blast_full_spe = $total_full_blast_TP / $total_full_blast_output;
my $blast_partial_sen = $total_partial_blast_TP / $total_TP;
my $blast_partial_spe = $total_partial_blast_TP / $total_partial_blast_output;
my $grasp_sen = $total_grasp_TP / $total_TP;
my $grasp_spe = $total_grasp_TP / $total_grasp_output;

print "$blast_full_sen  $blast_full_spe $blast_partial_sen  $blast_partial_spe  $grasp_sen  $grasp_spe\n";
