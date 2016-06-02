#!/usr/bin/perl -w
use strict;

my $inc_reads_folder = shift;
my $grasp_folder = shift;
my $blast_folder = shift;
my $score_cutoff = shift;

my $partial_cons_frac = 0.5;

foreach(<$inc_reads_folder/*.ref>) {
  my $inc_read_file = $_;
  $inc_read_file =~ /(.*)\/(.*)\.ref/;
  my $file_stem = $2;
  my $seq_ID = $file_stem;
  my $query_header;
  if(-e "$grasp_folder/$seq_ID.fa.0.reads" and -e "$blast_folder/$file_stem.frag.bst")  {
    # handle the included reads
    my $total_inc = 0;
    my %inc_hash;
    my %map_ID2NAME;
    my %map_NAME2ID;
    open my $IIN, "<$inc_read_file" or die "Cannot open included reads file: $!\n";
    while(<$IIN>) {
      chomp;
      my @decom = split /\s+/, $_;
      my $full_header = $decom[0];
      $full_header =~ /.*\/\d+\_(\d+)\_(\d+)/;
      my $read_len = ($2 - $1 + 1) / 3;
      $inc_hash{$full_header} = $read_len;
      $map_ID2NAME{$decom[1]} = $decom[0];
      $map_NAME2ID{$decom[0]} = $decom[1];
      #print "$full_header	$inc_hash{$full_header}\n";
      ++ $total_inc;
    }
    close $IIN;
    # handle the blast reads
    my %blast_full_hash;
    my %blast_partial_hash;
    my %blast_partial_cons_hash;
    my $total_blast_full = 0;
    my $total_blast_partial = 0;
    my $total_blast_partial_cons = 0;
    my $blast_full_overlap_with_inc = 0;
    my $blast_partial_overlap_with_inc = 0;
    my $blast_partial_cons_overlap_with_inc = 0;
    open my $BIN, "<$blast_folder/$file_stem.frag.bst" or die "Cannot open blast summary file: $!\n";
    while(<$BIN>) {
      chomp;
      if(/^#/)  {
        next;
      }
      my @decom = split /\s+/, $_;
      #print "-----	$decom[1]	$seq_ID\n";
      if($decom[11] > $score_cutoff)  {
        next;
      }
      $decom[1] =~ /.*\/\d+\_(\d+)\_(\d+)/;
      my $full_len = ($2 - $1 + 1) / 3;
      
      if(abs(($decom[10] - $decom[9] + 1) - $full_len) < 2)  {
        if(!(exists $blast_full_hash{$decom[1]}))  {
          if(exists $inc_hash{$decom[1]})  {
            ++ $blast_full_overlap_with_inc;
            #print "$decom[1]\n";
          }
          ++ $total_blast_full;
          $blast_full_hash{$decom[1]} = 1;
        }
      }
      if(($decom[7] - $decom[6] + 1) >= $full_len * $partial_cons_frac)  {
        if(!(exists $blast_partial_cons_hash{$decom[1]}))  {
          if(exists $inc_hash{$decom[1]})  {
            ++ $blast_partial_cons_overlap_with_inc;
          }
          ++ $total_blast_partial_cons;
          $blast_partial_cons_hash{$decom[1]} = 1;
        }
      } 
      #print "header blast:  $decom[1]\n";
      if(!(exists $blast_partial_hash{$decom[1]}))  {
        if(exists $inc_hash{$decom[1]})  {
          ++ $blast_partial_overlap_with_inc;
        }
        ++ $total_blast_partial;
        $blast_partial_hash{$decom[1]} = 1;
      } 
    }
    close $BIN;
    # handle the grasp reads
    my $total_grasp_max = 0;
    my $overlap_with_inc_max = 0;
    my $overlap_with_blast_full_max = 0;
    my $overlap_with_blast_partial_max = 0;
    my $total_grasp_min = 0;
    my $overlap_with_inc_min = 0;
    my $overlap_with_blast_full_min = 0;
    my $overlap_with_blast_partial_min = 0;
    open my $GIN, "<$grasp_folder/$seq_ID.fa.0.reads" or die "Cannot open grasp summary file: $!\n";
    while(<$GIN>) {
      chomp;
      if(/^\#/)  {
        next;
      }
      my @decom = split /\s+/, $_;
      #print "header grasp:  $decom[0]\n";
      if($decom[2] < $score_cutoff)  {
        if(exists $map_ID2NAME{$decom[1]})  {
          ++ $overlap_with_inc_max;
        }
        if(exists $map_ID2NAME{$decom[1]} and exists $blast_full_hash{$map_ID2NAME{$decom[1]}})  {
          ++ $overlap_with_blast_full_max;
        }
        if(exists $map_ID2NAME{$decom[1]} and exists $blast_partial_hash{$map_ID2NAME{$decom[1]}})  {
          ++ $overlap_with_blast_partial_max;
        }
        ++ $total_grasp_max;
      }
      #if($decom[2] < $score_cutoff)  {
      #  if(exists $inc_hash{$decom[0]})  {
      #    ++ $overlap_with_inc_min;
      #  }
      #  if(exists $blast_full_hash{$decom[0]})  {
      #    ++ $overlap_with_blast_full_min;
      #  }
      #  if(exists $blast_partial_hash{$decom[0]})  {
      #    ++ $overlap_with_blast_partial_min;
      #  }
      #  ++ $total_grasp_min;
      #}
      
      
    }
    close $GIN;
    # print out information

    print "$file_stem $total_inc $blast_full_overlap_with_inc $total_blast_full $blast_partial_cons_overlap_with_inc  $total_blast_partial_cons $blast_partial_overlap_with_inc $total_blast_partial $overlap_with_inc_max $total_grasp_max\n";
    #my $blast_full_sen = $blast_full_overlap_with_inc / $total_inc;
    #my $blast_full_spe = $blast_full_overlap_with_inc / $total_blast_full;
    #my $blast_partial_sen = $blast_partial_overlap_with_inc / $total_inc;
    #my $blast_partial_spe = $blast_partial_overlap_with_inc / $total_blast_partial;
    #my $grasp_sen = $overlap_with_inc / $total_inc;
    #my $grasp_spe = $overlap_with_inc / $total_grasp;
    #print "$blast_full_sen\t$blast_full_spe\t$blast_partial_sen\t$blast_partial_spe\t$grasp_sen\t$grasp_spe\n";
    #die;
  }
}

