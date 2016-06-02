#!/usr/bin/perl -w
use strict;

my $run_folder = shift;
my $score_cutoff = shift;

my $inc_reads_folder = $run_folder . '/Reference/';
my $grasp_folder = $run_folder . '/Grasp_remap/';
my $blast_folder = $run_folder . '/Short/';
my $psiblast_folder = $run_folder . '/Psi/';
my $fastm_folder = $run_folder . '/Fastm/';

#my $partial_cons_frac = 0.5;

foreach(<$inc_reads_folder/*.ref>) {
  my $inc_read_file = $_;
  $inc_read_file =~ /(.*)\/(.*)\.ref/;
  my $file_stem = $2;
  #$file_stem =~ /(.*)\.fa/;
  my $seq_ID = $file_stem;
  my $query_header;
  #print "$file_stem\n";
  #print "$seq_ID.fa.0.reads	$file_stem.blast.summary\n";
  #next;
  if(-e "$grasp_folder/$seq_ID.fa.reads" and -e "$blast_folder/$file_stem.frag.bst"
    and -e "$psiblast_folder/$file_stem.frag.psibst" and -e "$fastm_folder/$file_stem.frag.fastm")  {
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
      #print "$full_header      $inc_hash{$full_header}\n";
      ++ $total_inc;
    }
    close $IIN;
    # handle the blast reads
    my %blast_partial_hash;
    my $total_blast_partial = 0;
    my $blast_partial_overlap_with_inc = 0;
    open my $BIN, "<$blast_folder/$file_stem.frag.bst" or die "Cannot open blast summary file: $!\n";
    #print "$file_stem\n";
    while(<$BIN>) {
      chomp;
      if(/^#/)  {
        next;
      }
      my @decom = split /\s+/, $_;
      #print "-----	$decom[1]	$seq_ID\n";
      if($decom[1] =~ /$seq_ID/ || $decom[11] > $score_cutoff)  {
        next;
      }
      #$decom[1] =~ /(.*)\:\:(\d+)-(\d+)/;
      #print "header blast:  $decom[1]\n";
      if(!(exists $blast_partial_hash{$decom[1]}))  {
        #print "new read found\n";
        if(exists $inc_hash{$decom[1]})  {
          ++ $blast_partial_overlap_with_inc;
        }
        ++ $total_blast_partial;
        $blast_partial_hash{$decom[1]} = 1;
      } 
    }
    close $BIN;
    #die;
    # handle the psi-blast reads
    my %psiblast_partial_hash;
    my $total_psiblast_partial = 0;
    my $psiblast_partial_overlap_with_inc = 0;
    open my $PIN, "<$psiblast_folder/$file_stem.psibst" or die "Cannot open blast summary file: $!\n";
    while(<$PIN>) {
      chomp;
      if(/^#/)  {
        next;
      }
      if(!defined($_))  {
        print "empty string\n";
      }
      my @decom = split /\s+/, $_;
      if(scalar(@decom) < 11 || $decom[1] =~ /$seq_ID/ || $decom[10] > $score_cutoff)  {
        next;
      }
      #$decom[1] =~ /(.*)\:\:(\d+)-(\d+)/;
      #print "header blast:  $decom[1]\n";
      if(!(exists $psiblast_partial_hash{$decom[1]}))  {
        if(exists $inc_hash{$decom[1]})  {
          ++ $psiblast_partial_overlap_with_inc;
        }
        ++ $total_psiblast_partial;
        $psiblast_partial_hash{$decom[1]} = 1;
      }
    }
    close $PIN;
    # handle the grasp reads
    my $total_grasp = 0;
    my $overlap_with_inc = 0;
    open my $GIN, "<$grasp_folder/$seq_ID.fa.reads" or die "Cannot open grasp summary file: $!\n";
    while(<$GIN>) {
      chomp;
      if(/^\#/)  {
        next;
      }
      my @decom = split /\s+/, $_;
      #print "header grasp:  $decom[0]\n";
      if(exists $map_ID2NAME{$decom[1]} and $map_ID2NAME{$decom[1]} =~ /$seq_ID/)  {
        next;
      }
      if($decom[2] < $score_cutoff)  {
        if(exists $map_ID2NAME{$decom[1]})  {
          ++ $overlap_with_inc;
        }
        ++ $total_grasp;
      }
    }
    # handle the fastm reads
    my $total_fastm = 0;
    my $fastm_overlap_with_inc = 0;
    my %fastm_hash;
    open my $FIN, "<$fastm_folder/$file_stem.frag.fastm" or die "Cannot open grasp summary file: $!\n";
    while(<$FIN>) {
      chomp;
      if(/^\>\>(\S*)\s+\(/)  {
        my $id = $1;
        my $info = <$FIN>;
        chomp $info;
        my @decom = split /\s+/, $info;
        #print "@@@ $id--	$decom[12]\n";
        if($id =~ /$seq_ID/ || $decom[12] > $score_cutoff)  {
          next;
        }
        #foreach(keys %inc_hash)  {
        #  print "!!!	$_	$inc_hash{$_}\n";
        #}
        
        if(!(exists $fastm_hash{$id}))  {
          #print "first seen	$id	$inc_hash{$id}\n";
          if(exists $map_ID2NAME{$id})  {
            #print "good hit\n";
            ++ $fastm_overlap_with_inc;
          }
          ++ $total_fastm;
          $fastm_hash{$id} = 1;
        }
        #die;
      }
    }
    close $FIN;
    # print out information
    print "$file_stem $total_inc $blast_partial_overlap_with_inc $total_blast_partial $psiblast_partial_overlap_with_inc $total_psiblast_partial $fastm_overlap_with_inc $total_fastm $overlap_with_inc $total_grasp\n";
  }
  #die;
}

