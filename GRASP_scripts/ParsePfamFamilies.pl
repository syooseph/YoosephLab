#!/usr/bin/perl -w
use strict;

my $pfam_file = shift;  # PFAM alignment, in stockholme format
my $directory = shift;  # the output directory

if(!(-e "$directory"))  {
  mkdir "$directory" or die "Cannot create directory: $!\n";
}

open my $IN, "<$pfam_file" or die "Cannot open file: $!\n";

my @contents;
while(<$IN>) {
  chomp;
  push @contents, $_;
  if($_ eq '//')  {   # the end of a family
    if($contents[0] eq '# STOCKHOLM 1.0' && $contents[-1] eq '//')  {  # valid infomation
      my $i;
      my $family_name;    # name of the protein family
      my $family_ID;      # accession number of the protein family
      my %sequence_hash;  # hash map between the individual sequence ID and their sequences
      for($i = 0; $i < scalar(@contents); ++ $i) {
        my @decom = split /\s+/, $contents[$i];
        if($decom[0] eq '#=GF' and $decom[1] eq 'ID')  {    # recording family name
          $family_name = $decom[2];
        } elsif($decom[0] eq '#=GF' and $decom[1] eq 'AC') {    # recording accessing number
          $family_ID = $decom[2];
        } elsif($decom[0] eq '#=GC' and $decom[1] eq 'seq_cons') {
          $decom[2] =~ s/\.//g;
          $decom[2] =~ s/\W//g;
          $decom[2] = uc($decom[2]);
          $sequence_hash{"consensus"} .= $decom[2];
        } elsif(scalar(@decom) == 2 and !($decom[0] =~ /^\#/)) {    # recording sequence
          $decom[1] =~ s/\.//g;
          $sequence_hash{$decom[0]} .= $decom[1];
        }
      }
      # output the sequences
      if(!(-e "$directory/$family_ID"))  {
        mkdir "$directory/$family_ID" or die "Cannot create directory: $!\n";
      }
      foreach(keys %sequence_hash) {
        my $ori_seq_ID = $_;
        my $seq_ID = $_;
	      $seq_ID =~ s/\//\_/g;
        my $out_file_name = "$directory/$family_ID/$seq_ID.fa";
        open my $OUT, ">$out_file_name" or die "Cannot write file: $!\n";
        print $OUT ">$family_ID\:\:$family_name\:\:$seq_ID\n";
        print $OUT "$sequence_hash{$ori_seq_ID}\n";
        close $OUT;
      }
    }
    undef @contents;  # clear the current information block
  }
}

close $IN;
