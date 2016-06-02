#!/usr/bin/perl -w
use strict;

my $read_file = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/marine.core.faa";
my $annot_file = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/marine.core.amphora2.region.annot";
my $map_file = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/marine_sim_amphora2/read_remap.csv";
my $output_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/marine_sim_amphora2/Reference_Annot";

# load in the read file
my $num_reads = 0;
open my $RIN, "<$read_file" or die "Cannot open read file: $!\n";
my %id_hash;
my $index = 0;
while(<$RIN>) {
  chomp;
  if(/^>(.*)/)  {
    $id_hash{$1} = $index;
    ++ $index;
  }
}
close $RIN;

# load in the annotation file
my @annot_info;
open my $AIN, "<$annot_file" or die "Cannot open file: $!\n";
while(<$AIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  my @info = split /\-/, $decom[0];
  $decom[2] =~ /(\d+)\.\.(\d+)/;
  my $begin = $1; my $end = $2;
  my $recover_key = $info[0].'-'.$info[1].'-'.$info[2];
  if($decom[2] =~ /^complement/)  {
    push @annot_info, [$recover_key, $info[5], $end, $begin];
  }  else  {
    push @annot_info, [$recover_key, $info[5], $begin, $end];
  }
}
close $AIN;

# load in the mapping file info
my @map_info;
open my $MIN, "<$map_file" or die "Cannot open file: $!\n";
while(<$MIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  $decom[0] =~ /ref\|(NC_\d+)/;
  push @map_info, [$decom[0], $1, $decom[3], $decom[4]];
  #if($decom[3] > 515315 && $decom[4] < 515938)  {
  #  print "$_\n";
  #}
}
close $MIN;

#die;

my %read_hash;
# find pairwise correlation
my $i = 0;
my $j = 0;
foreach($i = 0; $i < scalar(@annot_info); ++ $i)  {
  print "Working on search file $annot_info[$i][0]...\n";
  for($j = 0; $j < scalar(@map_info); ++ $j)  {
    if($map_info[$j][1] eq $annot_info[$i][1])  {
      #if($map_info[$j][2] > 515315 && $map_info[$j][3] < 515938)  {
      #  print "****** $annot_info[$i][2]	$annot_info[$i][3]	$map_info[$j][2]	$map_info[$j][3]\n";
      #}
      if($map_info[$j][2] < $map_info[$j][3] && $annot_info[$i][2] < $annot_info[$i][3] &&
         !($map_info[$j][3] < $annot_info[$i][2] || $map_info[$j][2] > $annot_info[$i][3]))  {
        #print "read taken\n";
        push @{$read_hash{$annot_info[$i][0]}}, $map_info[$j][0];
      }  elsif($map_info[$j][2] > $map_info[$j][3] && $annot_info[$i][2] > $annot_info[$i][3] &&
         !($map_info[$j][2] < $annot_info[$i][3] || $map_info[$j][3] > $annot_info[$i][2]))  {
        #print "read taken\n";
        push @{$read_hash{$annot_info[$i][0]}}, $map_info[$j][0];
      }  else  {
        #print "read droppped\n";
      }
    }
  }
}

foreach(keys %read_hash)  {
  open my $OUT, ">$output_dir/$_.ref" or die "Cannot create file: $!\n";
  foreach(@{$read_hash{$_}})  {
    print $OUT "$_	$id_hash{$_}\n";
  }
  close $OUT;
}
