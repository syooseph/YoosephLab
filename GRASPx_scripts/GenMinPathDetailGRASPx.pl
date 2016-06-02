#!/usr/bin/perl -w
use strict;

my $kg_path_file = shift;	# the KEGG pathway file that contains gene ID and the correponding pathway
my $kg_ko_file = shift;		# the KEGG KO file that contains the gene ID and KO ID mapping
my $graspx_map = shift;		# the GRASPx mapping file
my $cutoff = shift;		# the E-value cutoff for parsing GRASPx results

# load the kg_path_file and construct a hash table {gene ID} => {pathway ID}
sub ParseKEGGPathwayFile($)  {
  my $file = shift;
  my %g2p_hash;	# gene ID to pathway ID
  open my $IN, "<$file" or die "Cannot open file: $!\n";
  while(<$IN>)  {
    chomp;
    my @decom = split /\s+/, $_;
    $g2p_hash{$decom[0]} = $decom[1];
  }
  close $IN;
  return %g2p_hash;
}

# load the kg_ko_file and construct a hash table {gene ID} => {KO ID}
sub ParseKEGGKOFile($)  {
  my $file = shift;
  my %g2k_hash; # gene ID to pathway ID
  open my $IN, "<$file" or die "Cannot open file: $!\n";
  while(<$IN>)  {
    chomp;
    my @decom = split /\s+/, $_;
    $g2k_hash{$decom[0]} = $decom[1];
  }
  close $IN;
  return %g2k_hash;
}

# parse the GRASPx mapping file and count the number of mapping to each gene
sub CountGRASPxMapping($$)  {
  my $file = shift;
  my $cutoff = shift;
  my %map_count_hash;	# gene ID to mapping count
  open my $IN, "<$file" or die "Cannot open file: $!\n";
  while(<$IN>)  {
    chomp;
    # the following steps are GRASPx-format specific
    # an example is:
    # 3 HWI-EAS404_103089860:5:100:10189:10815/2_1_100_-  sgo:SGO_1404||contig_35::0-251::5.95483e-64::246.514::628
    my @decom = split /\s+/, $_;
    my @m_info = split /\:\:/, $decom[2]; # the mappin information
    if($m_info[0] =~ /(.*)\|\|contig/ && $m_info[2] <= $cutoff)  {
      $map_count_hash{$1} += 1;
    }
  }
  close $IN;
  return %map_count_hash;
}

# load in information from the files
my %g2p_hash = ParseKEGGPathwayFile($kg_path_file);
my %g2k_hash = ParseKEGGKOFile($kg_ko_file);
#my %map_count_hash = CountGRASPxMapping($graspx_map, $cutoff);

open my $IN, "<$graspx_map" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  # the following steps are GRASPx-format specific
  # an example is:
  # 3 HWI-EAS404_103089860:5:100:10189:10815/2_1_100_-  sgo:SGO_1404||contig_35::0-251::5.95483e-64::246.514::628
  my @decom = split /\s+/, $_;
  my @m_info = split /\:\:/, $decom[2]; # the mappin information
  if($m_info[0] =~ /(.*)\|\|contig/ && $m_info[2] <= $cutoff)  {
    if(exists $g2k_hash{$1})  {
      $g2k_hash{$1} =~ /(\d+)/;
      print "$decom[0]	K$1\n";
    }
    #$map_count_hash{$1} += 1;
  }
}
close $IN;

