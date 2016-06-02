#!/usr/bin/perl -w
use strict;

# the script is used to format taxanomic mapping results into circos histogram format

my $file = shift;
my $out = shift;

my $taxa_mapping = "/usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Results/circos.karyotype.txt";

if(!defined $file || !defined $out)  {
  print "Usage: perl taxa_to_circos_hist.pl [TAXA_FILE] [OUTPUT]\n";
  print "Format of TAXA_FILE:\n";
  print "Lactobacillus casei     NODE_22_length_549_cov_14.0567_ID_43    NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      gi|823693350|ref|WP_047106176.1|        39.18   171     99      2       516     7       76      242     5e-25   109     109\nArcobacter butzleri     NODE_121_length_252_cov_1.25381_ID_241  NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      gi|822616330|ref|WP_046998362.1|        67.47   83      27      0       4       252     66      148     4e-33   127     127\nStreptococcus   NODE_131_length_243_cov_137.138_ID_261  NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      gi|568067975|ref|WP_024054795.1|        76.54   81      19      0       243     1       59      139     2e-33   128     128\nStreptococcus salivarius        NODE_41_length_432_cov_2.79045_ID_81    NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      gi|490286789|ref|WP_004182462.1|        100.00  144     0       0       432     1       81      224     7e-96   291     291\n";
  exit;
}

print "NOTE: using taxa-band mapping file $taxa_mapping\n";

my %taxa_map_hash;
my %taxa_count_hash;

# define presented taxa ID and mapping to bands
open my $TIN, "<$taxa_mapping" or die "Cannot open file: $!\n";
while(<$TIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  $taxa_map_hash{$decom[3]} = $decom[2];
  $taxa_count_hash{$decom[3]} = 0;
}
close $TIN;

# format the inputs
open my $OUT, ">$out" or die "Cannot create file: $!\n";
open my $IN, "<$file" or die "Cannot open file : $!\n";
while(<$IN>)  {
  chomp;
  my @decom = split /\t/, $_;
  $decom[0] =~ s/\s+/\_/g;
  my $begin = $taxa_count_hash{$decom[0]} * 10000;
  my $end = $begin + 10000 - 1;
  my $height = -log($decom[-2] + 1e-100)/log(10);
  print $OUT "$taxa_map_hash{$decom[0]}	$begin	$end	$height\n";
  ++ $taxa_count_hash{$decom[0]};
}
close $IN;
close $OUT;

