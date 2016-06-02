#!/usr/bin/perl -w
use strict;

my $tab_file = shift;	# expecting merged RPKM or count file that used to generate the heatmap

my $function_map = "/usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Data/antiSMASH_metadata.txt";
my $taxa_map = "/usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Results/SPAdes_contigs.cdhit.BLAST.NT.taxa";

# load in function map
my %f_hash;
open my $FIN, "<$function_map" or die "Cannot open file: $!\n";
while(<$FIN>)  {
  chomp;
  my @decom = split /\t/, $_;
  $f_hash{$decom[2]} = $decom[0];
}
close $FIN;

# load in taxa map
my %t_hash;
open my $TIN, "<$taxa_map" or die "Cannot open file: $!\n";
while(<$TIN>)  {
  chomp;
  my @decom = split /\t/, $_;
  $t_hash{$decom[1]} = $decom[0];
}
close $TIN;

# load the table file and generate labels
open my $IN, "<$tab_file" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  if(/^\#/)  {
    print "$_\n";
    next;
  }
  my @decom = split /\t/, $_;
  my @decom2 = split /\|\|/, $decom[0];
  $decom2[1] =~ /(ID\_\d+)/;
  my $name = $1;
  my $function = $decom2[0];
  my $categ = "unknown_categ";
  $categ = $f_hash{$function} if (exists $f_hash{$function});
  my $taxa = "unknown_taxa";
  $taxa = $t_hash{$decom2[1]} if (exists $t_hash{$decom2[1]});
  my $label = $name . "\|\|" . $function . "\|\|" . $categ . "\|\|" . $taxa;
  $label =~ s/ /\_/g;
  print "$label	$decom[1]	$decom[2]	$decom[3]	$decom[4]	$decom[5]	$decom[6]	$decom[7]	$decom[8]\n";
}
close $IN;
