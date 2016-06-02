# This is a perl script for computing the RPKM of each protein family from the HMM-GRASPx results
#
# Author: Cuncong Zhong (for comments and questions please email czhong at jcvi.org)

#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;

my $models;		# query HMMs (in HMMER3 format) for the protein/protein domain families
my $seq_db;		# the multi-fasta file for the short-peptide database
my $file;		# the recruited reads
my $col = 0;		# the column ID for the gene name
my $out;		# the output directory
my $label = "";		# the label to be append to each header

GetOptions (
  "hmm=s" => \$models,
  "seq=s" => \$seq_db,
  "file=s" => \$file,
  "col=i" => \$col,
  "out=s" => \$out,
  "label=s" => \$label
) or die "Error in command line arguments!\n";

if(!defined $models || !defined $seq_db || !defined $file || !defined $out)  {
  print "CountGeneAbundance.pl: a perl script that estimate the abundance level of each gene\n";
  print "Usage: perl CountGeneAbundance.pl --hmm=[QUERY_HMM] --seq=[SHORT-PEPTIDE_READS] --recruit=[LOCATION_OF_recruited.list] --out=[OUTPUT_FILE]\n";
  print "	--hmm:		the query HMMs, with a format as defined in HMMER3\n";
  print "	--seq:		the short-peptide reads, in FASTA format\n";
  print "	--file:		the location of the \"recruit.list\" file, should be found from output folder\n";
  print "	--col:		the column ID for the gene name\n";
  print "	--out:		output file\n";
  print "	--label:	the label for the header\n";
  exit;
}

# count the number of reads
my $total_reads = 0;
open my $SEQIN, "<$seq_db" or die "Cannot open sequence file: $!\n";
while(<$SEQIN>)  {
  ++ $total_reads;
}
close $SEQIN; 
$total_reads = $total_reads / 2;

# get the lengths of the profile
my %length;
open my $MIN, "<$models" or die "Cannot open query HMM file: $!\n";
my $name;
while(<$MIN>)  {
  if(/NAME\s+(.*)/)  {
    $name = $1;
  }  elsif(/LENG\s+(\d+)/)  {
    $length{$name} = $1;
  }
}
close $MIN;

# count the number of homolog reads for each gene
my %count;
open my $CIN, "<$file" or die "Cannot open recruited-reads file: $!\n";
while(<$CIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  ++ $count{$decom[$col]};
}
close $CIN;

# compute "RPKM" and output
open my $OUT, ">$out" or die "Cannot create file for output: $!\n";
if($label ne "")  {
  print $OUT "\#$label\_NAME	$label\_RAW_COUNT	$label\_RPKM\n";
}  else  {
  print $OUT "\#NAME	RAW_COUNT	RPKM\n";
}
foreach(sort keys %length)  {
  my $c = 0;
  $c = $count{$_} if exists $count{$_};
  my $rpkm = 1000000000 * $c / $total_reads / $length{$_};
  print $OUT "$_	$c	$rpkm\n";
  delete $count{$_};
}
close $OUT;

# check if all mapped reads belong to one of the queries
if(scalar keys %count != 0)  {
  print "Potential error: some mapped reads do not belong to any of the query model. Did you specify the correct input???\n";
}
