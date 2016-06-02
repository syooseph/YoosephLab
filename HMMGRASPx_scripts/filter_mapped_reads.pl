#!/usr/bin/perl -w
use strict;

my $peptide_file = shift;	# the complete peptide file output by FragGeneScan
my $bam_file = shift;		# the bam file that contains the mapped reads to reference
my $unmapped = shift;		# the resulting unmapped reads in multi-FASTA format

# first need to use Samtools to get the mapped reads from the bam file
print "Running samtools to convert BAM to SAM...\n";
system "/usr/local/packages/samtools/bin/samtools view $bam_file >filter_mapped.tmp.sam";
#system "cp $bam_file filter_mapped.tmp.sam";

# loads all mapped read IDs from the temporary sam file
print "Loading SAM file to find mapped reads...\n";
open my $SIN, "<filter_mapped.tmp.sam" or die "Cannot open file: $!\n";
my %mapped_reads;
while(<$SIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  #print "$decom[0]\n";
  $mapped_reads{$decom[0]} = 1;
}
close $SIN;

# reads in the peptide file and filter out unmapped reads
print "Filtering mapped reads...\n";
# unzip the peptide reads
system "zcat $peptide_file >filter_mapped.pep.fa";
#system "cp $peptide_file filter_mapped.pep.fa";
open my $PIN, "<filter_mapped.pep.fa" or die "Cannot open file: $!\n";
open my $OUT, ">$unmapped" or die "Cannot create file: $!\n";
while(<$PIN>)  {
  if(/^>(.*)/)  {
    my $header = $1;
    my @decom = split /\_/, $header;	# Note that this only works the Anna's data set, 
                                        # no guarantee for other file headers
    if(!(exists $mapped_reads{$decom[0]}))  {
      print $OUT ">$header\n";
      my $line = <$PIN>;
      print $OUT "$line";
    }
  }
}
close $PIN;
