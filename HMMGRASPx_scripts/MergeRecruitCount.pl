#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';

# tabulate the read counts from a list of recruit file (HMM-GRASPx format)
# pooled by genes

my $folder = shift;	# expected to process all files in the folder (keep only recruit files in)

my $current_dir = abs_path("./");
chdir $folder or die "Cannot change directory: $!\n";

#count the number of samples
my @ids;
foreach(<*>)  {
  push @ids, $_;
}

# prepare empty table
my @foo;
foreach(@ids)  {
  push @foo, 0;		# initialize the count of unseen genes to 0
}

# count samples one-by-one
my %gene_count;
my $n = 0;
foreach(<*>)  {
  open my $IN, "<$_" or die "Cannot open file: $!\n";
  while(<$IN>)  {
    chomp;
    my @decom = split /\s+/, $_;
    my @decom2 = split /\|\|/, $decom[2];
    #print "$decom2[1]	$n\n";
    if(!exists $gene_count{$decom2[1]})  {
      @{$gene_count{$decom2[1]}} = @foo;
    } 
    ++ $gene_count{$decom2[1]}->[$n];
  }
  close $IN;
  ++ $n;
}

chdir "$current_dir" or die "Cannot change directory: $!\n";

# print the results
print "@ids\n";
foreach(sort keys %gene_count)  {
  print "$_	@{$gene_count{$_}}\n";
}

