#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';

my $ID = shift;		# a list of Pfam IDs to simulate, one per line
my $pfam_info = shift;	# tab-delimited file, with columns meaning ID/Length/Name
my $out_file = shift;	# writing the output of simulated short reads

die "Output file already exists!\n" if -e "$out_file";

# load length information, skip profiles that are really short
my %long_profiles;
open my $IN, "<$pfam_info" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  $long_profiles{$decom[0]} = 1 if ($decom[1] >= 60);
}
close $IN;

my $pfam_full_dir = "/usr/local/projdata/0599/projects/SPA/czhong/Data/PfamFullSeq";

# get sequence one-by-one
open my $PIN, "<$ID" or die "Cannot open file: $!\n";
while(<$PIN>)  {
  chomp;
  my $pid = $_;
  next if !exists $long_profiles{$pid};
  my %file; 
  my $folder = $pfam_full_dir . '/' . $pid;
  #print "$folder\n";
  system "ls $folder/*.fa >tmp.ls";
  open my $LIN, "<tmp.ls" or die "Cannot open file: $!\n";
  while(<$LIN>)  {
    chomp;
    $file{$_} = 1;
  }
  close $LIN;
  my $n = 100;
  system "rm tmp.full" if -e "tmp.full";
  foreach(keys %file)  {
    system "cat $_ >>tmp.full";
    -- $n;
    last if $n <= 0;
  }
  # randomly generate reads 
  system "perl ~/ScriptRepo/FragmentSeqDB.pl --fasta=tmp.full --out=tmp.frag --length=32 --depth=10 --error=1";
  # incorporate the reads
  system "cat tmp.frag >>$out_file";
}
close $PIN;
