#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';

my $run_folder = shift;

$run_folder = abs_path($run_folder);
# get the project name
$run_folder =~ /.*\/(.*)/;
my $project_name = $1;

system "date";

# test if the folder contains "sim.fa" and "full.fa" files, if yes, delete them
if(-e "$run_folder/sim.fa")  {
  print "Simulated reads already exist... deleting the file...\n";
#  system "rm $run_folder/sim.fa";
}
if(-e "$run_folder/full.fa")  {
  print "Full length sequences already exist... deleting the file...\n";
#  system "rm $run_folder/full.fa";
}

# get the reads from the folder "Dataset"
my @files_to_search;
for(<$run_folder/Dataset/*>) {
  my $pfam_folder = $_;
  for(<$pfam_folder/*>) {
    push @files_to_search, $_;
    # record the full sequence
#    system "cat $_ >>$run_folder/full.fa";
    # simulate the short sequences
#    system "perl /usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Scripts/RandomSlice.pl $_ temp 32 20";
#    system "cat temp >>$run_folder/sim.fa";
  }
}

system "date";

# build blast database for sim.fa and full.fa
print "Building BLAST database...\n";
#system "/usr/local/packages/ncbi-blast+-2.2.28/bin/makeblastdb -in $run_folder/sim.fa -dbtype prot";
#system "/usr/local/packages/ncbi-blast+-2.2.28/bin/makeblastdb -in $run_folder/full.fa -dbtype prot";

my $large_num = 10000000000;

system "date";

# build core set as reference
if(!(-e "$run_folder/Coreset"))  {
#  system "mkdir $run_folder/Coreset";
}

my $counter = 0;
foreach(@files_to_search) {
  $_ =~ /.*\/(.*)/;
  my $f_name = $1;
  print "Runing BLAST on $f_name against full length data set...\n";
#  system "/usr/local/packages/ncbi-blast+-2.2.28/bin/blastp -query $_ -db $run_folder/full.fa -max_target_seqs 1000000 -outfmt 7 -num_threads 8 >$run_folder/Coreset/$f_name.N$counter.blast.summary";
  ++ $counter;
}

#exit;

system "date";

# construct the reads that are actually considered as the true hits
if(!(-e "$run_folder/Reference"))  {
#  system "mkdir $run_folder/Reference";
}

foreach(<$run_folder/Coreset/*.bst>) {
  my $file = $_;
  $file =~ /.*\/(.*)\.blast\.summary/;
  my $f_name = $1;
  print "Extracting reference hits for $f_name ...\n";
#  system "perl /usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Scripts/ListIncludedReads_new.pl $file $run_folder/Coverset $run_folder/sim.fa >$run_folder/Reference/$f_name.inc.reads";
#   system "perl /usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Scripts/ListIncludedReads.pl  >$run_folder/Reference_LenCutoff/$f_name.ref";
}

system "date";

#exit;

# running blast on the fragmented reads
if(!(-e "$run_folder/Short"))  {
#  system "mkdir $run_folder/Short";
}

$counter = 0;
foreach(@files_to_search) {
  $_ =~ /.*\/(.*)/;
  my $f_name = $1;
  print "Runing BLAST on $f_name against fragmented data set...\n";
#  system "/usr/local/packages/ncbi-blast+-2.2.28/bin/blastp -query $_ -db $run_folder/sim.fa -evalue $large_num -max_target_seqs 1000000 -outfmt 7 -num_threads 8 >$run_folder/Short/$f_name.N$counter.blast.summary";
  ++ $counter;
}

system "date";

# running grasp on the fragmented reads
if(!(-e "$run_folder/Grasp"))  {
#  system "mkdir $run_folder/Grasp";
}
if(!(-e "$run_folder/WorkSpace"))  {
#  system "mkdir $run_folder/WorkSpace";
}

# building index on the fragmented reads

print "Building Grasp indexing files...\n";
#system "/usr/local/depot/projects/SPA/czhong/Grasp/Build --work_space=$run_folder/WorkSpace $run_folder/sim.fa";

#exit;

$counter = 0;
foreach(@files_to_search) {
  $_ =~ /.*\/(.*)/;
  my $f_name = $1;
  print "Runing GRASP on $f_name against fragmented data set...\n";
  if(-e "$run_folder/Grasp/$f_name.0.reads")	{
    print "skipping...\n";
    next;
  }
#  system "/usr/local/depot/projects/SPA/czhong/Grasp/Search $_ $run_folder/sim.fa --e_value=$large_num --work_space=$run_folder/WorkSpace --result_space=$run_folder/Grasp";
  ++ $counter;
}

system "date";

#exit;

# running the benchmark
if(!(-e "$run_folder/Benchmark"))  {
#  system "mkdir $run_folder/Benchmark";
}
my $exponent;
for($exponent = -10; $exponent <= 1; ++ $exponent) {
  print "Trying e-value cutoff 1e$exponent...\n";
  #system "perl /usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Scripts/RawPerformancesCount.pl $run_folder 1e$exponent >$run_folder/Benchmark_tblastn_remap/performance.exp$exponent.count"
  system "perl /usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Scripts/BenchmarkWithBlast.pl $run_folder/Reference $run_folder/Grasp $run_folder/Short 1e$exponent >$run_folder/Benchmark_tblastn_local_noremap/performance.exp$exponent.count"
}

for($exponent = -10; $exponent <= 1; ++ $exponent) {
  system "perl /usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Scripts/ComputePerformanceFromCount.pl $run_folder/Benchmark_tblastn_local_noremap/performance.exp$exponent.count";
}

system "date";

