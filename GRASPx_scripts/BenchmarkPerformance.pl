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
#    system "perl /home/cczhong/Works/GuidedAssemble/Scripts/RandomSlice.pl $_ temp 32 20";
#    system "cat temp >>$run_folder/sim.fa";
  }
}

system "date";

# build blast database for sim.fa and full.fa
print "Building BLAST database...\n";
#system "makeblastdb -in $run_folder/sim.fa -dbtype prot";
#system "makeblastdb -in $run_folder/full.fa -dbtype prot";

my $large_num = 1000000000;

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
#  system "blastp -query $_ -db $run_folder/full.fa -num_descriptions 1000000 -num_alignments 1000000 -outfmt 7 -num_threads 2 >$run_folder/Coreset/$f_name.N$counter.blast.summary";
  ++ $counter;
}

system "date";

# construct the reads that are actually considered as the true hits
if(!(-e "$run_folder/Reference"))  {
#  system "mkdir $run_folder/Reference";
}

foreach(<$run_folder/Coreset/*.blast.summary>) {
  my $file = $_;
  $file =~ /.*\/(.*)\.blast\.summary/;
  my $f_name = $1;
  print "Extracting reference hits for $f_name ...\n";  
  #system "perl /usr/local/depot/projects/SPA/czhong/Works/GRASPx/Scripts/ListIncludedReads.pl $file $run_folder/Coreset $run_folder/sim.fa >$run_folder/Reference_BLASTP/$f_name.inc.reads";
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
  print "Runing BLASTP on $f_name against fragmented data set...\n";
#  system "blastp -query $_ -db $run_folder/sim.fa -evalue $large_num -num_descriptions 1000000 -num_alignments 1000000 -outfmt 7 -num_threads 2 >$run_folder/Short/$f_name.N$counter.blast.summary";
  ++ $counter;
}

system "date";

# running blast on the fragmented reads
if(!(-e "$run_folder/Psi"))  {
#  system "mkdir $run_folder/Psi";
}

$counter = 0;
foreach(@files_to_search) {
  $_ =~ /.*\/(.*)\.fa/;
  my $f_name = $1;
  print "Runing PSI-BLAST on $f_name against fragmented data set...\n";
  #system "/usr/local/packages/ncbi-blast+-2.2.28/bin/psiblast -query $_ -db $run_folder/sim.fa -evalue $large_num -num_alignments 1000000 -outfmt 7 -num_threads 2 -num_iterations 3 >$run_folder/Psi/$f_name.psibst";
  ++ $counter;
}

system "date";
#exit;

# running blast on the fragmented reads
if(!(-e "$run_folder/Fastm"))  {
#  system "mkdir $run_folder/Fastm";
}

$counter = 0;
foreach(@files_to_search) {
  $_ =~ /.*\/(.*)/;
  my $f_name = $1;
  print "Runing Fastm on $f_name against fragmented data set...\n";
  #system "/home/czhong/Tools/fasta-36.3.6/bin/fastm36 $_ $run_folder/sim.fa -s BL62 -C 100 -d $large_num -E 1000000 >$run_folder/Fastm/$f_name.N$counter.fastm.summary &";
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
#system "/home/cczhong/Codes/Grasp/Build --work_space=$run_folder/WorkSpace $run_folder/sim.fa";

#exit;

$counter = 0;
foreach(@files_to_search) {
  $_ =~ /.*\/(.*)/;
  my $f_name = $1;
  print "Runing GRASP on $f_name against fragmented data set...\n";
#  system "/home/cczhong/Codes/Grasp/Search $_ $run_folder/sim.fa --e_value=$large_num --work_space=$run_folder/WorkSpace --result_space=$run_folder/Grasp";
  ++ $counter;
}

system "date";
# building GRASPx index
print "Building GRASPx index...\n";
#system "/usr/local/depot/projects/SPA/czhong/GRASPx_current/grasp-build $run_folder/sim.fa --work_space=$run_folder/WorkSpace";

#exit;

$counter = 0;
foreach(@files_to_search) {
  $_ =~ /.*\/(.*)/;
  my $f_name = $1;
  print "Runing GRASPx on $f_name against fragmented data set...\n";
  #system "/usr/local/depot/projects/SPA/czhong/GRASPx_current/grasp-assemble $_ $run_folder/sim.fa $run_folder/GRASPx/$f_name.ctg --e_value=$large_num --work_space=$run_folder/WorkSpace\n";
  #system "/usr/local/depot/projects/SPA/czhong/GRASPx_current/grasp-map $run_folder/sim.fa $run_folder/GRASPx/$f_name.ctg $run_folder/GRASPx/$f_name.summary --work_space=$run_folder/WorkSpace --num_errors 0 --portion_mapped 0.99\n";
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
  #system "perl /usr/local/depot/projects/SPA/czhong/Works/GRASPx/Scripts/RawPerformancesCount_marineX.pl $run_folder 1e$exponent >$run_folder/Benchmark/performance.exp$exponent.count &";
}

for($exponent = -10; $exponent <= 1; ++ $exponent) {
  system "perl /usr/local/depot/projects/SPA/czhong/Works/GRASPx/Scripts/ComputePerformanceFromCount.pl $run_folder/Benchmark/performance.exp$exponent.count";
}

system "date";

