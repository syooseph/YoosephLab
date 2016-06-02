#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';
use Getopt::Long;

my $model_list;		# a list of HMM models to be searched
my $read_set;		# the sequences to be searched from
my $ref_seq;		# the reference sequences (unfragmented sequences)
my $work_dir;		# the working directory

# a set of assumed fixed paths (the databases and the executables)
my $HMMGRASPx_dir = "/usr/local/projdata/0599/projects/SPA/czhong/HMMGRASPx-release";
my $Script_dir = "/usr/local/projdata/0599/projects/SPA/czhong/Works/GRASPxp/Scripts";
my $BLAST_dir = "/usr/local/packages/ncbi-blast+/bin";
my $PFAM_dir = "/usr/local/projdata/0599/projects/SPA/czhong/Data/Pfam_HMM";
my $CDD_dir = "/usr/local/projdata/0599/projects/SPA/czhong/Data/CDD";
my $PFAM_SeedSeq_dir = "/usr/local/projdata/0599/projects/SPA/czhong/Data/PfamSeedSeq";
my $PFAM_FullSeq_dir = "/usr/local/projdata/0599/projects/SPA/czhong/Data/PfamFullSeq";
my $UProC_dir = "/home/czhong/Tools/uproc-1.2.0";

GetOptions (
  "hmm_list=s" => \$model_list,
  "seq=s" => \$read_set,
  "ref=s" => \$ref_seq,
  "work_dir=s" => \$work_dir
) or die("Error in command line arguments\n");

if(!defined $model_list || !defined $read_set || !defined $work_dir || !defined $ref_seq)  {
  print "	--hmm_list:	the list of PFAM IDs of interest (), one perl line\n";
  print "	--seq:		the FASTA file for short peptides\n";
  print "	--ref:		the reference (unfragmented) sequences where the short peptides are sampled\n";
  print "	--work_dir:	the working directory\n";
  exit();
}

# create directories
$work_dir = abs_path($work_dir);
if(! -e "$work_dir")  {
  mkdir $work_dir or die "Cannot create directory: $!\n";
}


# dump all HMMs into one file
print "Preparing sequences...\n";
open my $LIN, "<$model_list" or die "Cannot open file: $!\n";
my %unique_id;
while(<$LIN>)  {
  chomp;
  /PF(\d+)/;
  my $pfam_id = $1;
  next if exists $unique_id{$pfam_id};
  $unique_id{$pfam_id} = 1;
  if(-e "$PFAM_dir/PF$1.hmm" and -e "$CDD_dir/pfam$1.smp" and -d "$PFAM_SeedSeq_dir/PF$1" and -d "$PFAM_FullSeq_dir/PF$1")  {
    # prepare hmmsearch and HMM-GRASPx inputs
    system "cat $PFAM_dir/PF$1.hmm >>$work_dir/model.hmm";
    # prepare RPS-BLAST inputs
    #system "echo $CDD_dir/pfam$1.smp >>$work_dir/smp.list";
    # prepare UProC seed inputs
    #system "cat $PFAM_SeedSeq_dir/PF$1/* >>$work_dir/PFseed.fa";
    # prepare UProC full inputs
    #system "cat $PFAM_FullSeq_dir/PF$1/* >>$work_dir/PFfull.fa";
  }  else  {
    die "HMM model $_ does not exists!\n";
  }
}
close $LIN;

# run HMM-GRASPx pipeline
print "Running HMM-GRASPx...\n";
system "/usr/bin/time -v perl $HMMGRASPx_dir/scripts/RunHMMGRASPx.pl --hmm=$work_dir/model.hmm --seq=$read_set --home=$HMMGRASPx_dir --param=$HMMGRASPx_dir/settings/param --index=$HMMGRASPx_dir/index --out=$HMMGRASPx_dir/out";

# run hmmsearch to define the Ture-positive set
print "Defining True-positive set using hmmsearch...\n";
system "$HMMGRASPx_dir/thirdparty/hmmer-3/binaries/hmmsearch $work_dir/model.hmm $ref_seq >$work_dir/tp.raw";
system "perl $Script_dir/GetHmmerIdentifiedSeq.pl $work_dir/tp.raw $ref_seq >$work_dir/tp.hit.fa";
system "$HMMGRASPx_dir/bin/graspx-map $read_set $work_dir/tp.hit.fa $work_dir/tp.map --work_space=$HMMGRASPx_dir/index --num_errors=3 --portion_mapped=0.6";
system "perl $Script_dir/map_to_TP.pl $work_dir/tp.map >$work_dir/tp.ref";

# run benchmark script
system "perl $Script_dir/compute_graspxp_performance.pl $work_dir/tp.ref $read_set $HMMGRASPx_dir/out/recruited.list >$work_dir/graspxp.summary";
system "cp $HMMGRASPx_dir/out/* $work_dir/HMMGRASPx_out";

# run hmmsearch
print "Running hmmsearch...\n";
#system "/usr/bin/time -v $HMMGRASPx_dir/thirdparty/hmmer-3/binaries/hmmsearch --cpu 16 $work_dir/model.hmm $read_set >$work_dir/hmmsearch.raw";
#system "perl $Script_dir/hmmer3_to_graspxp_map.pl $work_dir/hmmsearch.raw 10 >$work_dir/hmmsearch.map";
system "perl $Script_dir/compute_graspxp_performance.pl $work_dir/tp.ref $read_set $work_dir/hmmsearch.map >$work_dir/hmmsearch.summary";

# run UProC on Pfam seed alignments
print "Running UProC seed...\n";
#system "perl $Script_dir/reformat_header.pl $work_dir/PFseed.fa >$work_dir/PFseed.formatted.fa";
#system "/usr/bin/time -v $UProC_dir/uproc-makedb $UProC_dir/model $work_dir/PFseed.formatted.fa $work_dir/uproc_db_seed";
#system "/usr/bin/time -v $UProC_dir/uproc-prot -p $work_dir/uproc_db_seed $UProC_dir/model $read_set -o $work_dir/uproc.seed.raw";
#system "perl $Script_dir/uprot_to_graspx_map.pl $work_dir/model.hmm $work_dir/uproc.seed.raw >$work_dir/uproc.seed.map";
system "perl $Script_dir/compute_graspxp_performance.pl $work_dir/tp.ref $read_set $work_dir/uproc.seed.map >$work_dir/uproc.seed.summary";

# run UProC on Pfam full alignments
print "Running UProC full...\n";
#system "perl $Script_dir/reformat_header.pl $work_dir/PFfull.fa >$work_dir/PFfull.formatted.fa";
#system "/usr/bin/time -v $UProC_dir/uproc-makedb $UProC_dir/model $work_dir/PFfull.formatted.fa $work_dir/uproc_db_full";
#system "/usr/bin/time -v $UProC_dir/uproc-prot -t 16 -p $work_dir/uproc_db_full $UProC_dir/model $read_set -o $work_dir/uproc.full.raw";
#system "perl $Script_dir/uprot_to_graspx_map.pl $work_dir/model.hmm $work_dir/uproc.full.raw >$work_dir/uproc.full.map";
system "perl $Script_dir/compute_graspxp_performance.pl $work_dir/tp.ref $read_set $work_dir/uproc.full.map >$work_dir/uproc.full.summary";

# run RPS-BLAST
print "Running RPS_BLAST...\n";
#system "$BLAST_dir/makeprofiledb -in $work_dir/smp.list -out $work_dir/rpsdb";
#system "/usr/bin/time -v $BLAST_dir/rpsblast -query $read_set -db $work_dir/rpsdb -outfmt 7 -evalue 0.001 >$work_dir/rps.raw";
#system "perl $Script_dir/rpsblast_to_map.pl $work_dir/smp.list $work_dir/rps.raw 0.001 >$work_dir/rps.map";
system "perl $Script_dir/compute_graspxp_performance.pl $work_dir/tp.ref $read_set $work_dir/rps.map >$work_dir/rps.summary";
