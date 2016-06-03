# This is a perl driver for HMM-GRASPx, which automatically runs the entire analysis pipeline.
# The pipeline consists of the HMM-guided assembly performed by HMM-GRASPx,
# followed by a post-processing step that uses hmmsearch (HMMER3 package) to recompute the statistics.
#
# Author: Cuncong Zhong (for comments and questions please email czhong at jcvi.org)

#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;

my $models;		# query HMMs (in HMMER3 format) for the protein/protein domain families
my $seq_db;		# the multi-fasta file for the short-peptide database
my $out_dir = "out";	# the output directory
my $hm_dir = ".";	# the home directory for HMM-GRASPx
my $idx_dir = "index";	# the indexing directory
my $paramf;		# the sample sheet that document important parameters
my $unique = 0;		# allow only unique contigs
my $help = 0;		# print help information


GetOptions (
  "hmm=s" => \$models,
  "seq=s" => \$seq_db,
  "out=s" => \$out_dir,
  "home=s" => \$hm_dir,
  "index=s" => \$idx_dir,
  "param=s" => \$paramf,
  "unique" => \$unique,
  "help" => \$help
) or die("Error in command line arguments\n");

if($help || !defined $models || !defined $seq_db || !defined $hm_dir)  {
  print "RunHMMGRASPx.pl: a perl driver for the HMM-GRASPx pipeline\n";
  print "Usage:  perl RunHMMGRASPx.pl --hmm=[QUERY_HMM] --seq=[SHORT-PEPTIDE_READS] --out=[OUT_DIR]\n";
  print "	--hmm:		the query HMMs, with a format as defined in HMMER3\n";
  print "	--seq:		the short-peptide reads, in FASTA format\n";
  print "	--out:		the output directory\n";
  print "	--home:		the home directory of HMM-GRASPx (optional, default \".\/\")\n";
  print "	--index:	the index directory containing pre-built indexes (optional, default \".\/index\")\n";
  print "	--param:	the parameter setting file (optional, see examples/param_example)\n";
  print "	--unique:	only consider unique contig regions\n";
  print "			(assembled by using and best aligned to the same protein family)\n";
  exit;
}

# #######################################################################################
# check the inputs and create directory according to file structure
# #######################################################################################
print "HMM-GRASPx: Verifying directory structure...\n";
if(!(-e "$models"))  {
  die "The query HMM file <$models> cannot be found, please check your input.\n";
}
if(!(-e "$seq_db"))  {
  die "The short-peptide read file <$seq_db> cannot be found, please check your input.\n";
}
$hm_dir = abs_path($hm_dir);
if(! (-d "$hm_dir") || ! (-d "$hm_dir/bin") ||
    ! (-d "$hm_dir/scripts") || ! -d ("$hm_dir/thirdparty/hmmer-3")
)  {
  die "The home directory <$hm_dir> does not exist or the program is incomplete. Abort.\n";
}

if(! (-d "$idx_dir"))  {
    mkdir "$idx_dir" or die "Failed to create indexing directory <$idx_dir>. Abort.\n";
}
if(! (-d "$out_dir"))  {
  mkdir "$out_dir" or die "Failed to create output directory <$out_dir>. Abort.\n";
}
# #######################################################################################

# #######################################################################################
# parsing the parameters input by the sample sheet
# #######################################################################################
print "HMM-GRASPx: Resolving parameters for the run...\n";
my $param_num_threads = 1;
my $param_seed_len;		# intentionally left as un-defined
my $param_overlap_len;		# intentionally left as un-defiend
my $param_asm_depth = 5;
my $param_band_size = 20;
my $param_dup_reads = "false";
my $param_extend = 1;
my $param_graspxp_pv = 0.05;
my $param_hmmsearch_ev = 0.01;
my $param_num_errors = 3;
my $param_map_portion = 0.6;
# reads in the parameter sheet if provided by the user
if(defined $paramf)  {
  open my $PARIN, "<$paramf" or die "Cannot open parameter setting file <$paramf>: $!. Abort\n";
  while(<$PARIN>)  {
    chomp;
    if(/\S+/ and !(/^\#/))  {	# contains information and is not comment
      my @decom = split /\=/, $_;
      if($decom[0] eq 'NUM_THREADS')  {
        $param_num_threads = $decom[1];
      }  elsif($decom[0] eq 'SEED_LEN')  {
        $param_seed_len = $decom[1];
      }  elsif($decom[0] eq 'OVERLAP_LEN')  {
        $param_overlap_len = $decom[1];
      }  elsif($decom[0] eq 'ASSEMBLY_DEPTH')  {
        $param_asm_depth = $decom[1];
      }  elsif($decom[0] eq 'ALIGNMENT_BAND_SIZE')  {
        $param_band_size = $decom[1];
      }  elsif($decom[0] eq 'ALLOW_DUP_READS')  {
        $param_dup_reads = $decom[1];
      }  elsif($decom[0] eq 'PROGRESSIVE_EXT')  {
        $param_extend = $decom[1];
      }  elsif($decom[0] eq 'HMMGRASPX_PVALUE')  {
        $param_graspxp_pv = $decom[1];
      }  elsif($decom[0] eq 'HMMSEARCH_EVALUE')  {
        $param_hmmsearch_ev = $decom[1];
      }  elsif($decom[0] eq 'MAP_NUM_ERRORS')  {
        $param_num_errors = $decom[1];
      }  elsif($decom[0] eq 'MAP_LEN_PORTION')  {
        $param_map_portion = $decom[1];
      }  else  {
        print "HMM-GRASPx warning: Unrecognized parameter setting: $decom[0]. Skipped.\n";
      }
    }
  }
  close $PARIN;
}
# writing parameter settings to log file
open my $LOGOUT, ">$out_dir/hmm_graspx.log" 
  or die "Cannot write log file <$out_dir/hmm_graspx.log>: $!. Abort.\n";
print $LOGOUT "NUM_THREADS=$param_num_threads\n";
if(defined $param_seed_len)  {
  print $LOGOUT "SEED_LEN=$param_seed_len\n";
}  else  {
  print $LOGOUT "SEED_LEN=6\n";
}
if(defined $param_overlap_len)  {
  print $LOGOUT "OVERLAP_LEN=$param_overlap_len\n";
}  else  {
  print $LOGOUT "OVERLAP_LEN=10\n";
}
print $LOGOUT "ASSEMBLY_DEPTH=$param_asm_depth\n";
print $LOGOUT "ALIGNMENT_BAND_SIZE=$param_band_size\n";
print $LOGOUT "ALLOW_DUP_READS=$param_dup_reads\n";
print $LOGOUT "PROGRESSIVE_EXT=$param_extend\n";
print $LOGOUT "HMMGRASPX_PVALUE=$param_graspxp_pv\n";
print $LOGOUT "HMMSEARCH_EVALUE=$param_hmmsearch_ev\n";
print $LOGOUT "MAP_NUM_ERRORS=$param_num_errors\n";
print $LOGOUT "MAP_LEN_PORTION=$param_map_portion\n";
close $LOGOUT;
# #######################################################################################

# #######################################################################################
# check if we need to run the indexing step of HMM-GRASPx
# #######################################################################################
my $seq_db_base = basename($seq_db);
if(! (-e "$idx_dir/$seq_db_base.gsa") || ! (-e "$idx_dir/$seq_db_base.hsm") ||
   ! (-e "$idx_dir/$seq_db_base.lcp") || ! (-e "$idx_dir/$seq_db_base.mcp") ||
   ! (-e "$idx_dir/$seq_db_base.sxt") || ! (-e "$idx_dir/$seq_db_base.rxt")
)  {
  print "HMM-GRASPx: Indexing file for the short-peptide database cannot be found...\n";
  print "HMM-GRASPx: Trying to build index...\n";
  my $run_tag = "";
  if(defined $param_seed_len)  {
    $run_tag .= "--seed_len=$param_seed_len ";
  }
  if(defined $param_overlap_len)  {
    $run_tag .= "--extension_len=$param_overlap_len"
  }
  system "$hm_dir/bin/graspx-build $seq_db --work_space=$idx_dir $run_tag";
}  elsif(defined $param_seed_len || defined $param_overlap_len)  {
  print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  print "HMM-GRASPx warning: Indexing files have been built (probably with a different parameter setting)\n";
  print "	Program continues by ignoring SEED_LEN and/or OVERLAP_LEN settings...\n";
  print "	Terminate the current program execution if you find re-indexing necessary...\n";
  print "	To re-build indexing files, please clear all files in $idx_dir and re-run the pipeline...\n";
  print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
}
# #######################################################################################

# #######################################################################################
# run the alignment/assembly module (graspxp-assemble) of HMM-GRASPx
# #######################################################################################
print "HMM-GRASPx: Searching/assembling short-peptide reads...\n";
system "$hm_dir/bin/graspxp-assemble $models $seq_db $out_dir/raw.contig.fa --work_space=$idx_dir --num_threads=$param_num_threads --assembly_depth=$param_asm_depth --band_size=$param_band_size --dup_reads=$param_dup_reads --progressive_extend=$param_extend --p_value=$param_graspxp_pv";  
# #######################################################################################

# #######################################################################################
# run hmmsearch (HMMER3) to recompute statistics and identify proper regions in contigs
# #######################################################################################
print "HMM-GRASPx: Recomputing statistics of the contigs using hmmsearch (HMMER3)...\n";
system "$hm_dir/thirdparty/hmmer-3/binaries/hmmsearch --domE $param_hmmsearch_ev $models $out_dir/raw.contig.fa >$out_dir/validate_hmm.results";
# #######################################################################################

# #######################################################################################
# parse the hmmsearch results and generating validated contigs
# #######################################################################################
print "HMM-GRASPx: Extracting validated contigs from hmmsearch (HMMER3) results...\n";
open my $HIN, "<$out_dir/validate_hmm.results"
  or die "Cannot open HMM result file <$out_dir/validate_hmm.results>: $!. Abort.\n";
my $fam_name;		# name of the family, like "AMP-binding", no space allowed
my $fam_id;		# ID of the family, like PF00010 etc.
my $target_name;	# the name of the reference HMM through refering which the contig is assembled
my %identified_regions;
while(<$HIN>)  {
  chomp;
  if(/^Query:\s+(\S+)/)  {
    $fam_name = $1;
  } elsif(/^Accession:\s+(PF\d+)\./)  {
    $fam_id = $1;
  } elsif(/^>>\s+(\S+)/)  {
    $target_name = $1;
  } elsif(/\!/ || /\?/)  {
    my @decom = split /\s+/, $_;
    my @info = ($fam_name, $fam_id, $decom[6], $decom[10], $decom[11]);
    if($decom[6] <= $param_hmmsearch_ev)  {
      push @{$identified_regions{$target_name}}, \@info;
    }
  }
}
close $HIN;

# loads in the sequence file
my %seq_hash;
open my $PIN, "<$out_dir/raw.contig.fa"
  or die "Cannot open raw contig file <$out_dir/raw.contig.fa>: $!. Abort\n";
while(<$PIN>)  {
  chomp;
  if(/^>(\S+)/)  {
    my $id = $1;
    my $line = <$PIN>;
    chomp $line;
    $seq_hash{$id} = $line;
  }
}
close $PIN;

sub GetOverlap($$$$)  {
  my $a = shift;
  my $b = shift;
  my $c = shift;
  my $d = shift;
  if($a <= $c and $b >= $d)  {
    return $d - $c + 1;
  }  elsif($c <= $a and $d >= $b)  {
    return $b - $a + 1;
  }  elsif($b < $c or $a > $d)  {
    return 0;
  }  elsif($a <= $c and $c <= $b and $b <= $d)  {
    return $b - $c + 1;
  }  elsif($c <= $a and $a <= $d and $d <= $b)  {
    return $d - $a + 1;
  }
  die "ERROR:   $a      $b      $c      $d\n";
}

# wring the validated contigs
open my $COUT, ">$out_dir/validated.contig.fa" 
  or die "Cannot write to file <$out_dir/validated.contig.fa>: $! Abort.\n";

my $n = 0;
foreach(keys %identified_regions)  {
  my $target_seq = $_;
  my @regions = @{$identified_regions{$target_seq}};
  @regions = sort {$a->[2] <=> $b->[2]} @regions;
  my @checked;
  foreach(@regions)  {
    my $region = $_;
    my $max_overlap = 0;
    for(my $i = 0; $i < scalar(@checked); ++ $i)  {
      my $overlap = GetOverlap($checked[$i][0], $checked[$i][1], $region->[3], $region->[4]);
      $max_overlap = $overlap if $overlap > $max_overlap;
    }
    if($max_overlap <= 10 and $region->[4] - $region->[3] + 1 >= 20)  {
      $target_seq =~ /\|\|(\S+)/;
      my $fam_name = $1;
      # if we allow ambiguous hit
      if(($unique and $fam_name eq $region->[0]) || !$unique)  {
        my $seq = substr($seq_hash{$target_seq}, $region->[3] - 1, $region->[4] - $region->[3] + 1);
        print $COUT ">contig_$n||$region->[0]||$region->[1]\n$seq\n";
        ++ $n;
      }
      my @r = ($region->[3], $region->[4]);
      push @checked, \@r;
    }
  }
}
close $COUT;
# #######################################################################################

# #######################################################################################
# mapping the original reads to the original contigs
# #######################################################################################
print "HMM-GRASPx: Mapping short-peptides to validated contigs...\n";
system "$hm_dir/bin/graspx-map $seq_db $out_dir/validated.contig.fa $out_dir/recruited.list --work_space=$idx_dir --num_errors=$param_num_errors --portion_mapped=$param_map_portion";
# #######################################################################################

# #######################################################################################
# reporting results
# #######################################################################################
print "HMM-GRASPx: Pipeline finished.\n";
print "HMM-GRASPx: Results can be found in:\n";
print "	<$out_dir/hmm_graspx.log>:		log file for this run\n";
print "	<$out_dir/raw.contig.fa>:		raw contigs generated by HMM-GRASPx\n";
print "	<$out_dir/validate_hmm.results>:	hmmsearch (HMMER3) results (for alignment viewing)\n";
print "	<$out_dir/validated.contig.fa>:		hmmsearch (HMMER3) validated contigs\n";
print "	<$out_dir/recruited.list>:		list of short-peptide reads recruited\n";
