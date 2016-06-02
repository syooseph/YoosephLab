#!/usr/bin/perl -w
use strict;

my $software = shift;		# can be either "FASTM", "BLASTP", or "PSIBLAST"
my $result_file = shift;	# the result file that contains all the identified sequence IDs
my $sequence_file = shift;	# the multi-fasta file that contains all sequences

# parsing the FASTM result file and find out the best hit of each sequence,
# and the corresponding E-values
sub GetBestHitInfoFASTM($)  {
  my $file = shift;
  open my $IN, "<$file" or die "Cannot open file: $!\n";
  my %hit_hash;
  my %evalue_hash;
  while(<$IN>)  {
    if(/^>>(\S+)\s+\(.*\)/)  {
      chomp;
      my $sid = $1;
      my $tmp = <$IN>;	# get the following row
      $tmp =~ /E\(.*\)\:\s+(\S+)/;
      my $evalue = $1;
      for(my $i = 0; $i < 4; ++ $i)  {
        $tmp = <$IN>;	# jump 4 rows
      }
      my @decom = split /\s+/, $tmp;
      my $hid = $decom[0];
      # record information
      if(!(exists $evalue_hash{$sid}) || $evalue_hash{$sid} > $evalue)  {
        $evalue_hash{$sid} = $evalue;
        $hit_hash{$sid} = $hid;
      }
    }
  }
  close $IN;
  return (\%hit_hash, \%evalue_hash);
}

# parsing the BLASTP result file and find out the best hit of each sequence,
# and the corresponding E-values
sub GetBestHitInfoBLASTP($)  {
  my $file = shift;
  my %hit_hash;
  my %evalue_hash;
  open my $IN, "<$file" or die "Cannot open file: $!\n";
  while(<$IN>)  {
    if(/^#/)  {
      next;
    }
    my @decom = split /\s+/, $_;
    if(scalar(@decom) < 12)  {
      next;
    }
    my $sid = $decom[1];
    # record information
    if(!(exists $evalue_hash{$sid}) || $evalue_hash{$sid} > $decom[10])  {
      $evalue_hash{$sid} = $decom[10];
      $hit_hash{$sid} = $decom[0];
    }
  }
  close $IN;
  return (\%hit_hash, \%evalue_hash);
}

# parsing the PSI-BLAST result file and find out the best hit of each sequence,
# # and the corresponding E-values
sub GetBestHitInfoPSIBLAST($)  {
  my $file = shift;
  my %hit_hash;
  my %evalue_hash;
  open my $IN, "<$file" or die "Cannot open file: $!\n";
  while(<$IN>)  {
    if(/^#/ || /^\/\//)  {
      next;
    }
    my @decom = split /\s+/, $_;
    if(scalar(@decom) < 12)  {
      next;
    }
    my $sid = $decom[1];
    if(!(exists $evalue_hash{$sid}) || $evalue_hash{$sid} > $decom[10])  {
      $evalue_hash{$sid} = $decom[10];
      $hit_hash{$sid} = $decom[0];
    }
  }
  close $IN;
  return (\%hit_hash, \%evalue_hash);
}

# loading the sequences from a fasta file and storing them in a hash table
sub LoadFASTASeqs($)  {
  my $file = shift;
  my %seq_hash;
  open my $IN, "<$file" or die "Cannot open file: $!\n";
  my $id;
  while(<$IN>)  {
    chomp;
    if(/^>(.*)/)  {
      $id = $1;
    }  else {
      $seq_hash{$id} .= $_;
    }
  }
  close $IN;
  return \%seq_hash;
}

my $h_ref;
my $e_ref;

if($software eq 'FASTM')  {
  ($h_ref, $e_ref) = GetBestHitInfoFASTM($result_file);
}  elsif($software eq 'BLASTP')  {
  ($h_ref, $e_ref) = GetBestHitInfoBLASTP($result_file);
}  elsif($software eq 'PSIBLAST')  {
  ($h_ref, $e_ref) = GetBestHitInfoPSIBLAST($result_file);
}  else  {
  die "Unsupported format.\n";
}

=pod
foreach(keys %{$h_ref})  {
  print "$_	$h_ref->{$_}	$e_ref->{$_}\n";
  my $a = length($_);
  print "$a\n";
}
=cut
#exit;
my $s_ref = LoadFASTASeqs($sequence_file);
=pod
foreach(keys %{$s_ref})  {
  print "$_	$s_ref->{$_}\n";
  my $a = length($_);
  print "$a\n";
}
=cut
#exit;

# format output
my $n = 0;
foreach(keys %{$h_ref})  {
  #print "$_	$h_ref->{$_}	$e_ref->{$_}	$s_ref->{$_}\n";
  if(exists $s_ref->{$_})  {
    print ">$h_ref->{$_}";
    print "||contig_$n";
    print "\:\:0-0\:\:";
    print "$e_ref->{$_}";
    print "\:\:0\:\:0\n";
    print "$s_ref->{$_}\n";
    ++ $n;
  }
}
