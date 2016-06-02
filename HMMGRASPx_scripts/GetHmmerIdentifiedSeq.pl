#!/usr/bin/perl -w
use strict;

my $hmm_full_file = shift;		# the hmmsearch result file
my $peptide_file = shift;		# the target peptide file in multi-fasta format
my $boundary = 1;			# consider exact boundaries when fetching sequences

open my $IN, "<$hmm_full_file" or die "Cannot open file: $!\n";

my $name;
my $id;
my $target;
my %identified_regions;
while(<$IN>)  {
  chomp;
  if(/^Query:\s+(\S+)/)  {
    $name = $1;
  } elsif(/^Accession:\s+(PF\d+)\./)  {
    $id = $1;
  } elsif(/^>>\s+(\S+)/)  {
    $target = $1;
  } elsif(/\!/ || /\?/)  {
    my @decom = split /\s+/, $_;
    #my @target_IDs = split /\|\|/, $target;
    #if($target_IDs[1] eq $name)  {
    my @info = ($name, $id, $decom[6], $decom[10], $decom[11]);
    if($decom[6] <= 0.01)  {
      push @{$identified_regions{$target}}, \@info;
    }
      #print "$name $id $target $decom[10] $decom[11]\n";
    #}
  }
}

# load the sequence in
my %seq_hash;
open my $PIN, "<$peptide_file" or die "Cannot open file: $!\n";
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
  die "ERROR:	$a	$b	$c	$d\n";
}

# get the domain-annotated sequence
my $n = 0;
foreach(keys %identified_regions)  {
  my $target_seq = $_;
  #print "$target_seq\n";
  my @regions = @{$identified_regions{$target_seq}};
  @regions = sort {$a->[2] <=> $b->[2]} @regions;
  my @checked;
  foreach(@regions)  {
    my $region = $_;
    #print "!!!	$region->[0]	$region->[1]	$region->[2]	$region->[3]	$region->[4]\n";
    my $max_overlap = 0;
    for(my $i = 0; $i < scalar(@checked); ++ $i)  {
      my $overlap = GetOverlap($checked[$i][0], $checked[$i][1], $region->[3], $region->[4]);
      $max_overlap = $overlap if $overlap > $max_overlap;
    }
    # if no signficant overlap with the already identified region
    if($max_overlap <= 10 and $region->[4] - $region->[3] + 1 >= 20)  {
      # get the sequence
      my $seq = substr($seq_hash{$target_seq}, $region->[3] - 1, $region->[4] - $region->[3] + 1);
      print ">contig_$n||$region->[0]||$region->[1]\n$seq\n";
      ++ $n;
      # mark the region
      my @r = ($region->[3], $region->[4]);
      push @checked, \@r;
    }
  }
}


=pod
my $n = 0;
foreach(@identified_regions)  {
  my @info = @{$_};
  if(exists $seq_hash{$info[2]})  {
    my $seq;
    if($boundary)  {
      $seq = substr($seq_hash{$info[2]}, $info[3] - 1, $info[4] - $info[3] + 1);
    }  else  {
      $seq = $seq_hash{$info[2]};
    }
    print ">contig_$n||$info[0]||$info[1]\n$seq\n";
    ++ $n;
  }
}
=cut

