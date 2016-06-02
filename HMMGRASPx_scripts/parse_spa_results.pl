#!/usr/bin/perl -w
use strict;

# parse the performance of SPA assembly for each homolog search tool


my $raw_reads = shift;			# the raw read file, used to get the lengths of the reads
my $hmmer_results = shift;		# the hmmsearch results that used to evaluate true contigs and domain
my $spa_dir = shift;			# the SPA directory, where the contig sequence and read place information is read

# load the read length information
my %len_hash;
open my $FIN, "<$raw_reads" or die "Cannot open file: $!\n";
my $n = 0;
while(<$FIN>)  {
  chomp;
  if(/^>/)  {
    my $seq = <$FIN>;
    chomp $seq;
    $len_hash{$n} = length($seq);
    ++ $n;
  }
}
close $FIN;

# parse the hmmer3 infomration
open my $HIN, "<$hmmer_results" or die "Cannot open file: $!\n";
my @identified_regions;
my %good_contigs;
my $cid;
while(<$HIN>)  {
  chomp;
  if(/^\>\>\s+(\S+)\s+/)  {
    $cid = $1;
  }
  if(/\!/ || /\?/) {
    my @decom = split /\s+/, $_;
    #if(!exists $domain_evalue{$cid} || $domain_evalue{$cid} > $decom[6])  {
    #  $domain_evalue{$cid} = $decom[6];
    #  $contig_begins{$cid} = $decom[10];
    #  $contig_ends{$cid} = $decom[11];
    #}
    my @info = ($cid, $decom[6], $decom[10], $decom[11]);
    if($decom[6] <= 0.01) {
      push @identified_regions, \@info;
      $good_contigs{$cid} = 1;
    }
  }
}
close $HIN;

#foreach(keys %domain_evalue)  {
#  print "$_	$domain_evalue{$_}	$contig_begins{$_}	$contig_ends{$_}\n";
#}

# use SPA recruitment results to compute benchmark
open my $CIN, "<$spa_dir/spa.fasta" or die "Cannot open file: $!\n";
my $total_contigs = 0;
while(<$CIN>)  {
  chomp;
  if(/^>/)  {
    ++ $total_contigs;
  }
} 
close $CIN;

sub GetOverlap($$$$)  {
  my $a = shift;	# begin of region 1
  my $b = shift;	# end of region 1
  my $c = shift;	# begin of region 2
  my $d = shift;	# end of region 2
  if($a > $d || $c > $b)  {
    # case of no overlap
    return 0;
  }  elsif($a <= $c and $b >= $d)  {
    # inclusion
    return $d - $c + 1;
  }  elsif($c <= $a and $d >= $b)  {
    return $b - $a + 1;
  } elsif($a <= $c and $b >= $c and $b <= $d)  {
    return $b - $c + 1;
  } elsif($c <= $a and $d >= $a and $d <= $b)  {
    return $d - $a + 1;
  } else  {
    die "Unexpected configuration! <$a, $b>:<$c, $d>\n";
  }
}

open my $PIN, "<$spa_dir/spa.place.txt" or die "Cannot open file : $!\n";
my %good_reads;
while(<$PIN>)  {
  if(/^ID\:(\d+)/)  {
    my $cid = 'spa' . $1;
    #next if !exists $domain_evalue{$cid};
    # skip formatting info
    for(my $i = 0; $i < 3; ++ $i)  {
      my $tmp = <$PIN>;
    }
    while(1)  {
      my $info = <$PIN>;
      last if ($info =~ /^\/\//);
      my @decom = split /\s+/, $info;
      my $ratio = 0;
      for(my $k = 0; $k < scalar(@identified_regions); ++ $k)  {
        if($cid eq $identified_regions[$k][0])  {
          my $overlap = GetOverlap($identified_regions[$k][2], $identified_regions[$k][3], $decom[3], $decom[3] + $len_hash{$decom[1]} - 1);
          $ratio = $overlap / $len_hash{$decom[1]};
          last if $ratio >= 0.6;
        }
      }
      #my $overlap = GetOverlap($contig_begins{$cid}, $contig_ends{$cid}, $decom[3], $decom[3] + $len_hash{$decom[1]} - 1);
      #my $ratio = $overlap / $len_hash{$decom[1]};
      #print "$ratio	$overlap	$contig_begins{$cid}, $contig_ends{$cid}, $decom[3], $decom[3] + $len_hash{$decom[1]} - 1\n";
      $good_reads{$decom[1]} = 1 if ($ratio > 0);
    }
  }
}
close $PIN;
my $num_good = scalar(keys %good_reads);
my $num_total = scalar(keys %len_hash);
my $good_ratio = $num_good / $num_total;

print "#good reads:	$num_good\n";
print "#good ratio:	$good_ratio\n";

my $num_good_contigs = scalar(keys %good_contigs);
my $good_contig_ratio = $num_good_contigs / $total_contigs;

print "#good contigs:	$num_good_contigs\n";
print "#good contig ratio:	$good_contig_ratio\n";
