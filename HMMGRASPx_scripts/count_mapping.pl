#!/usr/bin/perl -w
use strict;

my $contig = shift;	# the fasta file that defines the transcript
my $folder = shift;	# the folder contains the mapping results of GRASPx

# read in contig names
my %hash;
my %lenhash;
open my $IN, "<$contig" or die "Cannot open contig file: $!\n";
while(<$IN>)  {
  chomp;
  if(/^>(.*)/)  {
    $hash{$1} = 0;
    my $tmp = <$IN>;
    $lenhash{$1} = length($tmp) - 1;
  }
}
close $IN;

chdir "$folder" or die "Cannot change directory: $!\n";
my @hashes;
my @idname;
my @mapped;
my $index = 0;
foreach(<*.graspx.cdhit.e-10.map>)  {
  #print "loading $_...\n";
  push @idname, $_;
  open my $AIN, "<$_" or die "Cannot open map file: $!\n";
  foreach(keys %hash)  {
    $hashes[$index]->{$_} = 0;
  }
  #my %lhash;
  while(<$AIN>)  {
    chomp;
    my @decom = split /\s+/, $_;
    my @decom2 = split /\|\|/, $decom[2];
    if(exists $hashes[$index]->{$decom2[0]})  {
      ++ $hashes[$index]->{$decom2[0]};
      ++ $mapped[$index];
      if($hashes[$index]->{$decom2[0]} > $hash{$decom2[0]})  {
        $hash{$decom2[0]} = $hashes[$index]->{$decom2[0]};
      }
    }
  }
  close $AIN;
  ++ $index;
}

print "contig_id	";
foreach(@idname)  {
  print "$_	";
}
print "\n";

sub by_value  {
  $hash{$b} <=> $hash{$a};
}

my @normalize =  qw(10480524 8740916 11104094 11920968 26800436 6617138 14978858 18591974);
$index = 0;
foreach(sort by_value keys %hash)  {
  if($hash{$_} <= 0)  {
    last;
  }
  print "$_	";
  my $id = $_;
  for(my $i = 0; $i < scalar(@hashes); ++ $i) {
    #my $a = 1000000000 * $hashes[$i]->{$id} / ($mapped[$i] * $lenhash{$id});
    #print "$a	";
    my $a = int($hashes[$i]->{$id});
    print "$a	";
  }
  print "\n";
  ++ $index;
}
