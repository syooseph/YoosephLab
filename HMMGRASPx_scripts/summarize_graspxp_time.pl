#!/usr/bin/perl -w
use strict;

my $models = shift;	# the query hmm models, used to extract the profile length information
my $log_file = shift;	# expected to be the graspxp log file (or nohup.out)
#my $recruit = shift;	# graspxp output read recruitment file, used to normalize run time

# build the length mapping using protein family name as key
my %length_hash;
open my $IN, "<$models" or die "Cannot open file: $!\n";
my $name;
my $len;
while(<$IN>)  {
  chomp;
  my $line = $_;
  #print "$line\n";
  if($line =~ /^NAME\s+(\S+)/)  {
    $name = $1;
  }  elsif($line =~ /^LENG\s+(\S+)/)  {
    $len = $1;
    $length_hash{$name} = $len;
  }
}
close $IN;

#foreach(keys %length_hash)  {
#  print "$_	$length_hash{$_}\n";
#}
#die;

# summarize the runtime information
my %time_hash;
my %seed_hash;
open my $LIN, "<$log_file" or die "Cannot open file: $!\n";
while(<$LIN>)  {
  chomp;
  if(/^GRASPxp.*\[(.*)\].*\[(.*)\]/)  {
    my $name = $1;
    my $time = $2;
    my @decom = split /\:/, $time;
    $time_hash{$name} += $decom[0] * 3600 + $decom[1] * 60 + $decom[2];
  }  elsif(/^GRASPxp.*\[(.*)\].*\{(\d+)\} seeds identified/)  {
    $seed_hash{$1} = $2;
  }
}
close $LIN;

# count reads per query
#my %read_hash;
#open my $RIN, "<$recruit" or die "Cannot open file: $!\n";
#while(<$RIN>)  {
#  chomp;
#  my @decom = split /\s+/, $_;
#  my @decom2 = split /\|\|/, $decom[2];
#  ++ $read_hash{$decom2[1]};
#}
#close $RIN;

foreach(keys %time_hash)  {
  #print "$_\n";
  if(exists $length_hash{$_} and exists $seed_hash{$_})  {
    print "$_	$length_hash{$_}	$seed_hash{$_}	$time_hash{$_}\n";
  }
}
