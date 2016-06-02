#!/usr/bin/perl -w
use strict;

my $folder = shift;		# the folder expected to contain all individual pooled mapping results

my @headers;
my %count_hash;
foreach(<$folder/*.summary>)  {
  #print "$_\n";
  my $file = $_;
  $file =~ /.+\/(.+)/;
  my $file_stem = $1;
  push @headers, $file_stem;
}
my @init_count;
for(my $i = 0; $i < scalar(@headers); ++ $i)  {
  $init_count[$i] = 0;
}

my $index = 0;
foreach(<$folder/*.summary>) {
  # analyze the file
  open my $IN, "<$_" or die "Cannot open file: $!\n";
  while(<$IN>)  {
    if(/^\#/)  {
      next;
    }
    chomp;
    my @decom = split /\s+/, $_;
    my $ukey = $decom[0];
    if(!exists $count_hash{$ukey})  {
      @{$count_hash{$ukey}} = @init_count;
    }
    ${$count_hash{$ukey}}[$index] = $decom[1];
  }
  close $IN;
  ++ $index;
}

print "TAXA_GENE\t";
foreach(my $i = 0; $i < scalar(@headers) - 1; ++ $i)  {
  print "$headers[$i]\t";
}
print "$headers[-1]";
print "\n";
foreach(sort keys %count_hash)  {
  print "$_\t";
  foreach(my $i = 0; $i < scalar(@{$count_hash{$_}}) - 1; ++ $i)  {
    print "${$count_hash{$_}}[$i]\t";
  }
  print ${$count_hash{$_}}[-1];
  print "\n";
}
