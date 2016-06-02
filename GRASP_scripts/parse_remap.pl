#!/usr/bin/perl -w
use strict;

my $file = shift;
open my $IN, "<$file" or die "Cannot open file: $!\n";

while(<$IN>)  {
  chomp;
  if(/^\#\s+Query:\s+(.*)/)  {
    my $read_ID = $1;
    $read_ID =~ /^(.*\|)\_(\d+)\_(\d+)/;
    my $genome = $1;
    my $begin = $2;
    my $end = $3;
    # check how many reads are present
    my $temp = <$IN>;
    $temp = <$IN>;
    my $hits;
    if($temp =~ /^\#\s+Fields/)  {
      $hits = <$IN>;
    }  else  {
      $hits = $temp;
    }  
    $hits =~ /^\#\s+(\d+)\s+hits/;
    
    # if we have more than one hit
    #print "$1 hits\n";
    my $num_hits = $1;
    if($num_hits > 0)  {
      while($num_hits > 0)  {
        -- $num_hits;
        my $line = <$IN>;
        chomp $line;
        #print "LINE $line\n";
        if($line =~ /^\#/)  {
          #print "quit\n";
          last;
        }
        my @decom = split /\s+/, $line;
        #print "$genome $decom[1]\n";
        if($genome eq $decom[1] && $decom[8] >= $begin && $decom[8] <= $end && $decom[9] >= $begin && $decom[9] <= $end)  {
          my $strand;
          if($decom[8] < $decom[9])  {
            $strand = '+';
          }  else  {
            $strand = '-';
          }
          print "$read_ID  $genome  $strand $decom[8]  $decom[9]\n";
          last;
        }
      }
    }
  }
}
close $IN;
