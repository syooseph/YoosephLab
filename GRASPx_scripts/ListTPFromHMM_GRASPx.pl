#!/usr/bin/perl -w
use strict;
use Cwd;

my $pfamfile = shift;
my $contigfile = shift;
my $mapfile = shift;
my $outfolder = shift;

# loads in the pre-defined Pfam annotation for the genome
# returns the reference to a two-dimensional hash table
sub LoadPfamAnnot($)  {
  my $file = shift;
  my %hash;
  open my $IN, "<$file" or die "Cannot open file: $!\n";
  while(<$IN>)  {
    chomp;
    /(\S+)\s+pf\:(\S+)/;
    $hash{$2}{$1} = 1;
  }
  close $IN;
  return \%hash;
}

# loads in each contig and parse their E-values
# returns the reference to a hash table that maps contig ID to E-value
sub LoadContigEvalue($)  {
  my $file = shift;
  my %hash;
  open my $IN, "<$file" or die "Cannot open file: $!\n";
  while(<$IN>)  {
    chomp;
    if(/^>/)  {
      my @decom = split /\:\:/, $_;
      my $id = substr($_, 1, length($_) - 1);
      $hash{$id} = $decom[2];
    }
  }
  close $IN;
  return \%hash;
}

# loads in the mapping file
# returns the reference to a hash table that counts number of reads per contig
sub LoadMapCount($)  {
  my $file = shift;
  my %hash;
  open my $IN, "<$file" or die "Cannot open file: $!\n";
  while(<$IN>)  {
    chomp;
    my @decom = split /\s+/, $_;
    $hash{$decom[2]} ++;
  }
  close $IN;
  return \%hash;
}

# check the number of True Positive contigs from the given folder
# returns the reference to a hash table that contains the contig IDs
sub CheckHMMTPFolder($$)  {
  my $annot_file = shift;
  my $folder = shift;
  # load pfam annotation and construct hash table
  my $annot_hash_ref = LoadPfamAnnot($annot_file);
  # check each file in the folder
  my %tp_hash;
  my $current_dir = getcwd();
  chdir "$folder" or die "Cannot change directory: $!\n";
  foreach(<*.hmm.search>)  {
    /(\S+)\.hmm\.search/;
    my $pfam_id = $1;
    open my $IN, "<$_" or die "Cannot open file: $!\n";
    while(<$IN>)  {
      if(/^\s+\-\-\- full sequence \-\-\-/)  {
        my $tmp = <$IN>; $tmp = <$IN>; # skip two lines
        while(1)  {
          $tmp = <$IN>; chomp $tmp;
          #print "$tmp\n";
          if($tmp =~ /\-\-\-/ || !($tmp =~ /\S+/))  {
            last;
          }
          #if($tmp =~ /\-\-\-/)  {
          #  next;
          #}
          my @decom = split /\s+/, $tmp;
          $decom[9] =~ /^(\S+)\|\|/;
          if(exists $annot_hash_ref->{$pfam_id}{$1})  {
            $tp_hash{$decom[9]} = 1;
          }
        }
        last;
      }
    }
    close $IN;
  }
  chdir "$current_dir" or die "Cannot change directory: $!\n";
  return \%tp_hash;
}

my $tp_hash = CheckHMMTPFolder($pfamfile, $outfolder);
my $cg_hash = LoadContigEvalue($contigfile);
my $nt_hash = LoadMapCount($mapfile);

#foreach(keys %{$nt_hash}) {
#  print "$_	$nt_hash->{$_}\n";
#}
#exit;

my @num_tp;
my @num_all;
my @ecutoff = qw(1e-9 1e-7 1e-5 1e-3 1e-1 1e1);
foreach(keys %{$tp_hash})  {
  my $rid = $_;
  for(my $i = 0; $i < scalar(@ecutoff); ++ $i)  {
    if($cg_hash->{$rid} <= $ecutoff[$i] and exists $nt_hash->{$rid})  {
      $num_tp[$i] += $nt_hash->{$rid};
      # $num_tp[$i] ++;	# for non-grasp results
    }
  }
}

foreach(keys %{$cg_hash})  {
  my $rid = $_;
  for(my $i = 0; $i < scalar(@ecutoff); ++ $i)  {
    if($cg_hash->{$rid} <= $ecutoff[$i] and exists $nt_hash->{$rid})  {
      $num_all[$i] += $nt_hash->{$rid};
      #$num_all[$i] ++;	# for non-grasp results
    }
  }
}

for(my $k = 0; $k < scalar(@ecutoff); ++ $k)  {
  my $a = $num_tp[$k] / $num_all[$k];
  print "$num_tp[$k]	$num_all[$k]	$a\n";
}
