#!/usr/bin/env perl

BEGIN { push @INC, `echo -n \$SPA` }

use strict;
use warnings;
use Getopt::Long;
use Set::Scalar;
use Statistics::Lite qw(:all);
use POSIX qw(ceil);
use seq;

my ( $path_file, $gene_file, $lcutoff, $pcutoff, $k );
GetOptions( "path=s" => \$path_file,
			"gene=s" => \$gene_file,
			"length=f" => \$lcutoff,
			"cutoff=f" => \$pcutoff,
			"kmer=i" => \$k );

if ( !defined $path_file || !defined $gene_file ) {
	print "Usage:$0 -p path-file -g gene-file\n"; exit
}

$k = 3 if !defined $k;
$pcutoff = 0.8 if !defined $pcutoff;
$lcutoff = 0.5 if !defined $lcutoff;

my %paths = seq::loadFasta( $path_file );
my %genes = seq::loadFasta( $gene_file );

## remove neither amino-acid nor gap
foreach my $id ( keys %genes ) {
	$genes{$id} =~ s/[^\w-]//g;
}

my %Short;
foreach my $id ( keys %genes ) {
	$id =~ /_.*/;
	push @{$Short{$`}}, $id;
}


my $cid;
foreach my $id ( sort {$a<=>$b} keys %paths ) {
	#my $cid = $id;
	my (@Glens, @Strs, @Knum);

	my $plen = length($paths{$id});
	if ( !defined $Short{$id} ) {
		print STDERR join("\t", $id, $plen, 0), "\n";
		$cid = $id;
		next;
	}	
	
	foreach my $gid ( @{$Short{$id}} ) {
		my $strand = "+";
		#print join("\t", $id, $plen, scalar @{$Short{$id}}), "\t";
		$strand = "-" if ( $gid =~ /-$/ );
		
		my $glen = length($genes{$gid});
	
		## substring match
		my $match = ( $paths{$id} =~ /$genes{$gid}/ ) ? 1 : 0;
		
		my $count = 0;
		my $i = 0;
		while ( my $kmer = substr( $genes{$gid}, $i, 6 ) ) {
			last if length($kmer) < $k;
			$count++ if $paths{$id} =~ /$kmer/;
			$i++;
		}
		print STDERR join("\t", $id, $plen, scalar @{$Short{$id}}, $strand, $glen, $match, $count), "\n";
		push @Glens, $glen;
		push @Strs, $strand;
		push @Knum, $count;
			
	}

	if ( !defined $cid || $cid ne $id ) {
		$cid = $id;
		if ( valid( $plen, \@Glens, \@Strs, \@Knum ) ) {
			print "$cid\t1\n";
		} else {
			print "$cid\t0\n";
		}
		@Glens = @Strs = @Knum = ();
	}
}


sub valid 
{
	my ( $plen, $glens, $strs, $knums ) = @_;
	my $index = findMax($glens, $knums);

	my $filter = 0;
	for ( my $i = 0; $i < @$glens; $i++ ) {
		my $mink = minKmerCount( $glens->[$i] );
		if ( $knums->[$i] >= $mink ) {
			$filter = 1;
			last;
		}
	}
	
	return 0 if ! $filter;
	return 0 if $strs->[$index] eq '-';
	return 0 if sum(@$glens)/$plen < $lcutoff;
	return 1;
}

sub findMax
{
	my ( $lens, $knums ) = @_;
	my $max = 0;
	for ( my $i = 1; $i < @$lens; $i++ ) {
		if ($knums->[$i]/$lens->[$i] > $knums->[$max]/$lens->[$max] ) {
			$max = $i;
		}
	}
	return $max;
}

# sub findMax
# {
# 	my @lens = @_;
# 	my $max = 0;
# 	for ( my $i = 1; $i < @lens; $i++ ) {
# 		if ( $lens[$i] > $lens[$max] ) {
# 			$max = $i;
# 		}
# 	}
# 	return $max;
# }

sub minKmerCount
{
	my $seqlen = shift;
	my $mismatchs = ceil($seqlen*(1-$pcutoff));
	my $min_kmers = $seqlen - $k + 1 - $k * $mismatchs;
	return $min_kmers;
}
