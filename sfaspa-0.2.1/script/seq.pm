#!/usr/bin/env perl

BEGIN { push @INC, `echo -n \$SPA_SCRIPT` }

################################################################################
#
#
#
################################################################################

use strict;
use warnings;
use Set::Scalar;
use io;

package seq;


sub loadBlast {
	my ( $aln_file, $format ) = @_;
	
	return _load_tab_file( $aln_file ) if $format == 8;
}

##==============================================================================
## Input: BLAST table format file
## Output: alignment array
##==============================================================================
sub _load_tab_file {
	my ($aln_file) = @_;

	my $fh = io::openInput($aln_file);
	my @alignment = <$fh>; 
	chomp @alignment;
	close $fh;
	
	return @alignment;
}

##==============================================================================
## Input: HMMER v3.0 tab-delimited output
## Output: extract domain information
##==============================================================================
sub loadHmmer3 {
	my ( $file, $domain ) = @_;
	my $fh = io::openInput($file);
	return grep( /$domain/, <$fh> );
}


##==============================================================================
## Input: FASTA sequences (STDIN, ACSII, GZ, or BZIP)
## Output: Hash table (key:sequence ID, value: sequence)
##==============================================================================
sub loadFasta {
	my ($seq_file) = @_;
	my ( %seqs, $id, $seq );

	my $fh = io::openInput($seq_file);
	while ( <$fh> ) {
		chomp;
		if ( /^>([^\s]+)/ ) {
			$seqs{$id} = $seq if defined $seq;
			$id = $1; $seq = "";
		} else {
			$seq .= $_;
		}
	}
	$seqs{$id} = $seq; 
	close($fh);

	return %seqs;
}

##==============================================================================
## Input: FASTA sequences (STDIN, ACSII, GZ, or BZIP)
## Output: Arrays of sequence IDs and sequences
##==============================================================================
sub loadFastaIdsAndSeqs {
	my ($seq_file) = @_;
	my ( @Ids, @Seqs, $id, $seq );

	my $fh = io::openInput($seq_file);
	while ( <$fh> ) {
		chomp;
		if ( /^>([^\s]+)/ ) {
			if ( defined $seq ) { 
				push @Ids,  $id;
				push @Seqs, $seq;
			}
			$id = $1; $seq = "";
		} else {
			$seq .= $_;
		}
	}
	push @Ids,  $id;
	push @Seqs, $seq;

	close($fh);
	return ( \@Ids, \@Seqs );
}

sub loadFastaIds {
	my ($seq_file) = @_;
	my ( @Ids, $id );

	my $fh = io::openInput($seq_file);
	while ( <$fh> ) {
		chomp;
		if ( /^>([^\s]+)/ ) {
			push @Ids,  $1;
		}
	}
	close($fh);
	return @Ids;
}

##==============================================================================
## Input : FASTA sequences
## Output: A hash table (key:sequence, value:array of IDs)
##==============================================================================
sub __sequence2IdTable {
	my ($seq_file) = @_;
	my ( %seqs, $id, $seq );

	my $fh = io::openInput($seq_file);
	while ( <$fh> ) {
		chomp;
		if ( /^>([^\s]+)/ ) {
			push @{$seqs{$seq}}, $id if defined $seq;
			$id = $1; $seq = "";
		} else {
			$seq .= $_;
		}
	}
	push @{$seqs{$seq}}, $id;
	close($fh);
	return %seqs;
}

##==============================================================================
## Input : FASTA sequences
## Output: a hash table of unique sequences (key:sequence, value:sequence ids)
##==============================================================================
sub uniqueSequences {
	my ($seq_file) = @_;
	my %Seqs = __sequence2IdTable($seq_file);

	my @Dups;
	foreach my $seq ( keys %Seqs ) {
		push @Dups, $seq if @{$Seqs{$seq}} > 1;
	}
	delete @Seqs{@Dups};
	return %Seqs;
}


##==============================================================================
## Input : FASTA sequences
## Output: A hash table of redundant sequences (key:sequence, value:IDs)
##==============================================================================
sub redundantSequences {
	my ($seq_file) = @_;
	my %Seqs = __sequence2IdTable($seq_file);

	my @Uniq;
	foreach my $seq ( keys %Seqs ) {
		push @Uniq, $seq if @{$Seqs{$seq}} == 1;
	}
	delete @Seqs{@Uniq};
	return %Seqs;
}


##==============================================================================
## Input : A hash table of FASTA sequence and optional output file name
## OUtput: FASTA sequence 
## Note  : Pass hash by reference to avoid STDIN reopen error 
##   e.g.: seq::writeFasta(\%Seqs);
##==============================================================================
sub writeFasta {
	my ( $Seqs, $out_file ) = @_;

	my $ofh = io::openOutput($out_file);
	foreach my $sid ( sort keys %$Seqs ) {
		print $ofh join("\n", ">$sid", $Seqs->{$sid}), "\n";
	}
	close($ofh);
}

##==============================================================================
## Input: 
##     1. a set of sequence IDs
##     2. a sequence file
##     3. a output fileFASTA sequences (STDIN, ACSII, GZ, or BZIP)
## Output: a subset of sequences proveded with IDs
##==============================================================================
sub extractFastaBySeqId {
	my ( $Ids, $seq_file, $out_file, $memory, $mline, $column ) = @_;

	return __extractFastaFromRAM ($Ids, $seq_file, $out_file, $mline, $column) if $memory || ! defined $memory;
	return __extractFastaFromDisk($Ids, $seq_file, $out_file, $mline, $column);
}

sub __extractFastaFromRAM {
	my ( $Ids, $seq_file, $out_file, $mline, $column ) = @_;
	my %Seqs = loadFasta($seq_file);

	my $ofh = io::openOutput($out_file);
	foreach my $id ( sort keys %Seqs ) {
		#print $ofh join("\n", ">$id", $Seqs{$id}), "\n" if $Ids->has($id);
		__flushFastaSequence( $ofh, $id, $Seqs{$id}, $mline, $column ) if $Ids->has($id);
	}
	close($ofh);
}

sub __extractFastaFromDisk {
	my ( $Ids, $seq_file, $out_file, $mline, $column ) = @_;
	my ( $id, $oid, $seq );

	my $ofh = io::openOutput($out_file);
	my $ifh = io::openInput($seq_file);
	while ( <$ifh> ) {
		#chomp;
		if ( /^>([^\s]+)/ ) {
			chomp;
			#print $ofh join("\n", $oid, $seq), "\n" if defined $id && $Ids->has($id);
			if ( defined $id && $Ids->has($id) ) {
				__flushFastaSequence( $ofh, $oid, $seq, $mline, $column ) ;
			}
			$id = $1; $oid = $_; $seq = "";
		} else {
			$seq .= $_;
		}
	}
	#print $ofh join("\n", $oid, $seq), "\n" if $Ids->has($id);
	if ( defined $id && $Ids->has($id) ){
		__flushFastaSequence( $ofh, $oid, $seq, $mline, $column );
	}

	close($ifh);
	## omit close($ofh) intentionally
	## to avoid multiple close of STDOUT
}

sub __flushFastaSequence
{
	my ( $ofh, $sid, $seq, $mline, $column ) = @_;
	chomp $seq;

	print $ofh $sid, "\n";
	if ( ! defined $mline || $column < 1 ) {
		print $ofh $seq, "\n";
	} else {
		$seq =~ s/\n+//g; ## drop end of line characters.
		my $beg = 0;
		while () {
			if ( $beg+$column < length $seq ) {
				print $ofh substr( $seq, $beg, $column ), "\n";
				$beg += $column;
			}
			else { 
				print $ofh substr( $seq, $beg ), "\n";
				last;
			}
		}
	}
}


##==============================================================================
##
##==============================================================================
sub extractFastaByOrder {
	my ( $Snums, $seq_file, $out_file, $one_base, $mline, $column ) = @_;
	my ( $count, $id, $oid, $seq ) = ( 0, undef, undef, undef );
	
	$count -= 1 if !defined $one_base;

	my $ofh = io::openOutput($out_file);
	my $ifh = io::openInput($seq_file);
	while ( <$ifh> ) {
		#chomp;
		if ( /^>([^\s]+)/ ) {
			chomp;
			$count++ if defined $id;
			__flushFastaSequence( $ofh, $oid, $seq, $mline, $column ) if $Snums->has($count);
			$id = $1; $oid = $_; $seq = "";
		} else {
			$seq .= $_;
		}
	}

	$count++ if defined $id;
	__flushFastaSequence( $ofh, $oid, $seq, $mline, $column ) if $Snums->has($count);

	close($ifh);
	## omit close($ofh) intentionally
	## to avoid multiple close of STDOUT
}



sub dropFastaByOrder {
	my ( $Snums, $seq_file, $out_file, $one_base, $mline, $column ) = @_;
	my ( $count, $id, $oid, $seq ) = ( 0, undef, undef, undef );
	
	$count -= 1 if !defined $one_base;

	my $ofh = io::openOutput($out_file);
	my $ifh = io::openInput($seq_file);
	while ( <$ifh> ) {
		#chomp;
		if ( /^>([^\s]+)/ ) {
			chomp;
			$count++ if defined $id;
			__flushFastaSequence( $ofh, $oid, $seq, $mline, $column ) if ! $Snums->has($count);
			$id = $1; $oid = $_; $seq = "";
		} else {
			$seq .= $_;
		}
	}

	$count++ if defined $id;
	__flushFastaSequence( $ofh, $oid, $seq, $mline, $column ) if ! $Snums->has($count);

	close($ifh);
	## omit close($ofh) intentionally
	## to avoid multiple close of STDOUT
}


##==============================================================================
##
##==============================================================================
sub splitFasta {
	my ( $seq_file, $nseq ) = @_;
	
	my ($from, $last, $iter) = (0,0,0);
	my ( %Seqs, $id, $seq );

	my $ifh = io::openInput($seq_file);
	while ( <$ifh> ) {
		#chomp;
		#if ( /^>([^\s]+)/ ) {
		if ( /^>/ ) { ## use full ID field
			chomp;
			if ( defined $id && defined $seq ) {
				$Seqs{$id} = $seq;
			}
			if ( scalar keys %Seqs == $nseq ) {
				$from = $iter*$nseq + 1;
				$last = $iter*$nseq + $nseq;
				writeFasta( \%Seqs, "$seq_file.$from-$last" );
				%Seqs = ();
				$iter++;
			}
			$id = $'; $seq = "";
		} else {
			$seq .= $_;
		}
	}
	if ( defined $id && defined $seq ) {
		$Seqs{$id} = $seq; 
	}
	$from = $iter*$nseq + 1;
	$last = $iter*$nseq + (scalar keys %Seqs);
	writeFasta( \%Seqs, "$seq_file.$from-$last" );
	close($ifh);
}

sub splitFastaWithEvenSequences
{
	my ( $input, $npart, $outdir, $name, $mline, $column ) = @_;
	my @Ids   = loadFastaIds($input);
	my @Parts = divideByEvenSequences(\@Ids, $npart);
	for ( my $i = 0; $i < $npart; $i++ ) {
		my $output = $outdir."/".$name.$i;
		__extractFastaFromDisk(Set::Scalar->new(@{$Parts[$i]}), $input, $output, $mline, $column );
	}
}

sub splitFastaWithEvenSize
{
	my ( $input, $npart, $outdir, $name, $mline, $column ) = @_;
	my %Sizes = getSequenceLengthFromDisk($input);
	my @Ids   = loadFastaIds($input);
	my @Parts = divideByEvenSize(\%Sizes, \@Ids, $npart);
	
	for ( my $i = 0; $i < @Parts; $i++ ) {
		my $output = $outdir."/".$name.$i;
		__extractFastaFromDisk(Set::Scalar->new(@{$Parts[$i]}), $input, $output, $mline, $column );
	}
}

sub divideByEvenSequences
{
	my ( $Ids, $npart ) = @_;

	my @Parts;
	my @Apart; 
	my $binsize = @$Ids/$npart;
	my $nseq = 0;
	for ( my $i = 0; $i < @$Ids; $i++ ) {
		my $id = $Ids->[$i];
		if ( $nseq < $binsize ) {
			push @Apart, $id;
			$nseq++;
		}
		else {
			push @Parts, [@Apart];
			@Apart = ();
			push @Apart, $id;
			$nseq = 1;
		}
	}
	push @Parts, [@Apart];
	return @Parts;
}


sub divideByEvenSize
{
	my ( $Sizes, $Ids, $npart ) = @_;
	my $total = getTotalSequenceLength( values %$Sizes );
	my @Parts;
	my @Apart; 
	my $sum = 0;
	for ( my $i = 0; $i < @$Ids; $i++ ) {
		my $id = $Ids->[$i];
		if ( $sum <= $total/$npart ) {
			push @Apart, $id;
			$sum += $Sizes->{$id};
		} 
		else {
			push @Parts, [@Apart];
			@Apart = ();
			push @Apart, $id;
			$sum = $Sizes->{$id};
		}
	}
	push @Parts, [@Apart];
	return @Parts;
}


##==============================================================================
## Input: 
##     1. a set of sequence IDs
##     2. a sequence file
##     3. a output fileFASTA sequences (STDIN, ACSII, GZ, or BZIP)
## Output: a subset of sequences proveded with IDs
##==============================================================================
sub dropFasta {	
	my ( $Ids, $seq_file, $out_file, $mline, $column ) = @_;
	
	my ( $id, $oid, $seq );

	my $ofh = io::openOutput($out_file);
	my $ifh = io::openInput($seq_file);
	while ( <$ifh> ) {
		if ( /^>([^\s]+)/ ) {
			chomp;
			#print $ofh join("\n", $oid, $seq), "\n" if defined $id && ! $Ids->has($id);
			__flushFastaSequence( $ofh, $oid, $seq, $mline, $column ) if defined $id && !$Ids->has($id);
			$id = $1; $oid = $_; $seq = "";
		} else {
			$seq .= $_;
		}
	}
	__flushFastaSequence( $ofh, $oid, $seq, $mline, $column ) if ! $Ids->has($id);
	#print $ofh join("\n", $oid, $seq), "\n" if ! $Ids->has($id);

	close($ifh);
	## omit close($ofh) intentionally
	## to avoid multiple close of STDOUT
}

sub getSequenceLength {
	my ( $file ) = @_;
	my %Seqs = loadFasta($file);

	foreach my $id ( keys %Seqs ) {
		$Seqs{$id} = length $Seqs{$id};
	}
	return %Seqs;
}

sub getSequenceLengthFromDisk 
{
	my ( $file ) = @_;
	my %Sizes;
	my $inh = io::openInput($file);
	my ( $id, $seq );
	while ( <$inh> ) {
		chomp;	
		if ( /^>([^\s]+)/ ) {
			if ( defined $seq ) { 
				$Sizes{$id} = length $seq;
			}
			$id = $1; $seq = "";
		} else {
			$seq .= $_;
		}
	}
	$Sizes{$id} = length $seq;
	return %Sizes;
}

sub getTotalSequenceLength
{
	my @Sizes = @_;
	my $total = 0;
	foreach my $s ( @Sizes ) {
		$total += $s;
	}
	return $total;
}

##==============================================================================
## Input: FASTQ sequences (STDIN, ACSII, GZ, or BZIP)
## Output: Hash table (key:sequence ID, value: sequence & quality)
##==============================================================================
sub loadFastq {
	my ($seq_file) = @_;
	my ( %seqs, $id, $seq, $qual );

	my $fh = io::openInput($seq_file);
	while ( <$fh> ) {
		chomp;
		if ( /^@/ ) {
			if ( defined $seq && defined $qual ) {
				push @{$seqs{$id}}, $seq;
				push @{$seqs{$id}}, $qual;
			}
			$id = $';
			$seq = <$fh>; $seq =~ s/\s+$//;
		} #else {
		elsif ( /^\+/ ) {
			$qual = <$fh>; $qual =~ s/\s+$//;
		}
	}
	if ( !defined $seqs{$id} && defined $seq && defined $qual ) {
		push @{$seqs{$id}}, $seq;
		push @{$seqs{$id}}, $qual;
	}

	#$seqs{$id} = $seq; 
	#$seqs{$id} = \($seq, $qual) if defined $seq && defined $qual;
	#@{$seqs{$id}} = ($seq, $qual) if defined $seq && defined $qual;
	close($fh);

	return %seqs;
}

sub extractFastq {
	my ( $Ids, $seq_file, $out_file, $memory ) = @_;

	return __extractFastqFromRAM ($Ids, $seq_file, $out_file) if $memory || ! defined $memory;
	return __extractFastqFromDisk($Ids, $seq_file, $out_file);
}

sub __printFastqEntry {
	my ($ofh, $id, $seq, $qual) = @_;
	print $ofh join("", "@", $id), "\n";
	print $ofh $seq, "\n";
	#print $ofh join("", "+", $id), "\n";
	print $ofh "+\n";
	print $ofh $qual, "\n";
}

sub __extractFastqEntry {
	my ($ofh, $Seqs, $id) = @_;
	my $seq = $Seqs->{$id}->[0];
	my $qual = $Seqs->{$id}->[1];
	__printFastqEntry($ofh, $id, $seq, $qual);
}

sub __extractFastqFromRAM {
	my ( $Ids, $seq_file, $out_file ) = @_;
	my %Seqs = loadFastq($seq_file);

	my $ofh = io::openOutput($out_file);
	foreach my $id ( sort keys %Seqs ) {
		next unless $Ids->had($id);
		__extractFastqEntry($ofh, \%Seqs, $id);
	}
	close($ofh);
}

sub __extractFastqFromDisk {
	my ( $Ids, $seq_file, $out_file ) = @_;
	my ( $id, $seq, $aid, $qual );

	my $ofh = io::openOutput($out_file);
	my $ifh = io::openInput($seq_file);
	while ( <$ifh> ) {
		chomp;
		if ( /^@/ ) {
			$id = $'; $id =~ s/\s+$//;
			my $qid = $id;
			if ( $id =~ /\s+/ ) {
				$qid = $`;
			}
			$seq = <$ifh>; $seq =~ s/\s+$//;
			$aid = <$ifh>; ## id fieild for quality
			die "Format error:$id\t$aid" if $aid !~ /^\+/;
			$qual = <$ifh>; $qual =~ s/\s+$//;

			next unless $Ids->has($qid);

			__printFastqEntry($ofh, $id, $seq, $qual);
		}
	}
	close($ifh);
	## omit close($ofh) intentionally
	## to avoid multiple close of STDOUT
}

sub dropFastq {	
	my ( $Ids, $seq_file, $out_file ) = @_;
	
	my ( $id, $seq, $aid, $qual );

	my $ofh = io::openOutput($out_file);
	my $ifh = io::openInput($seq_file);
	while ( <$ifh> ) {
		chomp;
		if ( /^@/ ) {
			$id = $'; $id =~ s/\s+$//;
			my $qid = $id;
			if ( $id =~ /\s+/ ) {
				$qid = $`;
			}

			$seq = <$ifh>; $seq =~ s/\s+$//;
			$aid = <$ifh>; ## id fieild for quality
			die "Format error:$id\t$aid" if $aid !~ /^\+/;
			$qual = <$ifh>; $qual =~ s/\s+$//;

			next if $Ids->has($qid);

			__printFastqEntry($ofh, $id, $seq, $qual);
		}
	}
	close($ifh);
	## omit close($ofh) intentionally
	## to avoid multiple close of STDOUT
}


#==============================================================================
# Count k-mer occurrences.
# Add k-mer occurrences if there exists k-mer entry in a hash table
#------------------------------------------------------------------------------
# Input: 
#   1. k
#   2. sequence array
#   3. k-mer hash
# Output: k-mer hash table of k-mer frequency
#==============================================================================
sub countKmers {
	my ( $k, $Seqs, $Kmers ) = @_;
	foreach my $seq ( @$Seqs ) {
		$seq = uc $seq;
		$seq =~ s/[^A-Z]//g;
		my $i = 0;
		while ( my $kmer = substr( $seq, $i, $k ) ) {
			last if length $kmer < $k;
			
			$Kmers->{$kmer} = 0	 if ! defined $Kmers->{$kmer};
			$Kmers->{$kmer} += 1;
			$i++;
		}
	}
}

##==============================================================================
## Find k-mer occurrences in a sequence
## Append sequence ID and k-mer position to hash table
##------------------------------------------------------------------------------
## Input: 
##   1. k
##   2. sequence hash (key: sequece ID, value: sequence)
##   3. k-mer hash
## Output: k-mer hash table of k-mer positions
##==============================================================================
sub locateKmers {
	my ( $k, $Seqs, $Kmers ) = @_;

	foreach my $id ( keys %$Seqs ) {
		my $seq = $Seqs->{$id};
		$seq = uc $seq;
		$seq =~ s/[^A-Z]//g;

		my $i = 0;
		while ( my $kmer = substr( $seq, $i, $k ) ) {
			last if length $kmer < $k;
			
			push @{$Kmers->{$kmer}}, [ $id, $i ];
			$i++;
		}
	}
}

##==============================================================================
## Reverse complement
##==============================================================================
sub rc 
{
	my $seq = shift;
	$seq = reverse $seq;
	$seq =~ tr/ATUGCYRKMBDHV/TAACGRYMKVHDB/;
	return $seq;
}

1;
