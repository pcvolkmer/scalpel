package Utils;

###################################################################
# Utils
#
# Package with basic commons routines
# 
#  Author: Giuseppe Narzisi 
#    Date: December 11, 2013
#
###################################################################

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(runCmd leftNormalize uniq binarySearch inTarget bychrpos elapsedTime fisher_yates_shuffle Log10 genotype);
@EXPORT_OK = qw($findVariants $findDenovos $findSomatic $exportTool $bamtools $samtools $bcftools);

use strict;
use warnings;
use POSIX;
use FindBin qw($Bin);
use lib $Bin; # add $Bin directory to @INC

# programs via absolute path
our $findVariants = "$Bin/FindVariants.pl";
our $findDenovos  = "$Bin/FindDenovos.pl";
our $findSomatic  = "$Bin/FindSomatic.pl";
#our $exportTool   = "$Bin/ExportVariants.pl";
our $exportTool   = "$Bin/scalpel-export";
our $bamtools     = "$Bin/bamtools-2.3.0/bin/bamtools";
our $samtools     = "$Bin/samtools-1.1/samtools";
our $bcftools     = "$Bin/bcftools-1.1/bcftools";

# Run system command 
#####################################################

sub runCmd 
{ 
	my $info = $_[0];
	my $cmd  = $_[1];
	#my $VERBOSE = $_[2];
	
	my $retval = 0;

	#print STDERR "$info ($cmd)...\n" if $VERBOSE;

	my $rc = system($cmd); 
	#die "failed ($rc)\n" if $rc; 
	if ($rc) { 
		print STDERR "Command failure: $info ($cmd)...\n";
		
		$retval = -1;
		
		if ($? == -1) {
			print "failed to execute: $!\n";
		}
		elsif ($? & 127) {
			printf "child died with signal %d, %s coredump\n",
			($? & 127),  ($? & 128) ? 'with' : 'without';
		}
		else {
			printf "child exited with value %d\n", $? >> 8;
		}	
	}
	return $retval;
}

## Normalize indel location to leftmost position
#####################################################

sub leftNormalize {
	my $indel = ${$_[0]};
	my $delta = $_[1];
	my $REF = $_[2];
	my $FAIDX = $_[3];
	my $genome = $_[4];
	
	my $pos = $indel->{pos};
	my $type = $indel->{type};
	my $chr = $indel->{chr};
	my $indel_len = $indel->{len};

	## extract sequence around indel from reference file
	#pos = $pos+1; # convert from 0-based to 1-based coordinate system
	my $l = $pos-$indel_len-$delta;
 	my $r = $pos+$indel_len+$delta;
		
	my $reference = "";
	if ($FAIDX == 0) {
		$reference = substr($genome->{$chr}->{seq}, $l-1, $r-$l+1);
	}
	else {
		$reference = readpipe( "$samtools faidx $REF $chr:$l-$r | awk '\$0 !~ /^>/'" ); # 0-based coordinate system
	}
	$reference=~s/\R//g; #remove newlines
	#print "samtools faidx $REF $chr:$L-$R\n";
	#print STDERR "$reference\n";
 	my $p = $indel_len+$delta;
	my $new_pos = $pos;
	if($type eq "del") { # deletion
		my $left_seq  = substr($reference, 0, $p);
		my $right_seq = substr($reference, $p+$indel_len);
		my $true_indel_seq = $left_seq . $right_seq;
		$new_pos = $l+$p;
		$indel->{ref} = substr($reference, $p, $indel_len);
		
		## reposition the candidate variant to a canonical (leftmost) position
		my $left_hyplotype;
		my $right_hyplotype;
		my $new_indel_seq;
		for (my $sft = 1; $sft <= $delta; $sft++) {
			my $i = $p-$sft;
			$left_hyplotype  = substr($reference, 0, $i);
		  	$right_hyplotype = substr($reference, $i+$indel_len);
		  	$new_indel_seq = $left_hyplotype . $right_hyplotype;
		  	#print "$true_indel_seq";
		    #print "$new_indel_seq";
		   	if ($true_indel_seq eq $new_indel_seq) {
				#print "Repositioning!\n";
		  		$new_pos = $l+$i; # repositioning!
		  		$indel->{ref} = substr($reference, $i, $indel_len);
		  	}
		}
		
		if($new_pos != $pos) {
			#$num_indels_repositioned++;
			my $shift = $pos - $new_pos + 1;
			#print STDERR "Deletion ($chr:$pos) left-shifted of $shift bp\n";
		}	
		
		$indel->{seq} = ('-' x $indel_len);
		$indel->{pos} = $new_pos;
		#$indel->{end} = $indel->{start} + $indel->{len} - 1;
	}
	elsif($type eq "ins") {
	  	my $left_seq  = substr($reference, 0, $p);
		my $right_seq = substr($reference, $p);
		my $true_indel_seq = $left_seq . $indel->{seq} . $right_seq;
		$new_pos = $l+$p;
		
		## reposition the candidate variant to a canonical (leftmost) position
		my $left_hyplotype;
		my $right_hyplotype;
		my $new_indel_seq;
		my $ins;
		for (my $sft = 1; $sft <= $delta; $sft++) {
			my $i = $p-$sft;
			$left_hyplotype  = substr($reference, 0, $i);
		  	$right_hyplotype = substr($reference, $i);
		  	$ins = substr($true_indel_seq, $i, $indel_len);
		  	$new_indel_seq = $left_hyplotype . $ins . $right_hyplotype;
				
		  	#print "$true_indel_seq";
		    #print "$new_indel_seq";
		  	if ($true_indel_seq eq $new_indel_seq) {
		  		#print "Repositioning!\n";
		  		$new_pos = $l+$i; # repositioning!
		  		$indel->{seq} = $ins;
		  	}
		}
		
		if($new_pos != $pos) {
			#$num_indels_repositioned++;
			my $shift = $pos - $new_pos + 1;
			#print STDERR "Insertion ($chr:$pos) left-shifted of $shift bp\n";
		}
		
		$indel->{ref} = ('-' x $indel_len);
		$indel->{pos} = $new_pos;
		#$indel->{end} = $indel->{start};
	}
}

## return unique elements of input list
#####################################################

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}

# binary-search
##########################################
sub binarySearch {

	my $x = $_[0];	
	my $A = $_[1];
	my $l = $_[2];
	my $r = $_[3];
	
	while($l<=$r) {
		my $m = floor( ($l+$r)/2 );
		my $f = @$A[$m];
		if($f->{end} < $x) { $l = $m+1; }
		elsif($f->{start} > $x)  { $r = $m-1; }
		else { 
			return $m; 
		}
	}
	return -1;
}

## return true if mutation is in target (start or end
## positions intersect the target), false otherwise
#####################################################

sub inTarget {
	my $mutation = $_[0];
	my $exons = $_[1];
	
	my $chr = $mutation->{chr};
	my $pos = $mutation->{pos};
	my $len = $mutation->{len};
	
	my $result = "false";
	
	my @array = @{$exons->{$chr}};
	
	my $flag1 = binarySearch($pos, \@array, 0, (scalar @array)-1);
	my $flag2 = binarySearch($pos+$len, \@array, 0, (scalar @array)-1);
	
	if( ($flag1 != -1) || ($flag2 != -1) ) { $result = "true"; }
	
	#foreach my $exon (@{$exons->{$chr}}) {
		
	#	my $s = $exon->{start};
	#	my $e = $exon->{end};
		
	#	if( ($s <= $pos) && ($pos <= $e) ) { 
	#		$result = "true"; 
	#		last;
	#	}
	#}
	
	return $result;
}

# auxilary routine to sort variants (Lexicographically)
# by chromosome and then by position
##########################################

sub bychrpos {
	my ($a, $b, $h) = @_;
	my $result;
	if ( ($h->{$a}->{chr} =~ m/^-?\d+\z/) && ($h->{$b}->{chr} =~ m/^-?\d+\z/) ) {
		$result = ( ($h->{$a}->{chr} <=> $h->{$b}->{chr}) || ($h->{$a}->{pos} <=> $h->{$b}->{pos}) );
	}
	else {
		$result = ( ($h->{$a}->{chr} cmp $h->{$b}->{chr}) || ($h->{$a}->{pos} <=> $h->{$b}->{pos}) );
	}
	return $result;
}

# print total time elapsed 
##########################################

sub elapsedTime {
	my $time_taken = $_[0];
	my $tool = $_[1];

	my $hours = $time_taken / 3600;
	my $days = $hours / 24; 
	$hours = $hours % 24;
	my $seconds = $time_taken % 3600;
	my $minutes = $seconds / 60;
	$seconds    = $seconds % 60;

	$days = floor($days);
	$hours = floor($hours);
	$minutes = floor($minutes);
	$seconds = floor($seconds);

	print STDERR "$tool elapsed time: $days day(s) $hours hour(s) $minutes minute(s) $seconds second(s)\n";
}



# compute genotype info in VCF format (GT field)
##########################################

#0/0 - the sample is homozygous reference
#0/1 - the sample is heterozygous, carrying 1 copy of each of the REF and ALT alleles
#1/1 - the sample is homozygous alternate
#1/2 - the sample is heterozygous, carrying 1 copy of the REF and 2 ALT alleles

sub genotype {
	my $R = $_[0]; # reference coverage
	my $A = $_[1]; # alternative coverage
	my $O = $_[2]; # other coverage (>0 for multiple alleles)
	my $Z = $_[3]; # zygosity (het/hom)
	my $GT = $_[4]; # GT field on VCF format
	
	if($R>0) {
		if($A>0 && $O==0) { $Z = "het"; $GT = "0/1"; }
		elsif($A>0 && $O>0) { $Z = "het"; $GT = "0/1"; }
		elsif($A==0 && $O==0) { $Z = "hom"; $GT = "0/0"; }
		elsif($A==0 && $O>0) { $Z = "hom"; $GT = "0/1"; }
	}
	elsif($R==0) { 
		if($A>0 && $O==0) { $Z = "het"; $GT = "1/1"; }
		elsif($A>0 && $O>0) { $Z = "het"; $GT = "1/1"; }
		elsif($A==0 && $O==0) { $Z = "hom"; $GT = "?"; }
		elsif($A==0 && $O>0) { $Z = "hom"; $GT = "1/1"; }
	}
	return ($Z, $GT);
}


# fisher_yates_shuffle( \@array ) : 
# generate a random permutation of @array in place
##########################################

sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}

# log 10
##########################################

sub Log10 {
        my $n = shift;
        return log($n)/log(10);
}

1;