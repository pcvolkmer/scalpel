#!/usr/bin/env perl

###################################################################
# FindSomatic.pl
#
# Tool for detecting somatic mutations in normal/tumor pair 
# using microassembly
# 
#  Author: Giuseppe Narzisi 
#    Date: December 11, 2013
#
###################################################################

use warnings;
use strict;
use POSIX;
use FindBin qw($Bin);
use lib $Bin; # add $Bin directory to @INC
use Usage;
use Utils qw(:DEFAULT $findVariants $findDenovos $findSomatic $exportTool $bamtools $samtools $bcftools);
use HashesIO;
use SequenceIO;
use Parallel::ForkManager;
use List::Util qw[min max];
#use Math::Random qw(:all);
#use Sys::CPU;
#use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use File::Spec;
use Getopt::Long;

use MLDBM::Sync;                       # this gets the default, SDBM_File
use MLDBM qw(DB_File Storable);        # use Storable for serializing
use MLDBM qw(MLDBM::Sync::SDBM_File);  # use extended SDBM_File, handles values > 1024 bytes
use Fcntl qw(:DEFAULT);                # import symbols O_CREAT & O_RDWR for use with DBMs

$|=1;

# required arguments
my $BAMNORMAL;
my $BAMTUMOR;
my $BEDFILE;
my $REF;
my $USEFAIDX = 0;

my $defaults = getDefaults();  

# optional arguments
my $kmer = $defaults->{kmer};
my $cov_threshold = $defaults->{cov_threshold};
my $tip_cov_threshold = $defaults->{tip_cov_threshold};
my $radius = $defaults->{radius};
my $windowSize = $defaults->{windowSize};
my $max_reg_cov = $defaults->{max_reg_cov};
my $delta = $defaults->{delta};
my $map_qual = $defaults->{map_qual};
my $maxmismatch = $defaults->{maxmismatch};
my $min_cov = $defaults->{min_cov};
my $max_cov = $defaults->{max_cov};
my $covratio = $defaults->{covratio};
my $outratio = $defaults->{outratio};
my $WORK  = $defaults->{WORK};
my $MAX_PROCESSES = $defaults->{MAX_PROCESSES};
my $selected = $defaults->{selected};
my $dfs_limit = $defaults->{pathlimit};
my $outformat = $defaults->{format};
my $rand_seed = $defaults->{rand_seed};
my $maxNormalContamination = $defaults->{contaminant_fraction_normal};
my $intarget;
my $logs;
my $STcalling;

my $help = 0;
my $VERBOSE = 0;

my %candidates;
my %exons;
my %locations;
my %partitions;
my %genome;

my %normalCov;
my %tumorCov;

my %normalL2K;
my %tumorL2K;

my %normalSVs;
my %tumorSVs;
my %somaticSVs;
my %commonSVs;

my $tumor_name;
my $normal_name;

my %alnHashN; # hash of mutations from read alignments for normal
my %alnHashT; # hash of mutations from read alignments for tumor

my @members = qw/normal tumor/;

my $argcnt = scalar(@ARGV);
my $start_time = time;

#####################################################
#
# Command line options processing
#
GetOptions(
	'help!'    => \$help,
	'verbose!' => \$VERBOSE,
    
	# required parameters
    'normal=s' => \$BAMNORMAL,
    'tumor=s'  => \$BAMTUMOR,
    'bed=s'    => \$BEDFILE,
    'ref=s'    => \$REF,

	# optional paramters
    'kmer=i'       => \$kmer,
    'covthr=i'     => \$cov_threshold,
    'lowcov=i'     => \$tip_cov_threshold,
    'mincov=i'     => \$min_cov,
    'covratio=f'   => \$covratio,
    'radius=i'     => \$radius,
    'window=i'     => \$windowSize,
    'maxregcov=i'  => \$max_reg_cov,
    'step=i'       => \$delta,
    'mapscore=i'   => \$map_qual,
    'mismatches=i' => \$maxmismatch,
    'pathlimit=i'  => \$dfs_limit,
    'dir=s'        => \$WORK,
    'numprocs=i'   => \$MAX_PROCESSES,
    'coords=s'     => \$selected,
    'format=s'     => \$outformat,
    'seed=i'       => \$rand_seed,
    #'normalContamination=f' => \$maxNormalContamination,

	# ouptut parameters
	'intarget!'    => \$intarget,
	'logs!'   	   => \$logs,
	'STcalling!'   => \$STcalling,
	'outratio=f'   => \$outratio,

) or usageSomatic("FindSomatic.pl");

#####################################################
#
# Command line options processing
#
sub init()
{
	usageSomatic("FindSomatic.pl") if ($argcnt < 1);
	usageSomatic("FindSomatic.pl") if( $help );
	
	if( (!defined $BAMNORMAL) || (!defined $BAMTUMOR) || (!defined $BEDFILE) || (!defined $REF) ) { 
		print STDERR "Required parameter missing!\n";
		usageSomatic("FindSomatic.pl"); 
	}
	
	if (-e "$BEDFILE") { $USEFAIDX = 0; }
	elsif ($BEDFILE =~ m/(\w+):(\d+)-(\d+)/) {
		#parse region
		my ($chr,$interval) = split(':',$BEDFILE);
		my ($start,$end) = split('-',$interval);
		
		# use samtools faidx for regions smaller than 1Mb to speed up loading of the sequences
		if( ($end-$start+1) < 10000) { 
			#print STDERR "INFO: Region $bedfile is < 10Kb; Using samtools faidx to extract sequences.\n";
			$USEFAIDX = 1;
		}
		else { $USEFAIDX = 0; } 
	}
	
	mkdir $WORK if ! -r $WORK;
	
	if($MAX_PROCESSES == 1) { $MAX_PROCESSES = 0; }
}

#####################################################

sub printParams { 

	my $file = $_[0];
	
	print STDERR "-- Print parameters to $file\n";
	
	open PFILE, "> $file" or die "Can't open $file ($!)\n";
	
	print PFILE "Parameters: \n";
	print PFILE "-- K-mer size (bp): $kmer\n";
	print PFILE "-- coverage threshold: $cov_threshold\n";
	print PFILE "-- low coverage threshold: $tip_cov_threshold\n";
	#print PFILE "-- normal contamination fraction: $maxNormalContamination\n";
	print PFILE "-- size (bp) of the left and right extensions (radius): $radius\n";
	print PFILE "-- window-size size (bp): $windowSize\n";
	print PFILE "-- max coverage per region: $max_reg_cov\n";
	print PFILE "-- step size (bp): $delta\n";
	print PFILE "-- minimum mapping-quality: $map_qual\n";
	print PFILE "-- minimum coverage for exporting somatic mutation: $min_cov\n";
	print PFILE "-- minimum coverage ratio for sequencing error: $covratio\n";
	print PFILE "-- minimum coverage ratio for exporting mutation: $outratio\n";
	print PFILE "-- max number of mismatches for near-perfect repeats: $maxmismatch\n";
	print PFILE "-- limit number of sequence paths to: $dfs_limit\n";
	print PFILE "-- output directory: $WORK\n";
	print PFILE "-- max number of parallel jobs: $MAX_PROCESSES\n";
	print PFILE "-- BAM file normal: $BAMNORMAL\n";
	print PFILE "-- BAM file tumor: $BAMTUMOR\n";
	print PFILE "-- bedfile: $BEDFILE\n";
	print PFILE "-- reference file: $REF\n";
	print PFILE "-- file of selected locations: $selected\n";
	print PFILE "-- output format for variants: $outformat\n";
	print PFILE "-- random seed: $rand_seed\n";
	if($intarget) { print PFILE "-- output variants in target? yes\n"; }
	else { print PFILE "-- output variants in target? no\n"; }
	if($logs) { print PFILE "-- keep log files? yes\n"; }
	else { print PFILE "-- keep log files? no\n"; }
	if($STcalling) { print PFILE "-- call variant in normal with samtools? yes\n"; }
	else { print PFILE "-- call variant in normal with samtools? no\n"; }
	
	close PFILE;
}

## link input files
#####################################################

sub linkFiles {
	
	print STDERR "-- Linking input bamfile\n";
	
	my $bamfile = $_[0]; 
	my $dirname = $_[1];
	
	my $flag = 0;
	
	print STDERR "$bamfile\n";
	my $fdir = "$WORK/$dirname";
	
	# get absolute path
	my $bam_abs_path = File::Spec->rel2abs($bamfile);
	
	if (! -r $fdir) {
		$flag = 1;
		## link the input files
		mkdir $fdir;
	
		if (-e "$bam_abs_path") {
			symlink "$bam_abs_path", "$fdir/bamfile.bam";
		}
		else {
			print STDERR "Can't find BAM file: $bamfile\n";
			next;
		}
		
		if (-e "$bam_abs_path.bai") {
			symlink "$bam_abs_path.bai", "$fdir/bamfile.bam.bai";
		}
		else {
			print STDERR "Indexing BAM file...\n";
			runCmd("index bamfile", "$bamtools index -in $bam_abs_path");
			symlink "$bam_abs_path.bai", "$fdir/bamfile.bam.bai";
			#print STDERR "Can't find index BAM file: $bamfile.bai\n";
			next;
		}
	}
	
	return $flag;
}

## load mutations and coverage info
#####################################################

sub loadHashes {
	
	print STDERR "-- Loading mutations and coverage info\n";
		
	# clear hashes 
	for my $k (keys %normalSVs) { delete $normalSVs{$k}; }
	for my $k (keys %tumorSVs) { delete $tumorSVs{$k}; }
	for my $k (keys %normalCov) { delete $normalCov{$k}; }
	for my $k (keys %tumorCov) { delete $tumorCov{$k}; }
	
	# load databases of mutations
	if(-e "$WORK/normal/variants.db.dir") {
		print STDERR "load variants for normal...\n";
		loadDB("$WORK/normal/variants.db", \%normalSVs, \%exons, 0);
	}
	if(-e "$WORK/tumor/variants.db.dir") {
		print STDERR "load variants for tumor...\n";
		loadDB("$WORK/tumor/variants.db", \%tumorSVs, \%exons, 0);
	}
	if(-e "$WORK/normal/refcov.txt.gz") {
		print STDERR "load coverage for normal...\n";
		loadCov("$WORK/normal/refcov.txt", \%normalCov);
	}
	if(-e "$WORK/tumor/refcov.txt.gz") {
		print STDERR "load coverage for tumor...\n";
		loadCov("$WORK/tumor/refcov.txt", \%tumorCov);
	}
	if(-e "$WORK/normal/loc2keys.txt.gz") {
		print STDERR "load loc2key for for normal...\n";
		loadLoc2Key("$WORK/normal/loc2keys.txt.gz", \%normalL2K);
	}
	if(-e "$WORK/tumor/loc2keys.txt.gz") {
		print STDERR "load loc2key for for tumor...\n";
		loadLoc2Key("$WORK/tumor/loc2keys.txt.gz", \%tumorL2K);
	}
	
	#remove reference coverage files:
	if(-e "$WORK/normal/refcov.txt.gz") {
		print STDERR "remove coverage file for normal...\n";
		runCmd("remove coverage file", "rm $WORK/normal/refcov.txt.gz");
	}
	if(-e "$WORK/tumor/refcov.txt.gz") {
		print STDERR "remove coverage file for tumor...\n";
		runCmd("remove coverage file", "rm $WORK/tumor/refcov.txt.gz");
	}
}

## call variants on each family memeber
#####################################################

sub callSVs {
	
	my $bamfile = $_[0]; 
	my $dir = $_[1]; 
	
	print STDERR "-- Detect mutations on $dir\n";

	my $bam_abs_path = File::Spec->rel2abs($bamfile);

	my $command = "$findVariants ".
		"--bam $bam_abs_path ".
		"--bed $BEDFILE ".
		"--ref $REF ".
		"--kmer $kmer ". 
		"--mincov $min_cov ". 
		"--covthr $cov_threshold ". 
		"--lowcov $tip_cov_threshold ".
		"--covratio $covratio ".
		"--radius $radius ".
		"--window $windowSize ".
		"--maxregcov $max_reg_cov ".
		"--step $delta ".
		"--mapscore $map_qual ".
		"--mismatches $maxmismatch ".
		"--pathlimit $dfs_limit ".
		"--dir $WORK/$dir ".
		"--sample ALL ".
		"--numprocs $MAX_PROCESSES ".
		"--coords $selected ".
		"--format $outformat ".
		"--seed $rand_seed ".		
		"--cov2file";
	if($intarget) { $command .= " --intarget"; }	
	if($logs) { $command .= " --logs"; }	
	if($VERBOSE) { $command .= " --verbose"; }

	print STDERR "Command: $command\n" if($VERBOSE);
	my $retval = runCmd("findVariants", "$command");
	
	return $retval;
}

## Compute best state for the location of the mutation
#####################################################
sub computeBestState {
	
	print STDERR "-- Compute best state for each variant\n";
	
	my $vars;	
	my $covthr = 0; # min cov threshold to remove sequencing errors
	
	foreach my $ID (@members) {
		
		if($ID eq "normal") { $vars = \%normalSVs; }
		if($ID eq "tumor") { $vars = \%tumorSVs; }
		
		foreach my $currKey (keys %$vars) {
			#print STDERR "$currKey\n";
			
			next if($currKey eq "numPaths");
	
			my $var = $vars->{$currKey};
			my $chr = $var->{chr};
			my $pos = $var->{pos}; 
			my $loc = "$chr:$pos";
	
			my $refState;
			#my $mutState;
			my $altState;
			
			my $refCov;
			my $mutCov;
			my $altCov;
						
			foreach my $role (@members) {
				#print STDERR "$role: ";
				
				my $varh;	
				my $refhash;
				my $l2k;
				
				if($role eq "normal") { $varh = \%normalSVs; $refhash = \%normalCov; $l2k = \%normalL2K; }
				if($role eq "tumor") { $varh = \%tumorSVs; $refhash = \%tumorCov; $l2k = \%tumorL2K; }
								
				$refState->{$role} = 0; $altState->{$role} = 0;
				$refCov->{$role} = 0; $mutCov->{$role} = 0; $altCov->{$role} = 0;
				
				my $sum = 0;
				my $cnt = 0;
				if ( exists($refhash->{$loc}) ) { # if there is reference coverage at location
					## compute reference coverage at locus
					for (my $t = $pos; $t<($pos+$kmer+1); $t++) { # compute avg coverage over window of size K
						my $loc2 = "$chr:$t";
						if ( exists $refhash->{$loc2} && 
							!( exists($l2k->{$loc2}) && ( (scalar @{$l2k->{$loc2}->{keys}})!=0) ) ) {
							$sum += $refhash->{$loc2};
							$cnt++;
						}					
					}
				}
				my $RCov = $sum;
				if($cnt>0) { $RCov = floor($sum/$cnt); }
				if($RCov > $covthr) { $refState->{$role} = 1; $refCov->{$role} = $RCov; }
				else { $refState->{$role} = 0; $refCov->{$role} = 0; }
		
				## compute coverage of mut at locus
				my $ACov = 0;
				if( exists($varh->{$currKey}) ) {
					my $altMut = $varh->{$currKey};
					$ACov = $altMut->{mincov};
				}
				if($ACov > $covthr) { $altState->{$role} = 1; $mutCov->{$role} = $ACov; }
				else { $altState->{$role} = 0; $mutCov->{$role} = 0; }
				
				## compute alternative allele coverage
				my $OCov = 0;
				if( exists($l2k->{$loc}) ) {
					foreach my $altKey (@{$l2k->{$loc}->{keys}}) {
						if( exists($varh->{$altKey}) && ($altKey ne $currKey) ) {
							my $altMut = $varh->{$altKey};
							$OCov += $altMut->{mincov};
						}
					}
				}
				if($OCov > $covthr) { $altCov->{$role} = $OCov; }
				else { $altCov->{$role} = 0; }
			}
			
			# 2 2/0 1 (ref: N,T / qry: N,T)
			my @RS=(); my @AS=();
			my @RC=(); my @AC=(); my @OC=();
			foreach my $role (@members) { # the order or members must be N,T
				if( ($refState->{$role} == 1) && ($altState->{$role} == 1) ) {
					push @RS , 1;
					push @AS , 1;
				}
				elsif( ($refState->{$role} == 1) && ($altState->{$role} == 0) ) {
					push @RS , 2;
					push @AS , 0;
				}
				elsif( ($refState->{$role} == 0) && ($altState->{$role} == 1) ) {
					push @RS , 0;
					push @AS , 2;
				}
				elsif( ($refState->{$role} == 0) && ($altState->{$role} == 0) ) {
					push @RS , 0;
					push @AS , 0;
				}
				
				push @RC , $refCov->{$role};
				push @AC , $mutCov->{$role};
				push @OC , $altCov->{$role};
			}
			
			my $stateStr = "$RS[0]$RS[1]/$AS[0]$AS[1]";
			my $covStr   = "$RC[0] $RC[1]/$AC[0] $AC[1]/$OC[0] $OC[1]";
			
			$var->{bestState} = $stateStr;
			$var->{covState} = $covStr;
			
			# compute hom/het state and chi2score
			parseCovState($var, $ID);
			
			$vars->{$currKey} = $var;
		}
	}
}

## parse covState and compute hom/het state and chi2score
##########################################

sub parseCovState {

	my $sv = $_[0];
	my $role = $_[1];
	my $BS = $sv->{covState};
	
	my ($REF,$ALT,$OTH) = split("/",$BS);
	my @R = split(" ", $REF); # reference
	my @A = split(" ", $ALT); # alternative
	my @O = split(" ", $OTH); # other

	my $id = 0;	
	if($role eq "normal") { $id = 0; }
	if($role eq "tumor") { $id = 1; }	
	
	my $alleleCnt = 0;	
	if($R[$id] > 0) { $alleleCnt++; }
	if($A[$id] > 0) { $alleleCnt++; }
	if($O[$id] > 0) { $alleleCnt++; }
	
	my $totCov = $R[$id] + $A[$id] + $O[$id];
	$sv->{altcov} = $totCov - $sv->{mincov};

	if($alleleCnt > 1) { $sv->{zygosity} = "het"; }
	else { $sv->{zygosity} = "hom"; }
	
	my $covRatio = 0;
	my $chi2Score = "na";
	if ($totCov > 0) { 
		$covRatio = sprintf('%.2f', ($sv->{mincov} / $totCov) );
		$sv->{covratio} = $covRatio; 

		#Chi-squared Test for allele coverage
		my $o1 = $sv->{mincov}; # observed 1
		my $o2 = $sv->{altcov}; # observed 2
		my $e1 = $totCov/2; # expected 1
		my $e2 = $totCov/2; # expected 2

		my $term1 = (($o1-$e1)*($o1-$e1))/$e1;
		my $term2 = (($o2-$e2)*($o2-$e2))/$e2;

		$chi2Score = sprintf('%.2f', $term1 + $term2);		
	}
	$sv->{chi2score} = $chi2Score;	
}

## parse bestState and update mutation info
##########################################
sub parseBestState {

	my $sv = $_[0];
	my $BS = $sv->{bestState};
	
	$sv->{somatic} = "no";
	$sv->{inheritance} = "no";
	
	if( ($BS eq "12/10") || ($BS eq "02/20") ) { # only in normal
		$sv->{inheritance} = "no";
	}			
	else {
		my ($REF,$ALT) = split("/",$BS);
		my @R = split("", $REF); # reference
		my @A = split("", $ALT); # alternative
		my $n=0; my $t=1;
		
		if ( ($A[$n]>0) && ($A[$t]>0) ) { # both in normal and tumor
			if ($sv->{covratio} <= $maxNormalContamination) { # partially in normal (contamination)
				$sv->{inheritance} = "normal";
				$sv->{somatic} = "yes";
			}
			else {
				$sv->{inheritance} = "normal";
			}
		}
		elsif ( ($A[$n]==0) && ($A[$t]>0) ) { # only in tumor
			$sv->{inheritance} = "no";
			$sv->{somatic} = "yes";
		}
		else { # unknown state
			$sv->{inheritance} = "no"; 
		}
		#if($A[$t]>0) {  $sv->{id} = "tumor"; }
	}	
}

## parse covState and returns total
## coverage at mutation locus
##########################################
sub getTotalCov {

	my $sv = $_[0];
	my $role = $_[1];
	my $BS = $sv->{covState};
	
	my ($REF,$ALT,$OTH) = split("/",$BS);
	my @R = split(" ", $REF); # reference
	my @A = split(" ", $ALT); # alternative
	my @O = split(" ", $OTH); # other

	my $id = 0;	
	if($role eq "normal") { $id = 0; }
	if($role eq "tumor") { $id = 1; }	
	
	my $totCov = $R[$id] + $A[$id] + $O[$id];
	
	return $totCov;
}


## find somatic events
#####################################################

sub findSomaticMut {
	
	print STDERR "-- Finding somatic mutations\n";
	
	my $family = $_[0];
	
	if (-e "$WORK/somatic.db.dir") { runCmd("remove database files", "rm $WORK/somatic.db.*"); }
	if (-e "$WORK/common.db.dir") { runCmd("remove database files", "rm $WORK/common.db.*"); }
	
	my $somatic_dbm_obj = tie %somaticSVs, 'MLDBM::Sync', "$WORK/somatic.db", O_CREAT|O_RDWR, 0640; 
	my $common_dbm_obj = tie %commonSVs, 'MLDBM::Sync', "$WORK/common.db", O_CREAT|O_RDWR, 0640; 
	
	$somatic_dbm_obj->Lock;
	$common_dbm_obj->Lock;
	
	my $stats;
	$stats->{tumor_name} = $tumor_name;
	$stats->{normal_name} = $normal_name;
	
	$somaticSVs{stats} = $stats;
	$commonSVs{stats} = $stats;
		
	foreach my $k (keys %tumorSVs) {
				
		my $mut = $tumorSVs{$k};
		
		# parse bestState and update info
		parseBestState($mut);
				
		my $chr = $mut->{chr};
		my $pos = $mut->{pos};
		my $key = "$chr:$pos";
		
		#print STDERR "$k\t$mut->{bestState}\t$mut->{covState}\t$mut->{somatic}\t$normalCov{$key}\n";
		
		# skip mutation if normal not sampled
		my $totcov = $mut->{altcov} + $mut->{mincov};
		next if ($totcov < $min_cov);
		
		# skip mutation if reference for normal not well sampled
		next if (getTotalCov($mut,"normal") < $min_cov);
		#next if (!exists $normalCov{$key});
		#next if ($normalCov{$key} < $min_cov); # min cov requirement
		
		# skip mutation if present in normal mutations from read alignments
		#next if(exists $alnHashN{$k});		
		#if(exists $alnHashN{$k}) {
			#next;
			#my $aln_mut = $alnHashN{$k};
			#if ($aln_mut->{imf} > $maxNormalContamination) { next; }	
		#}
		
		# compute inheritance
		if($mut->{somatic} eq "yes") {
			$somaticSVs{$k} = $mut;
		}
		else { # non somatic
			if ($mut->{inheritance} eq "normal") {
				$commonSVs{$k} = $mut;
			}
		}
		$tumorSVs{$k} = $mut;
	}

	$somatic_dbm_obj->UnLock;
	$common_dbm_obj->UnLock;
}


## export family mutations to file
#####################################################

sub exportSVs {
	
	print STDERR "-- Exporting variants to file\n";
		
	# export somatic variants
	
	my $command_somatic = "$exportTool --somatic ".
		"--db $WORK/somatic.db ".
		"--bed $BEDFILE ".
		"--ref $REF ".
		"--output-format $outformat ".
		"--variant-type indel ". 
		"--min-alt-count-tumor $min_cov ". 
		"--min-vaf-tumor $outratio";
	if($intarget) { $command_somatic .= " --intarget"; }
	if ($outformat eq "annovar") { $command_somatic .= " > $WORK/somatic.indel.annovar"; }
	elsif ($outformat eq "vcf") { $command_somatic .= " > $WORK/somatic.indel.vcf"; }
	
	print STDERR "Command: $command_somatic\n" if($VERBOSE);
	runCmd("export somatic", $command_somatic);
	
	# export common variants
	
	my $command_comm = "$exportTool --somatic ".
		"--db $WORK/common.db ".
		"--bed $BEDFILE ".
		"--ref $REF ".
		"--output-format $outformat ".
		"--variant-type indel ". 
		"--min-alt-count-tumor $min_cov ". 
		"--max-alt-count-normal 1000000 ".
		"--max-vaf-normal 1.0 ".
		"--min-vaf-tumor $outratio";
	if($intarget) { $command_comm .= " --intarget"; }
	if ($outformat eq "annovar") { $command_comm .= " > $WORK/common.indel.annovar"; }
	elsif ($outformat eq "vcf") { $command_comm .= " > $WORK/common.indel.vcf"; }
	
	print STDERR "Command: $command_comm\n" if($VERBOSE);
	runCmd("export common", $command_comm);
}

## process single family
#####################################################

sub processBAM {

	my $bamfile = $_[0];
	my $outdir = $_[1];
	
	#my $toprocess = linkFiles("$bamfile", "$outdir"); # link input files
	#if($toprocess == 1) {
		my $retval = callSVs($bamfile, $outdir); # call mutations on each family
	#}
	#else { print STDERR "$outdir already processed!\n"; }
	
	#print STDERR "retval = $retval\n";
	if($retval < 0) { exit; }
	
	if($outdir eq "tumor") {
		$tumor_name = extractSM($WORK,$outdir);
		#print "SMtumor = $tumor_name\n";	
	}
	if($outdir eq "normal") {
		$normal_name = extractSM($WORK,$outdir);
		#print "SMnormal = $normal_name\n";	
	}
}

## call mutations from alignments
#####################################################

sub callMutFromAlignment {
	
	print STDERR "-- Calling mutations from alignments (samtools/bcftools)\n";
	
	my $cnt = 0;
	if($selected ne "null") {
		$cnt = loadCoordinates("$selected", \%exons, $VERBOSE);
	}
	else {
		$cnt = loadRegions("$BEDFILE", \%exons, $radius, $VERBOSE);		
	}
	# split region in $MAX_PROCESSES of eqaul size
	my $step = $cnt;
	if ($MAX_PROCESSES>0) {
		$step = ceil($cnt/$MAX_PROCESSES);
	}
	print STDERR "stepSize: $step\n";
	
	my $num_partititons=1; 
	my $count = 0;
	foreach my $k (keys %exons) { # for each chromosome
		foreach my $exon (@{$exons{$k}}) { # for each exon
			if($count>=$step) { 
				$num_partititons++; 
				push @{$partitions{$num_partititons}}, $exon; 
				$count = 1; # reset counter
			}
			else { 
				push @{$partitions{$num_partititons}}, $exon; 
				$count++;
			}
		}
	}
	
	# run parallel jobs
	my $pm = new Parallel::ForkManager($MAX_PROCESSES);
	
	# Setup a callback for when a child finishes up so we can get it's exit code
	$pm->run_on_finish( sub {
	    my ($pid, $exit_code, $ident, $exit_signal, $core_dump) = @_;
	    print STDERR "** $ident just got out of the pool ".
		"[PID: $pid | exit code: $exit_code | exit signal: $exit_signal | core dump: $core_dump]\n";
	});

	# Setup a callback for when a child starts
	$pm->run_on_start( sub {
	    my ($pid, $ident) = @_;
	    print STDERR "** $ident started, pid: $pid\n";
	});
	
	for(my $i = 1; $i <= $num_partititons; $i++) {
	
		my $pid = $pm->start($i) and next;

		my $size = scalar(@{$partitions{$i}});
		print STDERR "Alignment analysis $i [$size].\n";
		
		alignmentAnalysis($i);
		
		$pm->finish($i); # Terminates the child process
	}
	$pm->wait_all_children; # wait for all parallel jobs to complete
	

	# load genome from fasta file (if needed) and pass it to loadVCF()
	if($USEFAIDX == 0) {
		loadGenomeFasta($REF, \%genome);
	}
	
	# merge results
	for(my $i = 1; $i <= $num_partititons; $i++) {
		
		my $outvcfnormal = "$WORK/normal/aln.$i.vcf";
		#my $outvcftumor = "$WORK/tumor/aln.$i.vcf";
		
		loadVCF($outvcfnormal, \%alnHashN, $REF, $USEFAIDX, \%genome);
		#loadVCF($outvcftumor, \%alnHashT, $REF);
		
		runCmd("remove align file (normal)", "rm $outvcfnormal");
		#runCmd("remove align file (tumor)", "rm $outvcftumor");
	}
}

## update table of normal variants based on aligment analysis
## (only for loci with variant in the tumor
#####################################################

sub updateNormalVariants {

	print STDERR "-- Update supporting coverage in normal\n";

	foreach my $key (keys %tumorSVs) {
		next if($key eq "numPaths");
		
		if( (exists $alnHashN{$key}) && (!exists $normalSVs{$key}) ) { # add variant in normal
			$normalSVs{$key} = $alnHashN{$key};
		}
	}
}

## analize alignments
#####################################################
sub alignmentAnalysis {
	
	my $p = $_[0];	
	
	my $commandN = "";
	#my $commandT = "";
	
	my $outvcfnormal = "$WORK/normal/aln.$p.vcf";
	#my $outvcftumor = "$WORK/tumor/aln.$p.vcf";
	
	open FN, "> $outvcfnormal" or die "Can't open $outvcfnormal ($!)\n"; print FN "";
	#open FT, "> $outvcftumor" or die "Can't open $outvcftumor ($!)\n"; print FT "";
		
	foreach my $exon (@{$partitions{$p}}) { # for each exon			
			
		my $chr = $exon->{chr};
		my $start = $exon->{start};
		my $end = $exon->{end};
		
		## extend region left and right by radius
		my $left = $start-$radius;
		if ($left < 0) { $left = $start; }
		my $right = $end+$radius;
		
		my $REG = "$chr:$left-$right";
		#print STDERR "$REG\n";
			 
		$commandN = "$samtools mpileup -Q 0 -F 0.0 -r $REG -uf $REF $BAMNORMAL | $bcftools call -p 1.0 -P 0 -V snps -Am - | awk '\$0!~/^#/' >> $outvcfnormal";
		#$commandT = "$samtools mpileup -Q 0 -F 0.0 -r $REG -uf $REF $BAMTUMOR | $bcftools call -p 1.0 -P 0 -V snps -Am - | awk '\$0!~/^#/' >> $outvcftumor";
			
		runCmd("samtools/bcftools calling", "$commandN");
		#runCmd("samtools/bcftools calling", "$commandT");
	}
	close(FN);
	#close(FT);
}

## do the job
##########################################

init();
printParams("$WORK/parameters.txt"); # print parameters
processBAM($BAMNORMAL, "normal"); # process normal
processBAM($BAMTUMOR, "tumor"); # process tumor
loadHashes(); # load mutations and coverage info
if($STcalling) { 
	callMutFromAlignment(); 
	updateNormalVariants();
} # call mutations from alignment using samtools
computeBestState(); # compute bestState for each mutation
findSomaticMut(); # detect somatic mutations in tumor
exportSVs(); # export mutations to file

##########################################


my $time_taken = time - $start_time;

#if($VERBOSE) {
	elapsedTime($time_taken, "FindSomatic");
#}
