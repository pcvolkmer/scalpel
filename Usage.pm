package Usage;

###################################################################
# Usage
#
# Package with usage description of each operation mode 
# 
#  Author: Giuseppe Narzisi 
#    Date: December 11, 2013
#
###################################################################

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(getDefaults usageDenovo usageSingle usageSomatic usageExportSingle usageExportDenovo usageExportSomatic);

use strict;
use warnings;
use POSIX;

my %defaults = (
	kmer 				=> 25,
	cov_threshold 		=> 5,
	tip_cov_threshold	=> 2,
	covratio			=> 0.01,
	radius	 			=> 100,
	windowSize 			=> 400,
	delta				=> 100,
	map_qual			=> 1,
	maxmismatch 		=> 3,
	max_reg_cov			=> 10000,	

	outratio			=> 0.05,
	WORK 				=> "./outdir",
	MAX_PROCESSES 		=> 1,
	sample		 		=> "ALL",
	selected 			=> "null",
	pathlimit 			=> 1000000,
	format				=> "vcf",
	SVtype				=> "indel",
	version_num			=> "0.5.1 (beta), July 31 2015",
	
	minInsSize			=> 1,
	maxInsSize			=> 1000000000,
	minDelSize			=> 1,
	maxDelSize			=> 1000000000,
	
	# somatic parameters
	minAltCntTumor		=> 4, 
	maxAltCntNormal		=> 0,
	minVafTumor			=> 0.05,
	maxVafNormal		=> 0.0,
	minCovTumor			=> 4,
	maxCovTumor			=> 1000000000,
	minCovNormal		=> 10,
	maxCovNormal		=> 1000000000,
	minPhredFisher		=> 10,
	#maxPhredFisher		=> 1000000000,
	
	# denovo parameters
	minAltCntAffected 	=> 5,
	maxAltCntUnaffected => 0, 
	minVafAffected 		=> 0.05,
	maxVafUnaffected 	=> 0.0,
	minCovAffected 		=> 5, 
	maxCovAffected 		=> 1000000000,
	minCovUnaffected 	=> 10,
	maxCovUnaffected 	=> 1000000000,
	minchi2				=> 0,
	maxchi2				=> 20,
	
	# single paramters
	min_cov 			=> 5,
	max_cov				=> 1000000,
);

#####################################################
#
# export defaults
#
sub getDefaults { return \%defaults; }


#####################################################
#
# Message about the program and how to use it
#
sub usageSingle {	

my $name = $_[0];
print STDERR <<END;

usage: $name --bam <BAM file> --bed <BED file> --ref <FASTA file> [OPTIONS]
	
Detect indels in one single dataset (e.g., one individual).
	
OPTIONS:

    --help             : this (help) message
    --verbose          : verbose mode

  Required:
    --bam <BAM file>   : BAM file with the reference-aligned reads
    --bed <BED file>   : file with list of regions (BED format) in sorted order or single region in format chr:start-end (example: 1:31656613-31656883)		
    --ref <FASTA file> : reference genome in FASTA format (same one that was used to create the BAM file)

  Optional:
    --kmer <int>       : k-mer size [default $defaults{kmer}]
    --covthr <int>     : threshold used to select source and sink [default $defaults{cov_threshold}]
    --lowcov <int>     : threshold used to remove low-coverage nodes [default $defaults{tip_cov_threshold}]
    --covratio <float> : minimum coverage ratio for sequencing errors (default: $defaults{covratio})
    --radius <int>     : left and right extension (in base-pairs) [default $defaults{radius}]
    --window <int>     : window-size of the region to assemble (in base-pairs) [default $defaults{windowSize}]
    --maxregcov <int>  : maximum average coverage allowed per region [default $defaults{max_reg_cov}]
    --step <int>       : delta shift for the sliding window (in base-pairs) [default $defaults{delta}]
    --mapscore <int>   : minimum mapping quality for selecting reads to assemble [default $defaults{map_qual}]
    --pathlimit <int>  : limit number of sequence paths to [default $defaults{pathlimit}]
    --mismatches <int> : max number of mismatches in near-perfect repeat detection [default $defaults{maxmismatch}]
    --dir <directory>  : output directory [default $defaults{WORK}]
    --numprocs <int>   : number of parallel jobs (1 for no parallelization) [default $defaults{MAX_PROCESSES}]
    --sample <string>  : only process reads/fragments in sample [default $defaults{sample}]
    --coords <file>    : file with list of selected locations to examine [default $defaults{selected}]

  Output:
    --format           : export mutations in selected format (annovar | vcf) [default $defaults{format}]
    --intarget         : export mutations only inside the target regions from the BED file
    --logs             : keep log files

  Note 1: the list of detected indels is saved in file: OUTDIR/variants.*.indel.*
  where OUTDIR is the output directory selected with option "--dir" [default $defaults{WORK}]

  Note 2: use the export tool (scalpel-export) to export mutations using different filtering criteria

END
exit;
}

#####################################################
#
# Message about the program and how to use it
#
sub usageDenovo {
	
my $name = $_[0];
print STDERR <<END;

usage: $name --dad <BAM file> --mom <BAM file> --aff <BAM file> --sib <BAM file> --bed <BED file> --ref <FASTA file> [OPTIONS]

Detect de novo indels in a family of four individuals (mom, dad, aff, sib).

OPTIONS:

    --help             : this (help) message
    --verbose          : verbose mode

  Required:
    --dad <BAM file>   : father BAM file
    --mom <BAM file>   : mother BAM file
    --aff <BAM file>   : affected child BAM file
    --sib <BAM file>   : sibling BAM file
    --bed <BED file>   : file with list of regions (BED format) in sorted order or single region in format chr:start-end (example: 1:31656613-31656883)		
    --ref <FASTA file> : reference genome in FASTA format (same one that was used to create the BAM file)

  Optional:
    --kmer <int>       : k-mer size [default $defaults{kmer}]
    --covthr <int>     : threshold used to select source and sink [default $defaults{cov_threshold}]
    --lowcov <int>     : threshold used to remove low-coverage nodes [default $defaults{tip_cov_threshold}]
    --covratio <float> : minimum coverage ratio for sequencing errors (default: $defaults{covratio})
    --radius <int>     : left and right extension (in base-pairs) [default $defaults{radius}]
    --window <int>     : window-size of the region to assemble (in base-pairs) [default $defaults{windowSize}]
    --maxregcov <int>  : maximum average coverage allowed per region [default $defaults{max_reg_cov}]
    --step <int>       : delta shift for the sliding window (in base-pairs) [default $defaults{delta}]
    --mapscore <int>   : minimum mapping quality for selecting reads to assemble [default $defaults{map_qual}]
    --pathlimit <int>  : limit number of sequence paths to [default $defaults{pathlimit}]
    --mismatches <int> : max number of mismatches in near-perfect repeat detection [default $defaults{maxmismatch}]
    --dir <directory>  : output directory [default $defaults{WORK}]
    --numprocs <int>   : number of parallel jobs (1 for no parallelization) [default $defaults{MAX_PROCESSES}]
    --coords <file>    : file with list of selected coordinates to examine [default $defaults{selected}]
    --two-pass         : perform second pass of analysis to confirm candidate calls

  Output:
    --format           : export mutations in selected format (annovar | vcf) [default $defaults{format}]
    --intarget         : export mutations only inside the target regions from the BED file
    --logs             : keep log files

  Note 1: the list of de novo indels is saved in file: OUTDIR/denovos.*.indel.*
  where OUTDIR is the output directory selected with option "--dir" [default $defaults{WORK}]

  Note 2: use the export tool (scalpel-export) to export mutations using different filtering criteria

END
exit;
}

#####################################################
#
# Message about the program and how to use it
#
sub usageSomatic {

my $name = $_[0];
print STDERR <<END;

usage: $name --normal <BAM file> --tumor <BAM file> --bed <BED file> --ref <FASTA file> [OPTIONS]

Detect somatic indels in a tumor/normal pair

OPTIONS:

    --help                : this (help) message
    --verbose             : verbose mode

  Required:
    --normal <BAM file>   : normal BAM file
    --tumor  <BAM file>   : tumor BAM file
    --bed    <BED file>   : file with list of regions (BED format) in sorted order or single region in format chr:start-end (example: 1:31656613-31656883)	
    --ref    <FASTA file> : reference genome in FASTA format (same one that was used to create the BAM file)

  Optional:
    --kmer <int>          : k-mer size [default $defaults{kmer}]
    --covthr <int>        : threshold used to select source and sink [default $defaults{cov_threshold}]
    --lowcov <int>        : threshold used to remove low-coverage nodes [default $defaults{tip_cov_threshold}]
    --covratio <float>    : minimum coverage ratio for sequencing errors (default: $defaults{covratio})
    --radius <int>        : left and right extension (in base-pairs) [default $defaults{radius}]
    --window <int>        : window-size of the region to assemble (in base-pairs) [default $defaults{windowSize}]
    --maxregcov <int>     : maximum average coverage allowed per region [default $defaults{max_reg_cov}]
    --step <int>          : delta shift for the sliding window (in base-pairs) [default $defaults{delta}]
    --mapscore <int>      : minimum mapping quality for selecting reads to assemble [default $defaults{map_qual}]
    --pathlimit <int>     : limit number of sequence paths to [default $defaults{pathlimit}]
    --mismatches <int>    : max number of mismatches in near-perfect repeat detection [default $defaults{maxmismatch}]
    --dir <directory>     : output directory [default $defaults{WORK}]
    --numprocs <int>      : number of parallel jobs (1 for no parallelization) [default $defaults{MAX_PROCESSES}]
    --coords <file>       : file with list of selected coordinates to examine [default $defaults{selected}]
    --two-pass            : perform second pass of analysis to confirm candidate calls

  Output:
    --format              : export mutations in selected format (annovar | vcf) [default $defaults{format}]
    --intarget            : export mutations only inside the target regions from the BED file
    --logs                : keep log files
	
  Note 1: the list of somatic indels is saved in file: OUTDIR/somatic.*.indel.* 
  where OUTDIR is the output directory selected with option "--dir" [default $defaults{WORK}].
  
  Note 2: use the export tool (scalpel-export) to export mutations using different filtering criteria

END
exit;
}


#####################################################
#
# Message about this program and how to use it
#
sub usageExportSingle {
	
my $name = $_[0];
print STDERR <<END;

usage: $name --db <file> --bed <file> --ref <file> [OPTIONS]

OPTIONS:

    --help             : this (help) message
    --verbose          : verbose mode

  Required:
    --db <file>        : Database of mutations
    --bed <file>       : file with list of regions (BED format) in sorted order or single region in format chr:start-end (example: 1:31656613-31656883)
    --ref <FASTA file> : reference genome in FASTA format (same one that was used to create the BAM file)
  
  Optional:
    --output-format <text>   : output format for variants (annovar | vcf) [default $defaults{format}]
    --variant-type <text>    : mutation type (snp, del, ins, indel, all: everything) [default $defaults{SVtype}]
    --min-ins-size <int>     : minimum size of an insertion [default $defaults{minInsSize}]
    --max-ins-size <int>     : maximum size of an insertion [default $defaults{maxInsSize}]
    --min-del-size <int>     : minimum size of a deletion [default $defaults{minDelSize}]
    --max-del-size <int>     : maximum size of a deletion [default $defaults{maxDelSize}]
    --min-alt-count <int>    : minimum alternative count [default $defaults{min_cov}]
    --max-alt-count <int>    : maximum alternative count [default $defaults{max_cov}]
    --min-chi2-score <float> : minimum chi-square score [default $defaults{minchi2}]
    --max-chi2-score <float> : maximum chi-square score [default $defaults{maxchi2}]
    --min-vaf <float>        : minimum variant allele frequency (AlleleCov/TotCov) [default $defaults{outratio}]
    --intarget               : export mutations only inside the target regions from the BED file

  Supported output formats:
    1. annovar
    2. vcf
	
    NOTE: The database.db file can be found in the output directory for the single operation 
    mode or in the correspective subdirectories ("main" and "twopass' for denovo and soamtic modes).

END
exit;
}

#####################################################
#
# Message about this program and how to use it
#
sub usageExportSomatic {
	
my $name = $_[0];
print STDERR <<END;

usage: $name --db <file> --bed <file> --ref <file> [OPTIONS]

OPTIONS:

    --help             : this (help) message
    --verbose          : verbose mode

  Required:
    --db <file>        : Database of mutations
    --bed <file>       : file with list of regions (BED format) in sorted order or single region in format chr:start-end (example: 1:31656613-31656883)
    --ref <FASTA file> : reference genome in FASTA format (same one that was used to create the BAM file)
  
  Optional:
    --output-format <text>       : output format for variants (annovar | vcf) [default $defaults{format}]
    --variant-type <text>        : mutation type (snp, del, ins, indel, all: everything) [default $defaults{SVtype}]
    --min-ins-size <int>         : minimum size of an insertion [default $defaults{minInsSize}]
    --max-ins-size <int>         : maximum size of an insertion [default $defaults{maxInsSize}]
    --min-del-size <int>         : minimum size of a deletion [default $defaults{minDelSize}]
    --max-del-size <int>         : maximum size of a deletion [default $defaults{maxDelSize}]
    --min-alt-count-tumor <int>  : minimum alternative count in the tumor [default $defaults{minAltCntTumor}]
    --max-alt-count-normal <int> : maximum alternative count in the normal [$defaults{maxAltCntNormal}]
    --min-vaf-tumor <float>      : minimum variant allele frequency (AlleleCov/TotCov) in the tumor [default $defaults{minVafTumor}]
    --max-vaf-normal <float>     : maximum variant allele frequency (AlleleCov/TotCov) in the normal [default $defaults{maxVafNormal}]
    --min-coverage-tumor <int>   : minimum coverage in the tumor [default $defaults{minCovTumor}]
    --max-coverage-tumor <int>   : maximum coverage in the tumor [default $defaults{maxCovTumor}]
    --min-coverage-normal <int>  : minimum coverage in the normal [default $defaults{minCovNormal}]
    --max-coverage-normal <int>  : maximum coverage in the normal [default $defaults{maxCovNormal}]
    --min-phred-fisher <float>   : minimum fisher exact test score [default $defaults{minPhredFisher}]
    --intarget                   : export mutations only inside the target regions from the BED file

  Supported output formats:
    1. annovar
    2. vcf
	
    NOTE: The database.db file can be found in the output directory for the single operation 
    mode or in the correspective subdirectories ("main" and "twopass' for denovo and soamtic modes).

END
exit;
}

#####################################################
#
# Message about this program and how to use it
#
sub usageExportDenovo {
	
my $name = $_[0];
print STDERR <<END;

usage: $name --db <file> --bed <file> --ref <file> [OPTIONS]

OPTIONS:

    --help             : this (help) message
    --verbose          : verbose mode

  Required:
    --db <file>        : Database of mutations
    --bed <file>       : file with list of regions (BED format) in sorted order or single region in format chr:start-end (example: 1:31656613-31656883)
    --ref <FASTA file> : reference genome in FASTA format (same one that was used to create the BAM file)
  
  Optional:
    --output-format <text>           : output format for variants (annovar | vcf) [default $defaults{format}]
    --variant-type <text>            : mutation type (snp, del, ins, indel, all: everything) [default $defaults{SVtype}]
    --min-ins-size <int>             : minimum size of an insertion [default $defaults{minInsSize}]
    --max-ins-size <int>             : maximum size of an insertion [default $defaults{maxInsSize}]
    --min-del-size <int>             : minimum size of a deletion [default $defaults{minDelSize}]
    --max-del-size <int>             : maximum size of a deletion [default $defaults{maxDelSize}]
    --min-alt-count-affected <int>   : minimum alternative count in the affected sample [default $defaults{minAltCntAffected}]
    --max-alt-count-unaffected <int> : maximum alternative count in the unaffected samples [default $defaults{maxAltCntUnaffected}]
    --min-vaf-affected <float>       : minimum variant allele frequency (AlleleCov/TotCov) in the affected sample [default $defaults{minVafAffected}]
    --max-vaf-unaffected <float>     : maximum variant allele frequency (AlleleCov/TotCov) in the unaffected samples [default $defaults{maxVafUnaffected}]
    --min-coverage-affected <int>    : minimum coverage in the affected sample [default $defaults{minCovAffected}]
    --max-coverage-affected <int>    : maximum coverage in the affected sample [default $defaults{maxCovAffected}]
    --min-coverage-unaffected <int>  : minimum coverage in the unaffected samples [default $defaults{minCovUnaffected}]
    --max-coverage-unaffected <int>  : maximum coverage in the unaffected samples [default $defaults{maxCovUnaffected}]
    --min-chi2-score <float>         : minimum chi-square score [default $defaults{minchi2}]
    --max-chi2-score <float>         : maximum chi-square score [default $defaults{maxchi2}]
    --intarget                       : export mutations only inside the target regions from the BED file

  Supported output formats:
    1. annovar
    2. vcf
	
    NOTE: The database.db file can be found in the output directory for the single operation 
    mode or i