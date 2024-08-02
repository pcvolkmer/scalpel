package SequenceIO;

###################################################################
# SequenceIO
#
# Package for common IO routines
#
#  Author: Giuseppe Narzisi
#    Date: December 11, 2013
#
###################################################################

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(loadCoordinates loadRegions loadVCF loadExonsBed loadGenomeFasta parseHeader saveReadGroup extractCoords extractCoordsDB extractSM);

use strict;
use warnings;
use FindBin qw($Bin);
use lib $Bin; # add $Bin directory to @INC
use Utils;
use HashesIO;
use Digest::MD5 qw(md5_hex);

my $bamtools = "$Bin/bamtools/bin/bamtools";

## extract coordinates from file
#####################################################
sub extractCoords {

	my $in = $_[0];
	my $out = $_[1];

	print STDERR "Extract coordinates of candidates...";

	open IN, "< $in" or die "Can't open $in ($!)\n";
	open OUT, "> $out" or die "Can't open $out ($!)\n";

	my $cnt = 0;
	while (<IN>) {
		chomp;
		next if ($_ =~ /^#/); # skip comments

		my @A = split /\t/, $_;
		print OUT "$A[0]\t$A[1]\n";
		$cnt++
	}

	close IN;
	close OUT;

	return $cnt;
}

## extract coordinates from DB
#####################################################
sub extractCoordsDB {

	my $dbfile = $_[0];
	my $bedfile = $_[1];
	my $intarget = $_[2];
	my $out = $_[3];

	my %exons;
	my %variants;
	my %dups; # used to remove duplicated locations

	print STDERR "Extract coordinates of candidates...";

	loadRegions("$bedfile", \%exons, 0, 0);
	loadDB("$dbfile", \%variants, \%exons, $intarget);

	open OUT, "> $out" or die "Can't open $out ($!)\n";

	my $cnt = 0;
	foreach my $key (keys %variants) {
		my $mut = $variants{$key};

		next if($mut->{type} eq "snp"); # skip SNPs

		my $chr = $mut->{chr};
		my $pos = $mut->{pos};
		my $loc = "$chr\t$pos";

		if( !(exists $dups{$loc}) ) {
			print OUT "$loc\n";
			$cnt++;
			$dups{$loc} = 1;
		}
	}

	close OUT;

	return $cnt;
}


## load selected location (coordinates)
#####################################################
sub loadCoordinates {

	my $file = $_[0];
	my $exons = $_[1];
	my $VERBOSE = $_[2];

	print STDERR "Loading coordinates...";

	my $cntall_locs = 0;
	if ($file ne "null") {

		open SELECTED, "< $file" or die "Can't open $file ($!)\n";

		while (<SELECTED>) {
			chomp;
			next if ($_ =~ /^#/); # skip comments

			my ($chr, $pos) = split /\t/, $_, 2;
			#if ($chr =~ /^chr/) { $chr = substr($chr,3); }

			my $exon;
			$exon->{chr} = $chr;
			$exon->{start} = $pos;
			$exon->{end} = $pos;

			push @{$exons->{$chr}}, $exon;
			$cntall_locs++;
		}
		close SELECTED;

		#if($VERBOSE) {
			print STDERR "$cntall_locs locations.\n";
		#}
	}
	return $cntall_locs;
}

## load regions to process
#####################################################
sub loadRegions {

	my $bedfile = $_[0];
	my $exons = $_[1];
	my $radius = $_[2];
	my $VERBOSE = $_[3];

	my $cnt = 0;

	if (-e "$bedfile") {
		$cnt = loadExonsBed("$bedfile", $exons, $radius, $VERBOSE);
	}
	elsif ($bedfile =~ m/(\w+):(\d+)-(\d+)/) { # try to parse as region (e.g., 1:1234-4567)

		print STDERR "Loading targets region: $bedfile\n";

		#parse region
		my ($chr,$interval) = split(':',$bedfile);
		my ($start,$end) = split('-',$interval);

		my $exon;
		$exon->{chr} = $chr;
		$exon->{start} = $start;
		$exon->{end} = $end;

		push @{$exons->{$chr}}, $exon;

		$cnt = 1;
	}
	else {
		print STDERR "Error: unrecognized string in --bed option: $bedfile\n";
		exit;
	}

	return $cnt;
}

## load exons list
#####################################################
sub loadExonsBed {

	print STDERR "Loading targets from BED file...";

	my $file = $_[0];
	my $exons = $_[1];
	my $radius = $_[2];
	my $VERBOSE = $_[3];

	open EXONSLIST, "< $file" or die "Can't open $file ($!)\n";

	my $cntall_exons = 0;
	my $cntovl_exons = 0;

	my $prev_exon;
	$prev_exon->{chr} = 0;
	$prev_exon->{start} = 0;
	$prev_exon->{end} = 0;
	while (<EXONSLIST>)
	{
		chomp;
		next if ($_ =~ /^#/ || $_ =~ /^@/); # skip comments

		#my ($chr, $start, $end) = split /\t/, $_, 3;
		my @array = split /\t/, $_;
		my $chr = $array[0];
		my $start = $array[1];
		my $end = $array[2];
		#if ($chr =~ /^chr/) { $chr = substr($chr,3); }

		# exons coordinate are sorted in the bed file
		if ( ($prev_exon->{chr} eq $chr) && ($start <= $prev_exon->{end}) ) {
			#print STDERR "exons regions overlap: [$prev_exon->{start},$prev_exon->{end}] [$start,$end]\n";
			$prev_exon->{end} = $end;
			$cntovl_exons++;
		}
		else {
			my $exon;
			$exon->{chr} = $chr;
			$exon->{start} = $start;
			$exon->{end} = $end;

			## extend region left and right by radius
			#my $l = $start-$radius;
			#my $u = $end+$radius;

			push @{$exons->{$exon->{chr}}}, $exon;

			# update prev exon
			$prev_exon = $exon;

			$cntall_exons++;
		}
		#last if($cntall_exons > 10);
	}

	if($VERBOSE) {
		# print exons
		foreach my $k (keys %$exons) { # for each chromosome
			foreach my $exon (@{$exons->{$k}}) { # for each exon
				print STDERR "$exon->{chr}\t$exon->{start}\t$exon->{end}\n";
			}
		}
	}

	close EXONSLIST;

	#if($VERBOSE) {
		print STDERR "$cntall_exons targets (filtered $cntovl_exons overlapping).\n";
	#}

	return $cntall_exons;
}

# load VCF file
##########################################
sub loadVCF {

    print STDERR "Loading mutations from VCF file...";

    my $file = $_[0];
	my $hash = $_[1];
	my $REF  = $_[2];
	my $FAIDX = $_[3];
	my $genome = $_[4];

    open VCF, "< $file" or die "Can't open $file ($!)\n";

	my $cnt = 0;
    my $header = "";
    while ( <VCF> ) {
        chomp;
        #print "$_\n";
        if($_ =~ m/^#/) {
			$header = $_;
            #print "$header\n";
        }
        else {
			#2       154110840       .       ATTTTTT ATTTTT  3.45209 .       INDEL;IDV=1;IMF=0.0909091;DP=11;VDB=0.135767;SGB=-0.379885;MQSB=0.763675;MQ0F=0;AC=0;AN=2;DP4=6,4,1,0;MQ=56     GT:PL   0/0:0,17,82
            my @fields = split("\t", $_);
			my $chr = $fields[0];
			my $pos = $fields[1];
			my $ref = $fields[3];
			my $alts = $fields[4];
			my $info = $fields[7];

			# extract supporting coverage fror the allele
			# IDV (=Maximum number of reads supporting an indel) and IMF (=Maximum fraction of reads supporting an indel)
			my @I = split(';', $info);

			my $idv = 0;
			my $imf = 0;
			#my $dp  = 0;
			foreach (@I) {
				my ($tag,$value) = split ('=',$_);
				if($tag eq "IDV") { $idv = $value; }
				if($tag eq "IMF") { $imf = $value; }
				#if($tag eq "DP")  { $dp = $value; }
			}

			# process each alternative allele
			my @A = split(',', $alts);

			foreach (@A) {
				#$cnt++;

				my $alt = $_;

				my $type;
				my $newref;
				my $newalt;
				my $prevbpalt;
				my $prevbpref;
				my $len;
				if(length($ref)==1 && length($alt)==1 ) {
					$type = "snp";
					$len = 1;
					$newref = $ref;
				 	$newalt = $alt;
				}
				elsif(length($ref) > length($alt)) {
					$type = "del";
					#my $idx = index($ref, $alt);
					#$newref = substr($ref,length($alt));
					#$len = length($newref);
					$len = length($ref)-length($alt);
					$newalt = ('-' x $len);
					$newref = substr($ref,1,$len);
					$prevbpref = substr($ref,0,1);
					$prevbpalt = $prevbpref;
					$pos+=1;
				}
				elsif(length($ref) < length($alt)) {
					$type = "ins";
					#my $idx = index($alt, $ref, 1);
					#$newalt = substr($alt,length($ref));
					#$len = length($newalt);
					$len = length($alt)-length($ref);
					$newalt = substr($alt,1,$len);
					$prevbpalt = substr($alt,0,1);
					$prevbpref = $prevbpalt;
					$newref = ('-' x $len);
					$pos+=1;
				}

				my $mut;
				$mut->{chr}  = $chr;
				$mut->{pos}  = $pos;
				$mut->{ref}  = $newref;
				$mut->{seq}  = $newalt;
				$mut->{type} = $type;
				$mut->{len}  = $len;
				$mut->{idv}  = $idv;
				$mut->{imf}  = $imf;
				$mut->{prevbpref} = $prevbpref;
				$mut->{prevbpalt} = $prevbpalt;
				$mut->{avgcov} = $idv;
				$mut->{mincov} = $idv;
				$mut->{status} = "ok";
				$mut->{zygosity} = "na";
				$mut->{altcov} = 0;
				$mut->{inheritance} = "na";

				leftNormalize(\$mut, 300, $REF, $FAIDX, $genome);
				#print STDERR "$mut->{chr}\t$mut->{pos}\t$mut->{ref}\t$mut->{seq}\t$mut->{type}\t$mut->{len}\n";

				#my $ref_encoded = md5_hex($mut->{ref});
				#my $qry_encoded = md5_hex($mut->{seq});

				#my $key = sprintf("%s:%d:%s:%d:%s:%s", $mut->{chr}, $mut->{pos}, $mut->{type}, $mut->{len}, $ref_encoded, $qry_encoded);
				my $key_long = sprintf("%s:%d:%s:%d:%s:%s", $mut->{chr}, $mut->{pos}, $mut->{type}, $mut->{len}, $mut->{ref}, $mut->{seq});
				my $key = md5_hex($key_long);

				if( !(exists $hash->{$key}) ) {
					$hash->{$key} = $mut;
					$cnt++;
				}
			}
        }
    }
    close VCF;

	print STDERR "$cnt mutations.\n";
}

# load the genome in fasta format
##########################################
sub loadGenomeFasta {

    print STDERR "Loading genome from FASTA file...";

    my $file = $_[0];
    my $genome = $_[1];

    open FASTAFILE, "< $file" or die "Can't open $file ($!)\n";

    my $header = "";
    my $SDNA = "";
    my $chr = "";
	my $cnt = 0;
    while ( <FASTAFILE> ) {
        chomp;
        #print "$_\n";
        if($_ =~ m/>/) {
			$cnt++;
            $header = $_;
            #print "$header\n";
            if($chr ne "") {
                #print "$chr\n";
                $genome->{$chr}->{seq} = $SDNA;
            }

            $header =~ s/^>//; # remove ">"
            $header =~ s/\s+$//; # remove trailing whitespace

            my ($label, $tmp) = split / /, $header, 2;
			$chr = $label;
			#if ($label =~ /^chr/) { $chr = substr($label,3); } # update chromosome label
            $SDNA = "";# clear out old sequence
        }
        else {
            s/\s+//g; # remove whitespace
            $SDNA .= $_;
        }
    }
    close FASTAFILE;

    if($chr ne "") { # handle last sequence
        #print "$chr\n";
        $genome->{$chr}->{seq} = $SDNA;
    }

	#if($VERBOSE) {
		print STDERR "$cnt sequences.\n";
	#}
}


## extract SM info from BAM header
#####################################################
sub extractSM {

	my $PREFIX = $_[0];
	my $SAMPLE = $_[1];
	my $headerfile = "header.txt";
	#my $headerfile = "header.".$SAMPLE.".txt";

	my $SM = "";

	## extract header from SAM file
	#runCmd("extract SAM header", "$bamtools header -in $PREFIX/bamfile.bam > $PREFIX/$headerfile");

	## read the list of read groups per sample from header
	open HEADER, "< $PREFIX/$SAMPLE/$headerfile" or die "Can't open $PREFIX/$SAMPLE/$headerfile ($!)\n";

	while (<HEADER>) {
		chomp;
		next if($_ eq ""); # skip empty string
		next if($_ !~ /^\@/); # skip non-header line

		my @records = split /\t/, $_;
	  	if ($records[0] eq "\@RG") {
			#@RG     ID:READGROUP    SM:SAMPLE

	  		#print join(' ', @records), "\n";

	  		my $sm;
	  		my $id;
	  		foreach my $tag (@records) {
	  			my @data = split /:/, $tag;
	  			#print join(':', @data), "\n";
	  			if($data[0] eq "SM") { $SM = $data[1]; }
	  			#if($data[0] eq "ID") { $id = $data[1]; }
	  		}
	  	}
	}
	close HEADER;

	return $SM;
}

## process SAM header to extract read groups info
#####################################################
sub parseHeader {

	my $PREFIX = $_[0];
	my $readgroups = $_[1];
	my $headerfile = "header.txt";

	# erease readgroup information from previous family
	for (keys %$readgroups) { delete $readgroups->{$_}; }

	## extract header from SAM file
	runCmd("extract SAM header", "$bamtools header -in $PREFIX/bamfile.bam > $PREFIX/$headerfile");

	## read the list of read groups per sample from header
	open HEADER, "< $PREFIX/$headerfile" or die "Can't open $PREFIX/$headerfile ($!)\n";

	while (<HEADER>) {
		chomp;
		next if($_ eq ""); # skip empty string
		next if($_ !~ /^\@/); # skip non-header line

		my @records = split /\t/, $_;
	  	if ($records[0] eq "\@RG") {
			#@RG     ID:READGROUP    SM:SAMPLE

	  		#print join(' ', @records), "\n";

	  		my $sm;
	  		my $id;
	  		foreach my $tag (@records) {
	  			my @data = split /:/, $tag;
	  			#print join(':', @data), "\n";
	  			if($data[0] eq "SM") { $sm = $data[1]; }
	  			if($data[0] eq "ID") { $id = $data[1]; }
	  		}
	  		push @{$readgroups->{$sm}->{rg}}, $id;
	  	}
	}
	close HEADER;
}

## save read groups to file
#####################################################
sub saveReadGroup {

	my $sample = $_[0];
	my $PREFIX = $_[1];
	my $readgroups = $_[2];
	my $rgfile = $_[3];

	open SAMPLE, "> $PREFIX/$rgfile" or die "Can't open $PREFIX/$rgfile ($!)\n";

	for my $sm ( keys %$readgroups ) {

		if($sample eq "ALL") {
			#print $sample.": ";
			foreach my $rg ( @{$readgroups->{$sm}->{rg}}) {
				#	print $rg.", ";
				print SAMPLE "$rg\n";
			}
		}
		elsif($sm eq $sample) {
			#print $sample.": ";
			foreach my $rg ( @{$readgroups->{$sm}->{rg}}) {
				#	print $rg.", ";
				print SAMPLE "$rg\n";
			}
		}
	}
	close SAMPLE;
}

1;
