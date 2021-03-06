#!/usr/bin/perl

BEGIN {
	use FindBin '$RealBin';
	push @INC, $RealBin;
}

use Getopt::Std; %opt = ();

# defaults
$adapters = "$RealBin/SciMETv2_Adapters.fa";
#$adapters = "$RealBin/sciMETv2_adapters.fa";
$trimmomatic = "trimmomatic";
$min_RL = 30;
$threads = 1;
$r2_trim = 10;

getopts("O:A:1:2:T:m:t:e:", \%opt);

$die = "

sciMETv2_trim.pl (options) -O output_prefix -1 read1.fq.gz -2 read2.fq.gz

Runs trimmomatic with sciMET parameters.

Options:

-A   [STR]   Adapter fastq (def = $adapters)
-O   [STR]   Output Prefix
-1   [STR]   Read1 fastq
-2   [STR]   Read2 fastq
-T   [STR]   Trimmomatic jar file path (def = $trimmomatic)
-m   [INT]   Min read length (def = $min_RL)
-t   [INT]   Threads to use (def = $threads)
-e   [INT]   Trim bases from the end of read 2 after adapter trim (def = $r2_trim)
             If set to 0 will not trim any additional bases.

";

if (!defined $opt{'O'}) {die "\nERROR: Specify output prefix as -O\n$die"};
if (!defined $opt{'1'} && !defined $opt{'2'}) {die "\nERROR: Specify at least one read as -1 or -2\n$die"};
if (defined $opt{'A'}) {$adapters = $opt{'A'}};
if (-e "$adapters") {print STDERR "Adapter file confirmed!\n"} else {die "ERROR: Adapter file $adapters does not exist!\n$die"};
if (defined $opt{'T'}) {$trimmomatic = $opt{'T'}};
if (defined $opt{'m'}) {$min_RL = $opt{'m'}};
if (defined $opt{'t'}) {$threads = $opt{'t'}};
if (defined $opt{'e'}) {$r2_trim = $opt{'e'}};

if (defined $opt{'1'}) {
	$trim_R1 = "$trimmomatic SE -threads $threads $opt{'1'} $opt{'O'}.trimmed.R1.fq.gz ILLUMINACLIP:$adapters:2:30:10 MINLEN:$min_RL >> $opt{'O'}.trimmed.R1.log 2>> $opt{'O'}.trimmed.R1.log";
	print STDERR "Trimming read 1. Command:\n\t$trim_R1\n";
	system("$trim_R1");
}

if (defined $opt{'2'}) {
	$trim_R2 = "$trimmomatic SE -threads $threads $opt{'2'} $opt{'O'}.trimmed.R2.fq.gz ILLUMINACLIP:$adapters:2:30:10 MINLEN:$min_RL >> $opt{'O'}.trimmed.R2.log 2>> $opt{'O'}.trimmed.R2.log";
	print STDERR "Trimming read 2. Command:\n\t$trim_R2\n";
	system("$trim_R2");
	if ($r2_trim > 0) {
		$trim_R2e = "cutadapt --cores $threads -u $r2_trim -m $min_RL -o $opt{'O'}.trimmed_e.R2.fq.gz $opt{'O'}.trimmed.R2.fq.gz";
		print STDERR "Trimming additional $r2_trim from end of read 2. Command:\n\t$trim_R2e\n";
		system("$trim_R2e");
		system("rm -f $opt{'O'}.trimmed.R2.fq.gz && mv $opt{'O'}.trimmed_e.R2.fq.gz $opt{'O'}.trimmed.R2.fq.gz");
	}
}

exit;
