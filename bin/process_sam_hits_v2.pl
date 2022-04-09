#!/usr/bin/perl -w

# Contact: Dinh Diep
# Version 1.5

use strict;
use warnings;

my $TMP_DIR;
my $qual_base = 33;
my $encoded_offset = 45;
my $min_mapq = 10;

my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
$rcTable{'N'}='N';
$rcTable{'R'}='Y';
$rcTable{'Y'}='R';
$rcTable{'M'}='K';
$rcTable{'K'}='M';
$rcTable{'S'}='S';
$rcTable{'W'}='W';

my $missed_guesses = 0;
my $total_mreads = 0;
my $total_bases = 0;
my $total_mbases = 0;

&main;
exit;

sub main(){
    my $map_file_CT = $ARGV[0];
    my $map_file_GA = $ARGV[1];
    my $out_file = $ARGV[2];
    my $unmapped_file = $ARGV[3];
    sortsam($map_file_CT, $map_file_GA, $out_file, $unmapped_file);
    print "Total mapped sequences\t$total_mreads\n";
    print "Total mapped bases\t$total_mbases\n";
    print "Total wrong guesses\t", $missed_guesses, "\n";
}

sub processhit{
    my $array = shift;
    my @f = @{$array};
    
    my @qual = split "", $f[10];
    my @encoded_seq = split "", $f[9];
    for(my $i = 0; $i < scalar(@qual); $i++){
       if(ord($qual[$i]) > 78){
           if($encoded_seq[$i] =~ m/T/i){
              $encoded_seq[$i] = 'C';
           }else{
              $encoded_seq[$i] = 'G';
           }
           my $old_qual = ord($qual[$i]) - $encoded_offset;
           $qual[$i] = chr($old_qual);
       }
    }
    $f[9] = join("", @encoded_seq);
    $f[10] = join("", @qual);

    $f[2] =~ s/_(CT|GA)_converted$//;    
    ###-----------------BEGIN deal with CIGAR---------------------###
    my @CIGARS;
    # I - insertion to the reference
    # D - deletion to the reference
    # N - skipped region from the reference
    # S - soft clip on the read
    # H - hard clip on the read
    # P - padding - silent deletion from the padded reference sequence


    $total_mbases+=length($f[9]);
    return @f[0...10];
}

sub sortsam{
    my $map_file_CT = shift;
    my $map_file_GA = shift;
    my $out_file = shift;
    my $unmapped_file = shift;
    open(SAM_IN, "cat $map_file_CT $map_file_GA | sort -k1,1d |") || die("Error in opening $map_file_CT and $map_file_GA.");
    system("samtools view -H $map_file_CT | sed 's/_CT_converted//g' > $out_file");
    open(SAM_OUT, ">>$out_file") || die("Error writing to file $out_file.\n");
    open(FASTQ_OUT, ">$unmapped_file") || die("Error writing to file $unmapped_file.\n");
    my $last_line = 0;
    my @last_fields;
    my $last_cnt = 0;
    while(my $line =  <SAM_IN>){
        my @fields = split(/\t/, $line);
        ##print $fields[0], "\n";
        next if($fields[0] =~ /^@/);
        if(!$last_line){
	   $last_line = $line;
	   @last_fields = @fields;
	   $last_cnt = 0;
	   next;
        }
        if($fields[0] eq $last_fields[0]){
	   my $score = $fields[4];
	   my $last_score = $last_fields[4];
	   if($score > $last_score){
	       #2nd line is a better hit
	       $last_line = $line;
	       @last_fields = @fields;
	       $last_cnt = 0;
	   }elsif($score == $last_score){
	       #check if one is from alternate reference (alternate references are not listed in fai file)
	       #two equivalent good hits, increment last_cnt
	       $last_cnt++;
	   }else{
	       #1st line is a better hit, do nothing
	   }
        }else{	
	   if($last_cnt ne 0 || $last_fields[4] < $min_mapq){ # not a unique best hit
               print FASTQ_OUT "@", $last_fields[0], "\n", $last_fields[9], "\n+\n", $last_fields[10], "\n";
	       undef $last_line;
	       undef @last_fields;
	       next;
	   }
	   if(my @hit = processhit(\@last_fields)){
	       print SAM_OUT join("\t", @hit), "\n";
               $total_mreads++;
	   }else{
               if($last_fields[1] & 16){
                   # reads mapped in reverse
                   my $orig_seq = revComp($last_fields[9]);
                   my $orig_qual = reverse($last_fields[10]);
                   print FASTQ_OUT "@", $last_fields[0], "\n", $orig_seq, "\n+\n", $orig_qual, "\n";
               }else{
                   print FASTQ_OUT "@", $last_fields[0], "\n", $last_fields[9], "\n+\n", $last_fields[10], "\n";
               }
           }
	   $last_line = $line;
	   @last_fields = @fields;
        }
    }
    #print the last line
    if($last_line and $last_cnt eq 0 and $last_fields[4] > $min_mapq){
        if(my @hit = processhit(\@last_fields)){
	   print SAM_OUT join("\t", @hit), "\n";
           $total_mreads++;
        }else{
           if($last_fields[1] & 16){
              # reads mapped in reverse
              my $orig_seq = revComp($last_fields[9]);
              my $orig_qual = reverse($last_fields[10]);
              print FASTQ_OUT "@", $last_fields[0], "\n", $orig_seq, "\n+\n", $orig_qual, "\n";
           }else{
              print FASTQ_OUT "@", $last_fields[0], "\n", $last_fields[9], "\n+\n", $last_fields[10], "\n";
           }
        }
    }else{
        if($last_fields[1] & 16){
           # reads mapped in reverse
           my $orig_seq = revComp($last_fields[9]);
           my $orig_qual = reverse($last_fields[10]);
           print FASTQ_OUT "@", $last_fields[0], "\n", $orig_seq, "\n+\n", $orig_qual, "\n";
        }else{
           print FASTQ_OUT "@", $last_fields[0], "\n", $last_fields[9], "\n+\n", $last_fields[10], "\n";
        }
    }
    close(SAM_IN);
    close(SAM_OUT);
    close(FASTQ_OUT);
}

sub revComp{
    my $seq = shift;
    my $rcSeq='';
    for(my $i=0; $i<length($seq); $i++){
        $rcSeq = $rcTable{uc(substr($seq,$i,1))} . $rcSeq;
    }
    return $rcSeq;
}
