#!/usr/bin/perl -w

use strict;

my $total_bases = 0;
my $c_or_g_bases = 0;
while(my $line = <STDIN>){
  chop($line);
  my @f = split(/\t/, $line);
  $total_bases += length($f[9]);
  my $string = $f[9];
  my $count = 0;
  if($f[1] & 16){
    # reverse string
    $count = () = $string =~ /g|G/g;
  }else{
    $count = () = $string =~ /c|C/g;
  }
  $c_or_g_bases += $count;
}

print "$ARGV[0]\t$total_bases\t$c_or_g_bases\n";

