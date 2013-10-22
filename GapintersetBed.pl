#!/usr/bin/env perl
###################################
# Author: Jiang (River) Li
# Email:  riverlee2008@gmail.com
# Date:   Wed Sep 18 12:30:23 2013
###################################
use strict;
use warnings;

my ($abed,$bbed,$gap,$out) = @ARGV;
my $usage="Usage: perl $0 bed1 bed2 gap out\n";
$usage.="Details: will extended bed2 up/downstream gap and then invoke\n";
$usage.="intersectBed -a bed1 -b bed2 -wa |sort -k1,1 -n2,2k |uniq >out\n";

die $usage if (@ARGV !=4);

my $rand = rand();
my $bed1 = "$rand.tmp1.bed";
my $bed2 = "$rand.tmp2.bed";
# sort bed1
`sort -k1,1 -k2,2n $abed > $bed1`;

# extend bed2
open(OUT,">$bed2") or die $!;
open(IN,"sort -k1,1 -k2,2n $bbed|") or die $!;
while(<IN>){
    s/\r|\n//g;
    my($chr,$start,$end) = split "\t";
    $start=$start-$gap;
    $end=$end+$gap;
    $start=1 if ($start<1);
    print OUT join "\t",($chr,$start,$end."\n");
}
close IN;
close OUT;

`intersectBed -a $bed1 -b $bed2 -wa |sort -k1,1 -k2,2n |uniq >$out`;


unlink $bed1 if (-e $bed1);
unlink $bed2 if (-e $bed2);



