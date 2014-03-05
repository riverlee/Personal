#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email:  riverlee2008@gmail.com
# Date:   Wed Mar  5 10:25:51 2014
###################################
use strict;
use warnings;

# After applying annovaReAnnoHash.pl, try to get the functional variants

die "getFunctional.pl merged.vcf_add_anno.txt\n" if (@ARGV!=1);

my $in = shift @ARGV;
open(IN,$in) or die $!;
open(OUT,">${in}_functional.txt") or die $!;
while(<IN>){
    if(/^#/){
        print OUT $_;
        next;
     }
     my @a = split "\t";
     my $details=pop @a;
     my $class=pop @a;
     my $gene=pop @a;
     my $region=pop @a;
     if($class eq "stoploss SNV" || $class eq "stopgain SNV" ||
        $class eq "frameshift insertion" || $class eq "frameshift deletion" ||
        $class eq "nonsynonymous SNV"){
        print OUT $_;
     }
}

