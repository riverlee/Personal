#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email: riverlee2008@gmail.com
# Date: Mon Aug  5 10:53:42 2013
###################################
use strict;
use warnings;

# Read the VCF file and get the annotation from anova geno
# Note the order of vcf file is the same as those for the annovo input
#
die "$0 in.vcf in.variant_function in.exonic_variant_function\n" if (@ARGV !=3);

my ($invcf,$inv,$inexon) = @ARGV;

my @variants=();
my %exon=();
# Read the in.variant_function file
open(INV,$inv) or die $!;
while(<INV>){
    my($region,$details,@a) = split "\t";
    my $key = join "\t",($region,$details);
    push @variants,$key;
}
close INV;


# Read the in.exonic_variant_function file
open(IN,$inexon) or die $!;
while(<IN>){
    my($line,$type,$detail,@a) = split "\t";
    $line=~s/line//g;
    $line=$line-1;  # use 0-based for index
    $exon{$line}=join "\t",($type,$detail);
}
close IN;


open(IN,$invcf) or die $!;
open(OUT,">${invcf}_add_anno.txt") or die $!;

my $row=0;
while(<IN>){
    next if (/^##/);
    if(/^#/){
        print OUT $_;
        next;
    }
    s/\r|\n//g;
    print OUT join "\t", ($_,$variants[$row]);
    if(exists($exon{$row})){
        print OUT "\t",$exon{$row};
    }
    print OUT "\n";
    $row++;
}
close IN;
close OUT;
