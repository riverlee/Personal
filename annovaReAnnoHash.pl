#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Emaii: riverlee2008@gmail.com
# Date: Mon Aug  5 10:53:42 2013
###################################
use strict;
use warnings;

# Read the VCF file and get the annotation from anova geno
# Note the order of vcf file is the same as those for the annovo input

# Wed Feb 26 16:20:20 PST 2014
# Use hash instead of array to map
# 

die "$0 in.vcf in.variant_function in.exonic_variant_function\n" if (@ARGV !=3);

my ($invcf,$inv,$inexon) = @ARGV;

my %variants=();
my %exon=();

# Read the in.variant_function file

open(INV,$inv) or die $!;
while(<INV>){
    my($region,$details,$chr,$start,$end,$ref,$alt,@others) = split "\t";
    my $key = join "\t",($region,$details);
    my $kk = join "\t",($chr,$start,$end,$ref,$alt);
    $variants{$kk} = $key;
}
close INV;


# Read the in.exonic_variant_function file
open(IN,$inexon) or die $!;
while(<IN>){
    my($line,$type,$detail,$chr,$start,$end,$ref,$alt,@a) = split "\t";
    my $kk = join "\t",($chr,$start,$end,$ref,$alt);
    $exon{$kk}=join "\t",($type,$detail);
}
close IN;


open(IN,$invcf) or die $!;
open(OUT,">${invcf}_add_anno.txt") or die $!;

while(<IN>){
    next if (/^##/);
    if(/^#/){
        s/\r|\n//g;
        print OUT $_;
        print OUT "\t";
        print OUT join "\t",("region","gene","class","details\n");
        next;
    }
    s/\r|\n//g;

    my ($region,$gene,$class,$details)=("NA","NA","NA","NA");

    # Take the following by  from line 1215
    my ($chr, $start, $ID, $ref_allele, $mut_allele, $quality_score, $filter, $info, $format, $sample) = split "\t";
    my $end;

    $ref_allele = uc $ref_allele;
    $mut_allele = uc $mut_allele;

    if ($mut_allele eq '.') {           #no variant call was made at this position
        next;
    }

    my $mut_allele2;
    if ($mut_allele =~ m/([^,]+),([\w,]+)/) {   #there could be more than two alternative alleles
        $mut_allele = $1;
        $mut_allele2 = $2;
    }

    if(length($ref_allele)==1 && length($mut_allele)==1) {  	### output snv
        $end = $start;
    } elsif (length($ref_allele) > 1 || length($mut_allele) > 1) {  ### output indel
        #example VCF4 records below:
        #20      2       .       TCG     T       .       PASS    DP=100
        #Chr1    5473    .       AT      ATT     23.5    .       INDEL;DP=16;AF1=0.5;CI95=0.5,0.5;DP4=4,2,3,1;MQ=42;PV4=1,0.41,0.042,0.24
        #Chr1    6498    .       ATTTT   ATTTTT  53.5    .       INDEL;DP=9;AF1=1;CI95=1,1;DP4=0,0,5,3;MQ=28
        if(length($ref_allele) > length ($mut_allele)) { 		# deletion or block substitution
            my $head = substr($ref_allele, 0, length ($mut_allele));
            if ($head eq $mut_allele) {
                my $startbk = $start;
                $start = $start+length($head);
                $end = $startbk+length($ref_allele)-1;

                $ref_allele = substr ($ref_allele, length ($mut_allele));
                $mut_allele = "-";
            } else {
                $end = $start + length($ref_allele)-1;
            }
        } elsif(length($mut_allele) >= length ($ref_allele)) { 		# insertion or block substitution
            my $head = substr ($mut_allele, 0, length ($ref_allele));
            if ($head eq $ref_allele) {
                $start = $start+length($ref_allele)-1;
                $end = $start;
                $mut_allele = substr ($mut_allele, length ($ref_allele));
                $ref_allele = "-";
            } else {
                $end = $start + length($ref_allele)-1;
            }
        }
    }

    my $kk = join "\t",($chr,$start,$end,$ref_allele,$mut_allele);

    if(exists($variants{$kk})){
        ($region,$gene) = split "\t",$variants{$kk};
    }
    if(exists($exon{$kk})){
        ($class,$details) = split "\t",$exon{$kk};
    }

    print OUT join "\t",($_,$region,$gene,$class,$details); 
    print OUT "\n";
}
close IN;
close OUT;


