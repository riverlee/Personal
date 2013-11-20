#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email:  riverlee2008@gmail.com
# Date:   Wed Nov 20 13:59:43 2013
###################################
use strict;
use warnings;

# Add SNP ID to the ID columns
my $script=$0;
$script=~s/.*\/|.*\\//g;
my $usage="$script in.vcf variant.avinput.hg19_snp138_dropped\n";

die $usage if (@ARGV !=2);

my ($in,$snpdrop) = @ARGV;

die "$in is not exists \n" if (! -e $in);
die "$snpdrop is not exists \n" if (! -e $snpdrop);

# Read snp id
my %snp;
open(IN,$snpdrop) or die $!;
while(<IN>){
    s/\r|\n//g;
    my ($db,$id,$chr,$start,$end,$ref,$alt) = split "\t";
    my $k=join "\t",($chr,$start);
    $snp{$k} = $id;
}
close IN;


# Make a copy
`cp $in ${in}.bk`;
open(IN,"${in}.bk") or die $!;
open(OUT,">$in") or die $!;
while(<IN>){
    if(/^#/){
        print OUT $_;
    }else{
        my($chr,$pos,$id,@a) = split "\t";
        my $k=join "\t",($chr,$pos);
        $id=$snp{$k} if (exists($snp{$k}));
        print OUT join "\t",($chr,$pos,$id,@a);
    }
}
close IN;
close OUT;

`rm ${in}.bk`;


