#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email: riverlee2008@gmail.com
# Date: Wed Mar 12 11:32:57 2014
###################################
use strict;
use warnings;

# http://sift.bii.a-star.edu.sg/

### SIFT
## SNV INPUT 
# chromosome,coordinate,oientation,alleles,user comment(optional)
#Format Example 1 (for NCBI Human assembly version 37): RESIDUE BASED COORDINATE SYSTEM (comma separated) 
#1,100382265,1,C/G 
#1,100380997,1,A/G 

my $usage = "vcf2SIFT.pl in.vcf\n";
die $usage if (@ARGV !=1);
open(IN,"$ARGV[0]") or die $!;
open(SNV,">SNV_SIFT_Input.txt") or die $!;
open(INDEL,">INDEL_SIFT_Input.txt") or die $!;


while(<IN>){
    next if (/^#/);
    my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format) = split "\t";
    
    my $ori=join "|",("#$chr",$pos,$id,$ref,$alt);
    $chr=modifiyChr($chr);

    # Determine SNV or Indel
    my $snv=1;
    if(length($ref)>1 || ($alt!~/,/ && length($alt)>1)){
        $snv=0; # Indel
    }

    if($snv){
        print SNV join ",",($chr,$pos,"1","$ref/$alt\n");
    }else{
       if(length($ref)==1){
            # Insertion  
            print INDEL join ",",($chr,$pos,$pos,"1","-/".substr($alt,1),$ori."\n");
       }else{
            # Deletion
            print INDEL join ",", ($chr,$pos,$pos+length($ref)-1,"1",substr($ref,1)."/-","$ori\n");
       }
    }
}
sub modifiyChr{
    my $chr = shift @_;
    $chr=~s/chr//gi;
    $chr=uc($chr);
    if($chr eq "M"){
        $chr="MT";
    }
    return $chr;
}
