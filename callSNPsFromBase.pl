#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email:  riverlee2008@gmail.com
# Date:   Tue Nov 19 13:07:28 2013
###################################
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;

# Custom SNP caller based the frequency of A/T/C/G
# 1) First you run samtools mpileup 
# 2) Then you apply pileup2base_no_strand.pl to get frequency of A/T/C/G
# Check https://github.com/riverlee/pileup2base/

my $infile=""; # base file from pileup2base_no_strand.pl output
my $depth=200; # Above $depth we use $depth*$ratio_cutoff as cutoff
my $absolute_cutoff=10; # When depth below $depth, use $absolute_cutoff as cutoff
my $ratio_cutoff=0.05;  # 
my $output="";
my $help=0;
my $name="";

unless(
    GetOptions(
        "i=s"=>\$infile,
        "o=s"=>\$output,
        "d=i"=>\$depth,
        "a=i"=>\$absolute_cutoff,
        "r=f"=>\$ratio_cutoff,
        "n=s"=>\$name,
        "h"=>\$help
    )
){
    print $!,"\n";
    usage(1);
}

check();
if($name eq ""){
    $name=basename($infile);
    print $name,"\n";
}

## Main Program
open(IN,$infile) or die $!;
open(OUT,">$output") or die $!;
my $header="##fileformat=VCFv4.1\n";
$header.="##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n";
$header.="##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"# bases of A,T,C,G\">\n";
$header.="##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
$header.="##FORMAT=<ID=DP2,Number=2,Type=Float,Description=\"Alleles for ref,alt\">\n";
$header.=join "\t",("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",$name);
$header.="\n";
print OUT $header;
<IN>;
info("Start SNP calling ...");
my $count=0;
while(<IN>){
#chr  loc ref A   T   C   G
    $count++;
    if($count % 5000 ==0){
        print ".";
    }
    s/\r|\n//g;
    my($chr,$loc,$ref,$A,$T,$C,$G) = split "\t";
    my $t=$A+$T+$C+$G;
    my $info_dp=$t;
    my $info_dp4="$A,$T,$C,$G";
    my $info="$info_dp;$info_dp4";

    my %tmp;
    $tmp{A}=$A;$tmp{C}=$C;$tmp{T}=$T;$tmp{G}=$G;

    my $cutoff=$absolute_cutoff; # default )use absolute_cutoff
    if($t>$depth){
        $cutoff=$t*$ratio_cutoff;
    }
    my $format="GT:DP2";
    
    foreach my $tt (sort {$tmp{$b}<=>$tmp{$a}} keys %tmp){
        next if ($tt eq uc($ref));
        if($tmp{$tt}>$cutoff){
            #determine genotype 0/1 when %alt in [0,75%]
            my $geno="0/1";
            if($tmp{$tt}/$t > 0.75){
                $geno="1/1";
            }
            print OUT join "\t",($chr,$loc,".",uc($ref),$tt,255,".",$info,$format,"$geno:$tmp{uc($ref)},$tmp{$tt}\n");
        }
    }
}
print "\n";
info("Finished");
close IN;
close OUT;

sub usage{
    my ($flag) = @_;
    my $p = $0;
    $p=~s/.*\/|.*\\//g;
    print <<USAGE;
Usage:perl $p -i infile -o outfile -d depth -a absolute_cutoff -r ratio_cutoff
      -i input file
      -o output file
      -d depth used to select alternative allele count cutoff
         if >depth, use depth*ratio_cutoff as cutoff,
         otherwise, use absolute_cutoff as cutoff. default is 200.
      -a absolute alternative allele count cutoff, default is 10.
      -r the ratio of alternative allele as cutoff,default is 0.05
      -n name in the vcf column, default is use the input filename
      -h help out this help

USAGE
    exit 1 if ($flag);
}

sub check{
    my $message="";
    if($infile eq ""){
        $message.="infile (-i ) is not defined\n";
    }elsif(! -e $infile){
        $message.="'$infile' is not exists\n";
    }
    if($output eq ""){
        $message.="output (-o ) is not defined\n";
    }

    if($help || $message ne ""){
        print $message;
        usage(1);
    }

}

sub info{
    my ($str,$flag) = @_;
    print "[",scalar(localtime),"] $str ";
    if($flag){ 
        # No newline
    }else{
        print "\n";
    }
}
