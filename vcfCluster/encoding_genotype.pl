#!/usr/bin/perl
use strict;
use warnings;

=head1 description
 Take an vcf input file and coding the genotype as 0/0 (./.) as 0, 0/1 as 1 and 1/1 as 2
=cut

my $usage=<<USAGE;
Usage: perl $0 vcf output
take an vcf file and encode the genotype to 
0/0 => 0
0/1 => 1
1/1 => 2
the output format is
chr pos ref alt qual genotypes

USAGE

#die "perl $0 vcf output\n" if(@ARGV !=2 || !-e $ARGV[0]);
die $usage if(@ARGV !=2 || !-e $ARGV[0]);

#genotype encoding (skip those on x,y, and m)
#0/0 to 0
#0/1 to 1
#1/1 to 2
#
open(VCF,$ARGV[0]) or die $!;
open(OUT,">$ARGV[1]") or die $!;
while(<VCF>){
    next if (/^##/);
    s/\r|\n//g;
    if(/^#/){
        my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@samples) = split "\t";
        print OUT join "\t",($chr,$pos,$ref,$alt,$qual,@samples);
        print OUT "\n";
        next;
    }

    my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@samples) = split "\t";
    next if ($chr!~/\d/); #only chromosome 1-22
    print OUT join "\t",($chr,$pos,$ref,$alt,$qual);
    foreach my $g (@samples){
        my $v=0;
        $v=1 if($g=~/0\/1|1\/0/);
        $v=2 if($g=~/1\/1/);
        print OUT "\t",$v;
    }
    print OUT "\n";
}

close OUT;
close VCF;
