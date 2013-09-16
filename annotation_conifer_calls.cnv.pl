#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email:  riverlee2008@gmail.com
# Date:   Wed Sep  4 10:17:19 2013
###################################
use strict;
use warnings;

# Taken Conifer CNV result as input, use bedtools to get the interval located genes
#
my $usage="perl $0 conifer.cnv.call.txt hg19.gtf output.txt\n";

my ($incnv,$genome,$output) = @ARGV;

die $usage if(@ARGV!=3 || ! -e $incnv);

# First step, convert cnv to bed format
info("Convert Conifer CNV result to bed format");
`sed -n '2,\$p' $incnv |awk 'BEGIN{OFS="\\t"} {print \$2,\$3,\$4,\$1,\$5}' >${incnv}.tmp.bed`;

# Running annotatePeaks.pl
info("Running intersectBed (bedtools)");
`intersectBed -a  ${incnv}.tmp.bed -b  $genome -wa -wb > ${incnv}.tmp2`;

# Reformat the result
my %result;
open(IN,"${incnv}.tmp2") or die $!;
while(<IN>){
    #chr10	116213601	117198403	RW7	dup	chr10	unknown	exon	116213601	116213870	.	-	.	gene_id "Myrfl"; gene_name "Myrfl"; p_id "P4872"; transcript_id "NM_001033333"; tss_id "TSS4260";
   s/\r|\n//g;  
   my ($chr,$start,$end,$sample,$dupdel,$tchr,$source,$class,$tstart,$tend,undef,$strand,undef,$details) = split "\t";
    my $g="";
    if($details=~/gene_name "(.*?)"/){
        $g=$1;
    }
    my $k=join "\t",($chr,$start,$end,$sample,$dupdel);
    $result{$k}->{$g}=1;
}

close IN;
info("Write out ");
open(OUT,">$output") or die $!;
foreach my $k (sort keys %result){
    my $str="";
    foreach my $ g (sort keys %{$result{$k}}){
        $str.="$g;";
    }
    print OUT $k,"\t",$str,"\n";
}
close OUT;

END{
    unlink "${incnv}.tmp.bed" if (-e "${incnv}.tmp.bed");
    unlink "${incnv}.tmp2" if (-e "${incnv}.tmp2");
}


sub info{
    my $str = shift;
    print "[",scalar(localtime),"] $str\n";
}
