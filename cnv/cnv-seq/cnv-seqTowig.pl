#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email:  riverlee2008@gmail.com
# Date:   Tue Aug 27 16:05:17 2013
###################################
use strict;
use warnings;

my $usage=<<USG;
perl $0 cnv-seq.cnv out.wig [name]

if name is provided, will added the track in the wig output file

USG

die $usage if(@ARGV <=2);
my ($in,$out,$name) = @ARGV;

open(IN,$in) or die $!;
open(OUT,">$out") or did $!;
if(defined($name)){
    print OUT "track type=wiggle_0 name=$name description=$name\n";
}
<IN>;
my %dat;
while(<IN>){
    s/\r|\n|"//g;
    my($chr,$start,$end,$test,$ref,$position,$log,@a) = split "\t";
    next if ($log=~/inf|NA/i);
    # convert to integer
    $position=sprintf("%d",$position);
    $dat{$chr}->{$position}=$log;
}
close IN;
foreach my $chr (sort keys %dat){
    print OUT "variableStep chrom=$chr\n";
   foreach my $p (sort {$a<=>$b}  keys %{$dat{$chr}}){
        print OUT $p,"\t",$dat{$chr}->{$p},"\n";
   }
}
close IN;
