#!/usr/bin/env perl
###################################
# Author: Jiang (River) Li
# Email:  riverlee2008@gmail.com
# Date:   Thu May 30 15:08:23 2013
###################################
use strict;
use warnings;

# Fetch the significant interaction pairs
my($infile,$cutoff,$output) = @ARGV;
if(@ARGV !=3){
    die "perl $0 in.cor.matrix.csv 0.8 output.txt\n";
}

open(IN,$infile) or die $!;
my $firstline=<IN>;
$firstline=~s/\r|\n|"//g;
my @header=split ",",$firstline;
shift(@header);
my %pairs;
my $count=0;
while(<IN>){
    s/\r|\n|"//g;
    my($g,@cor) = split ",";
    for(my $i=$count+1;$i<@cor;$i++){
        if(abs($cor[$i])>$cutoff){
            my($a,$b) = sort ($g,$header[$i]);
            $pairs{$a}->{$b}=$cor[$i];
        }
    }
    $count++;
}
close IN;

# Write OUT ;
open(OUT,">$output") or die $!;
foreach my $a (sort keys %pairs){
    foreach my $b(sort keys %{$pairs{$a}}){
        print OUT join "\t",($a,$b,$pairs{$a}->{$b});
        print OUT "\n";
    }
}
close OUT;
 
