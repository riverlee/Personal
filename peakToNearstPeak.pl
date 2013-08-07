#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email: riverlee2008@gmail.com
# Date: Mon Aug  5 19:58:54 2013
###################################
use strict;
use warnings;
use Getopt::Long;

my ($inbed,$searchbed,$output) = ("","","");
my $help=0;

unless(GetOptions("i=s"=>\$inbed,
           "d=s"=>\$searchbed,
           "o=s"=>\$output,
           "h"=>\$help)){
    help();
    die $!;
 }

help(1) if ($help);
help(1) if($inbed eq "" && $searchbed eq "" && $output eq "");
die ("$inbed not exists\n") if(! -e $inbed);
die ("$searchbed not exists\n") if(! -e $searchbed);

info("Loading backgroud peaks coordinates");
my %background;
open(IN,"$searchbed") or die $!;
while(<IN>){
    s/\r|\n//g;
    my($chr,$start,$end) = split "\t";
    my $center = int (($start+$end)/2);
    my $key="$chr:$start-$end";
    $background{$chr}->{$center}=$key;
}
close IN;

info("Assign to the nearest peaks");
my $total=0;
open(IN,"$inbed") or die $!;
open(OUT,">$output") or die $!;
my $maxdis=10000000;
while(<IN>){
    s/\r|\n//g;
    my($chr,$start,$end) = split "\t";
    my $center = int (($start+$end)/2);
    $total++;
    unless($total %100){
        info("Reading $total");
    }
    if(!exists($background{$chr})){
        print OUT join "\t",($chr,$start,$end,$maxdis,"NA\n");
    }else{
        my @pos = sort keys %{$background{$chr}};
        my %tmp;
        foreach my $p (@pos){
            my $diff = abs($p-$center);
            $tmp{$diff}=$p;
        }
        my @pp = sort {$a<=>$b} keys %tmp;
        my $pos = $tmp{ $pp[0]};
        print OUT join "\t",($chr,$start,$end,$pp[0],$background{$chr}->{$pos});
        print OUT "\n";
    }
}
close IN;

close OUT;
sub help{
    print <<HELP;
Usage:  peakToNeastPeak.pl -i <in.peak.bed> -s <search.peak.bed> -o <output> 
                -i  input peak file in bed format,will only use chr,start,end columns
                -d  another peak file which will assign in.peak.bed to.
                -o  output file name
                -h  display this help message
HELP
    exit(1) if (shift @_);
}


sub info{
    my $str=shift;
    print "[",scalar(localtime),"] $str \n";
}



