#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email: riverlee2008@gmail.com
# Date: Wed Jul 24 11:09:08 2013
###################################
use strict;
use warnings;
use Getopt::Long;

# Given a bed format file (alignment) and center position file (chr\t pos).
# Get the signal around this center
#

my ($alnfile,$centerfile,$output,$extend,$window,$normalized)=("","","",8000,10,1);
my $help=0;
my $iscenterbed=0; # default the input center file is in format chr\tpos, if bed format, get the center by start+end/2
my $isalnbam=0;    # default the alnment is in bed format. 
my $readlen=36;   # default assume read length is 36bp
unless(GetOptions("a=s"=>\$alnfile,
           "c=s"=>\$centerfile,
           "o=s"=>\$output,
           "e=i"=>\$extend,
           "w=i"=>\$window,
           "n=i"=>\$normalized,
           "b=i"=>\$iscenterbed,
           "bam=i"=>\$isalnbam,
           "len=i"=>\$readlen,
           "h"=>\$help)){
    help();
    die $!;
 }

help(1) if ($help);
help(1) if($alnfile eq "" && $centerfile eq "" && $output eq "");
die ("$alnfile not exists\n") if(! -e $alnfile);
die ("$centerfile not exists\n") if(! -e $centerfile);


my %aln;

# First read the center file and store the position to aln
info("Initialize the signal matrix ...");
open(FIRST,$centerfile) or die $!;
while(<FIRST>){
    s/\r|\n//g;
    my ($chr,@a) = split "\t";
    my $center="";
    if($iscenterbed){
        my ($start,$end) = @a; # bed format
        $center = int (($start+$end)/2);
    }else{
        $center = $a[0];
    }
    $chr=~s/chr//g;  #remove chr
    my $extend_start = $center-$extend/2;
    my $extend_end = $extend_start+$extend;
 
    # using a 10bp window to get the signal for each of this 10bp region across the 8kb region.
    for(my $i=$extend_start;$i<$extend_end;$i++){
        #calculate the signal
        $aln{$chr}->{$i}=0;
    }
}
close FIRST;

#Load alignment
info("Load alignment");
my $totalreads=0;
my $flag_read_unmapped = 0x0004;

if($isalnbam){
    if($alnfile=~/bam$/i){
        open(IN,"samtools view $alnfile |") or die $!;        
    }else{
        open(IN,"$alnfile") or die $!;
    }
    while(<IN>){
        next if (/^@/);
        my($rname,$flag,$chr,$pos,$mapQ,$cigar,@others) = split "\t";
        next if ($flag & $flag_read_unmapped);
        $totalreads++;

        unless($totalreads %100000){
            info("Reading $totalreads");
        }

        $chr=~s/chr//g;  #remove chr
        my $start = $pos;
        my $end = $start+$readlen;
#        ($start,$end)=sort {$a<=>$b} ($start,$end); #in case $end < $start
        for (my $i=$start;$i<=$end;$i++){
            $aln{$chr}->{$i}++ if(exists($aln{$chr}->{$i}));  #the depth for each base
        }
    }
}else{
    # alignment is in bed format
    open(IN,"$alnfile") or die $!;

    while(<IN>){
        s/\r|\n//g;
        $totalreads++;
        unless($totalreads %100000){
            info("Reading $totalreads");
        }
    
        my($chr,$start,$end) = split "\t";
        $chr=~s/chr//g;  #remove chr
        ($start,$end)=sort {$a<=>$b} ($start,$end); #in case $end < $start
        for (my $i=$start;$i<=$end;$i++){
            $aln{$chr}->{$i}++ if(exists($aln{$chr}->{$i}));  #the depth for each base
        }
    }
    close IN;
}

# read the peak region
my @signal;
my $index=0;
open(PEAK,$centerfile) or die $!;
info("Loading peakfile '$centerfile' and calculate the signal matrix");
while(<PEAK>){
    s/\r|\n//g;
    my ($chr,@a) = split "\t";
    my $center="";
    if($iscenterbed){
        my ($start,$end) = @a; # bed format
        $center = int (($start+$end)/2);
    }else{
        $center = $a[0];
    }
    $chr=~s/chr//g;  #remove chr
    info("Doing $index ");
   
    my $extend_start = $center-$extend/2;
    my $extend_end = $extend_start+$extend;
    
    # using a 10bp window to get the signal for each of this 10bp region across the 8kb region.
    for(my $i=$extend_start;$i<$extend_end;$i+=$window){
        #calculate the signal
        my $intensity=0;
        for (my $j=$i;$j<$i+$window;$j++){
            if(exists($aln{$chr}->{$j})){
                $intensity+=$aln{$chr}->{$j};
            }
        }
        $intensity=$intensity/$window;
        if($normalized){
            $intensity=($intensity/$totalreads)*1e9
        }
        push @{$signal[$index]->{'regions'}}, {'signal'=>$intensity};
        #push @{$signal[$index]}
    }
    $index++;
}
close PEAK;

#Write out
info("Write out signal matrix");
open(OUT,">$output") or die $!;
for(my $i=0;$i<@signal;$i++){
    my $jj=$i+1;
    print OUT "peak${jj}";
    for (my $j=0;$j<@{$signal[$i]->{'regions'}};$j++){
       print OUT "\t",$signal[$i]->{'regions'}->[$j]->{'signal'};
    }
    print OUT "\n";
}

close OUT;




sub help{
    print <<HELP;
Usage: getCenteredSig.pl -a <alignfile> -c <centerfile> -o <output> -e [extend bp] -w [window size] -n/-non
                -a  alignment file in bed format
                -c  center position file, format is "chr    pos"
                -b  is the center input file in bed format, thus will estimate the center by (start+end)/2, defalut 0, values is 0 or 1
                -o  output file name
                -e  extend length around the center, extend by upstream extend/2 bp and downstream extend/2bp, default 8000bp
                -w  window size to get the signal, default is 10,(1 means get the average coverage for each base around the extended region)
                -n  normalized the signal to RPKM or use the raw count, values is 0 or 1
                -bam whether input alignment file is in bam/sam format, values is 0 or 1,default 0
                -len the average read length, default is 36
                -h  display this help message
HELP
    exit(1) if (shift @_);
}


sub info{
    my $str=shift;
    print "[",scalar(localtime),"] $str \n";
}



