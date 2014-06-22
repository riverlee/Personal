#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email: riverlee2008@gmail.com
# Date: Fri Jun 13 11:24:06 2014
###################################
use strict;
use warnings;

## Base on refseq transcript table which is downloaded from UCSC, annotate peaks to 4 levels (Promoter (-2kb to +1kb), gene body, Enhancer (-2kk to -100kb) and intergenic regions).
#

my ($inbed,$refseq,$outfile) = @ARGV;

## golobal variants
our $promoterup=-2000;
our $promoterdown=1000;

our $enhancerup=-100000;
our $enhancerdown=-2000;


## 1) load refseq information
print "[",scalar(localtime),"] Load RefSeq ...\n";
my %hash;
open(IN,$refseq) or die $!;
<IN>;
while(<IN>){
    chomp;
    my ($bin,$name,$chr,$strand,$txstart,$txEnd,$cdsstart,$cdsend,$exoncount,$exonstarts,$exonends,$score,$name2,@a) = split "\t";
    $hash{$chr}->{$name} = [$txstart,$txEnd,$strand,$name2];
}
close IN;

print "[",scalar(localtime),"] Write out ..\n";
open(IN,$inbed) or die $!;
open(OUT,">$outfile") or dir $!;

while(<IN>){
    chomp;
    my ($chr,$start,$end) = split "\t";
    my $center = int(($start + $end)/2);

    my ($t,$c) = annotate($chr,$center,\%hash);
    print OUT join "\t",($chr,$start,$end,$t,$c);
    print OUT "\n";
       
}

sub annotate{
    my ($chr,$pos,$ref) = @_;
    my $type="Intergenic";
    my $comment="";
    my $flag=0;
    my %tmphash1;  # downstream
    my %tmphash2;  # upstream
     my %promoter;
        my %enhancer;
        my %genebody;
    foreach my $tid (sort keys %{$ref->{$chr}}){
        my ($start,$end,$strand,$gene) = @{$ref->{$chr}->{$tid}};
        my $tss=$start;
        $tss = $end if($strand eq "-");
        # If Promoter, then return
        my $dist = $pos-$tss;
        $dist = -1*$dist if ($strand eq "-");

       

        if($dist >= $promoterup && $dist <= $promoterdown){
            $type = "Promoter";
            $comment="$tid|$gene|$dist";
            $flag=1;
            $promoter{$dist} = [$tid,$gene];

        }elsif($dist>=$enhancerup && $dist<=$enhancerdown){
            $type="Enhancer";
            $comment="$tid|$gene|$dist";
            $flag=1;
            $enhancer{$dist}=[$tid,$gene];
        }elsif($pos>=$start && $pos <=$end){
            $type="GeneBody";
            $comment="$tid|$gene|$dist";
            $flag=1;
            $genebody{$dist}=[$tid,$gene];
        }else{
            if($dist>=0){
             $tmphash1{$dist}=[$tid,$gene];
            }else{
                $tmphash2{$dist}=[$tid,$gene];
            }
        }
    }

    if(!$flag){
        #Try to get the neighbor genes;
        # get downstream;
        my @down=sort {$a <=> $b } keys %tmphash1;
        my @up=sort {$b<=>$a} keys %tmphash2;

        my $up="NA";
        my $down="NA";
        if(scalar(@up)>0){
            $up=join "|",(@{$tmphash2{$up[0]}},$up[0]);
        }
        if(scalar(@down)>0){
            $down=join "|",(@{$tmphash1{$down[0]}},$down[0]);
        }
        $comment="$up;$down";
    }else{
        if(keys %promoter){
            my @tmp=sort {abs($a) <=>abs($b)} keys %promoter;
            $type="Promoter";
            $comment = join "|",(@{$promoter{$tmp[0]}},$tmp[0]);
        }elsif(keys %enhancer){
            my @tmp=sort {abs($a) <=>abs($b)} keys %enhancer;
            $type="Enhancer";
            $comment = join "|",(@{$enhancer{$tmp[0]}},$tmp[0]);
        }else{
            my @tmp=sort {abs($a) <=>abs($b)} keys %genebody;
            $type="GeneBody";
            $comment = join "|",(@{$genebody{$tmp[0]}},$tmp[0]);
        }
    }

    return ($type,$comment);
}
