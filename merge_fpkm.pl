#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email: riverlee2008@gmail.com
# Date: Sun May 12 13:43:44 2013
###################################
use strict;
use warnings;

my $usage="perl $0 dir output\n"; 

my ($dir,$output) = @ARGV;
die $usage if(@ARGV!=2 || !-d $dir || $output eq "");

chdir $dir or die $!;

opendir (DIR,".") or die $!;
my @folders=grep {$_ ne "." && $_ ne ".." && -d $_} readdir DIR;
#print join "\n",@folders;
my %fpkm;
my %samples;
foreach my $f (@folders){
    print "Reading $f ... Is cufflinks ouput ";
    if(-e "$f/genes.fpkm_tracking"){
        print " Yes\n";
        $samples{$f}++;
        open(IN,"$f/genes.fpkm_tracking") or die $!;
        <IN>;
        while(<IN>){
            s/\r|\n//g;
            my($trackiggid,$classcode,$nearestrefid,$geneid,$geneshortname,$tssid,$locus,$length,$coverage,$fpkm,$fpkmlo,$fpkmhi,$fpkmstatus) = split "\t";
            $fpkm{$geneid}->{'samples'}->{$f}=$fpkm;
            $fpkm{$geneid}->{'locus'}=$locus;
        }
        close IN;
    }else{
        print " No\n";
    }
}

# write out
open(OUT,">$output") or die $!;
my @samples=sort keys %samples;
print OUT join "\t",("Gene","Locus",@samples);
print OUT "\n";
foreach my $g (sort keys %fpkm){
    print OUT join "\t",($g,$fpkm{$g}->{'locus'});
    foreach my $s (@samples){
        print OUT "\t",$fpkm{$g}->{'samples'}->{$s};
    }
    print OUT "\n";
}
close OUT;
