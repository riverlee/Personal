#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email: riverlee2008@gmail.com
# Date: Wed Aug  7 16:35:41 2013
###################################
use strict;
use warnings;

# Take a fastq file which containing mulitple sequences.
# Write them in individual one in fasta format

my ($in) = @ARGV;

die "Usage: fa2individualfa.pl mm9.fa\n" if (@ARGV<1 || ! -e $in);
my $header="";
my $seq="";
open(IN,"$in") or die $!;
while(<IN>){
    if(/^>/){
        print $_;
        s/\r|\n|>//g;
        if($header eq ""){
            $header=$_;
        }else{
            writeout($header,$seq);
            $header=$_;
        }
        $seq="";
    }else{
        $seq.=$_;
    }
}

writeout($header,$seq);

sub writeout{
    my($h,$seq) = @_;
    open(OUT,">${h}.fa") or die $!;
    print OUT ">$h\n";
    print OUT $seq;
    close OUT;
}
