#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email: riverlee2008@gmail.com
# Date: Sun Aug  4 22:35:29 2013
###################################
use strict;
use warnings;

# Read samtools flagstat output in a folder and then merge the result together
# Just take the total reads and mapped reads


my ($dir) = @ARGV;

die "$0 flagstatFolder \n" unless ($dir);

opendir(DIR,$dir) or die $!;

my @files = grep {$_ ne "../" && $_ ne "./" && -f "$dir/$_"} readdir DIR;
#my @files = readdir DIR;
print join "\t",("#Sample","Total","Mapped","Rate\n");
foreach my $f (@files){
#    print $f,"\n";
    my $ff="$dir/$f";
    open(IN,$ff) or next;
    my $totalline=<IN>;
    next if($totalline !~/QC\-passed reads/);
    <IN>;
    my $mappedline=<IN>;
    close IN;

    my $total=0;
    my $mapped=0;

    if($totalline=~/(^\d+)/){
        $total=$1;
    }
    if($mappedline=~/(^\d+)/){
        $mapped=$1;
    }

    print join "\t",($f,$total,$mapped,$mapped/$total);
    print "\n";
}




