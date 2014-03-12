#!/usr/bin/perl
use strict;
use warnings;

my @files = <*.bw>;
foreach my $f (@files){
    my $n=$f;
    $n=~s/\.bw//g;
    $n=~s/_norm_sorted_rmdup//g;
    my $str="track type=bigWig name=\"$n\" description=\"$n\" color=25,25,225 visibility=full bigDataUrl=https://s3.amazonaws.com/ATAC_mouse/$f maxHeightPixels=32:32:32";

    print $str,"\n";
}

