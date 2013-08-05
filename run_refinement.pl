#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
my $currentdir=getcwd;

#1) Markduplicate
#2) Realignment
#3) Recalibration

my $refseq="/home/jiangli/reference/hg19/hg19.fa";
my $gatk="/home/jiangli/bin/seq/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar";
my $markdup="/home/jiangli/bin/seq/picard-tools-1.92/MarkDuplicates.jar";
my $dbsnp="/home/jiangli/reference/dbsnp135/dbsnp_135.hg19.vcf";
my $java="/home/jiangli/bin/jre1.7.0_21/bin/java";

my $pbs= <<PBS;
#!/bin/bash
#Beginning of PBS bash script
#PBS -M jiangli\@stanford.edu
#Status/Progress Emails to be sent
#PBS -m bae
#Email generated at b)eginning, a)bort, and e)nd of jobs
#PBS -l mem=25000mb
#Total job memory required (specify how many megabytes)
#PBS -l walltime=48:00:00
#You must specify Wall Clock time (hh:mm:ss) [Maximum allowed 30 days = 720:00:00]
#PBS -l nodes=1:ppn=1
#run 8 threads on 1 node
#PBS -q normal
#PBS -V
cd /oasis/scratch/\$USER/temp_project
PBS


my @files=<*_sort.bam>;
foreach my $f (@files){
	my $id = "";
	if($f=~/(.*)\.bam/){
		$id=$1;
	}
	
	my $out = "run_refinement_${id}.sh";
	open(OUT,">$out") or die $!;
    #PBS header
    print OUT $pbs;
    print OUT "#PBS -N $id\n";
    print OUT "cd $currentdir\n";
    print OUT "source ~/.bashrc\n";
    
	#markdup
	print OUT "echo Running markdup at \`date\`\n";
	my $comm="$java -Xmx15g -Djava.io.tmpdir=../tmp -jar $markdup INPUT=$f OUTPUT=${id}_markdup.bam METRICS_FILE=${id}_Metric.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT TMP_DIR=../tmp\n";
	print OUT  $comm;
    print OUT "samtools index ${id}_markdup.bam\n";

	# RealignerTargetCreator
	print  OUT "\necho Running RealignerTargetCreator  \`date\`\n";
	$comm="$java -Xmx15g -jar $gatk -T RealignerTargetCreator -R $refseq -I ${id}_markdup.bam -o ${id}.intervals\n";
	print OUT  $comm;
	# Realignment
	print OUT "\n echo Running realignment \`date\`\n";
	$comm="$java -Xmx15g -jar $gatk -T IndelRealigner -R $refseq -I ${id}_markdup.bam -targetIntervals ${id}.intervals -o ${id}_markdup_realign.bam\n";
	print OUT  $comm;
	print OUT "samtools index ${id}_markdup_realign.bam\n";

	# BaseRecalibrator
	print OUT "\n echo Running BaseRecalibrator \`date\`\n";
	$comm="$java -Xmx15g -jar $gatk -T BaseRecalibrator -R $refseq -I ${id}_markdup_realign.bam -o ${id}_recal_data.grp -knownSites $dbsnp --default_platform Illumina\n";
	print OUT $comm; 
	print OUT "\necho Running PrintReads \`date\`\n";
	$comm="$java -Xmx15g -jar $gatk -T PrintReads -R $refseq -I ${id}_markdup_realign.bam -o ${id}_markdup_realign_recal.bam -BQSR ${id}_recal_data.grp \n";
	print OUT $comm;
	print OUT "samtools index ${id}_markdup_realign_recal.bam\n";
}

