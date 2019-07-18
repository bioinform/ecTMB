#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

my $exon_file = "/net/kodiak/volumes/river/shared/users/yaol12/software/GSMuta_data/GSMuta/ensembl_81_exon_pos.hg38.txt";
my $output = "./";
my $help = 0; 
my $fasta;

GetOptions('e=s'              => \$exon_file,
           'o=s'              => \$output,
           'f=s'              => \$fasta,
           'h'                => \$help,  
	);


unless ($help == 0){
        print <<EOF;

        Usage :  perl annotate.exome.changes.pl -e exon_file -f fasta -o output_dir
                        -e : exon_file (3 fields, tab delimited; chr, start, end) 
                        -f : path to fasta file
                        -o : [ optional ] output directory, default is './'

EOF
        exit;
}


if ( !-d $output){
	`mkdir $output`;
}

my $outbed="${output}/extended.bed";
my $outvcf="${output}/output.vcf";


open(OUT, ">$outvcf")||die;
print "-- Generating VCF ... \n";

`awk '{OFS="\t"} {print \$1,\$2-3,\$3+2}' $exon_file  > $outbed`;
my @line = split /\n+/,`awk '{OFS="\t"}{print \$1,\$2-1,\$3+1}' $outbed | fastaFromBed -fi $fasta -bed stdin -fo stdout`;

print "Retrieve reference is done ... \n";

my %ref;
my $start;
my $end; 
foreach my $line (@line){
  if ($line =~ />/){
    ($start, $end) = (split />|:|-/,$line)[2..3];
  }else{
    chomp $line;
    my @string = split //,$line;
    for my $j (0 ..$#string){
      $ref{($j+$start)} = $string[$j];
    }
  }
    
}

open(A, $outbed)||die;
while(<A>){
  my ($chr,$start,$end) = split /\s+/,$_;
  for my $k($start ..($end-1)){
    my $ref_allele = uc($ref{$k});
    my $context = join("", uc($ref{($k-1)}), $ref_allele, uc($ref{($k+1)}));
    $context =~ tr/ACGT/1432/;
    
    if ($ref_allele ne "A"){
      print OUT $chr,"\t",$k+1,"\t", ".\t", $ref_allele,"\tA\t.\t.\t$context\n";
    }
    if ($ref_allele ne "C"){
      print OUT $chr,"\t",$k+1,"\t", ".\t", $ref_allele,"\tC\t.\t.\t$context\n";
    }
    if ($ref_allele ne "G"){
      print OUT $chr,"\t",$k+1,"\t", ".\t", $ref_allele,"\tG\t.\t.\t$context\n";
    }
    if ($ref_allele ne "T"){
      print OUT $chr,"\t",$k+1,"\t", ".\t", $ref_allele,"\tT\t.\t.\t$context\n";
    }
  }
}
undef %ref;
undef @line;


  