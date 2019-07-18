#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;

my $input;
my $outDir = "./";

GetOptions('i=s'              => \$input,
           'o=s'              => \$outDir,
        );

unless ($input){
        print <<EOF;

        Usage :  perl process.vep.output.pl -i input_dir -o output
                        -i : input directory containing the MAF output files
                        -o : outDir - [ optional ], default is './'

EOF
        exit;
}

if (! -e "$input"){
	exit;
}


open(OUT,">$outDir/processed_MAF.tsv") || die;

### A,C, G,T 

my $start = 0;
my $chr_old;
my $gene_old;
my $protein_old;
my $context_old;
my %saw;



open(IN,"$input")||die;
while(<IN>){
	next if /^#/;
	chomp $_;
	
	my ($chr,$pos,$alt,$context, $cons, $gene, $protein) = (split /\t/, $_)[ 5,6,13,0,9,48,54];
	# $chr =~ s/^chr//;

    print STDOUT "$cons\t" ;
    my $a = &effect($cons);
    print STDOUT "$a\n";



	$cons = &effect($cons);
	if ($cons eq "NA"){
		next;
	}	
	if ($pos ne $start){
		if ($start != 0){
            print OUT join("\t",$chr_old,$start, $saw{"A"}||'0', $saw{"C"} ||'0', $saw{"G"} || '0', $saw{"T"} || '0', $gene_old, $protein_old,$context_old),"\n";
        }
        undef %saw;
        $saw{$alt} = $cons;
    }else{
        $saw{$alt} = $cons;
    }
	$start = $pos;
    $chr_old = $chr;
    $gene_old = $gene;
    $protein_old = $protein;
    $context_old = $context;

	if (eof(IN)){
        print OUT join("\t",$chr_old,$start, $saw{"A"}||'0', $saw{"C"} ||'0', $saw{"G"} || '0', $saw{"T"} || '0', $gene_old, $protein_old,$context_old),"\n";
	}
}
close IN;





sub effect{
	my ($effect) = @_;
    if ($effect =~ /Silent/ ){
        return "0";
    }elsif($effect =~ m(Missense_Mutation|In_Frame_Del |Nonsense_Mutation|In_Frame_Ins |Frame_Shift_Del |Frame_Shift_Ins)){
        return "1"
    }else{
        return "NA"
    }
}
    


#  0 - silent
#  1 - nonsilent



