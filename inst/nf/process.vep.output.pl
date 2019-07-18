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
                        -i : input directory containing the VEP output files
                        -o : outDir - [ optional ], default is './'

EOF
        exit;
}

if (! -e "$input"){
	exit;
}


open(OUT,">$outDir/processed_vep.tsv") || die;

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
	
	my ($chr,$pos,$alt,$context, $cons, $gene, $protein) = (split /\t|;|\|/, $_)[0,1,4,7,9,12,22];
	# $chr =~ s/^chr//;

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
       	return "5" if ($effect =~ /splice_acceptor_variant|splice_donor_variant|transcript_ablation/);    		  
	return "2" if( $effect =~ /stop_gained/ );
    	return "3" if( $effect =~ /stop_lost/ );
    	return "4" if( $effect =~ /initiator_codon_variant|start_lost/ );
    	return "1" if( $effect =~ /missense_variant|conservative_missense_variant|rare_amino_acid_variant/ );
	return "NA" if ($effect =~ /transcript_amplification/);
    	return "0" if( $effect =~ /synonymous_variant|stop_retained_variant/ );
    	return "NA" if ( $effect =~ /splice_region_variant/ );
	return "0" if ($effect =~ /incomplete_terminal_codon_variant/);
	return "1" if ($effect =~ /protein_altering_variant|coding_sequence_variant/);
	return "NA";		
}


#  0 - silent
#  1 - miss-sense
#  2 - nonsense
#  3 - nonstop
#  4 - TSS
#  5 - splice 
