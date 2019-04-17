#!/usr/bin/perl -w

my $maf = $ARGV[0];
my $ref = $ARGV[1];
my $out = $ARGV[2];


my $chr;
my $id;
my $pos;
my %ref;
my @line = split /\n+/,`grep -v "Ensembl" $maf | awk '{OFS="\t"}{print \$2,\$3-2,\$3+1}' | fastaFromBed -fi $ref -bed stdin -fo stdout`;
foreach my $line (@line){
	if ($line =~ />/){
		($chr, $pos) = (split />|:|-/,$line)[1..2];
                $id = join("",$chr,":",$pos);
	}else{
		 chomp $line;
                 $ref{$id} =$line;
	}
}


open(O,">$out")||die;
open(IN,$maf)||die;
while(<IN>){
	if (/^Ensembl_gene_id/){
		chomp $_;
		s/\s+$//g;
		print O $_,"\t","Context\n";
	}else{
		chomp $_;
		if (/In_frame/ || /Frame_shift/){
			print O $_,"\tINDEL\n";
		}else{
			my @temp = split /\t/,$_;
			$id = join("",$temp[1],":",($temp[2]-2));
			print O $_,"\t",&convert(join("",$ref{$id},$temp[11])),"\n";
		}
	}
}

sub convert{
	my ($input) = @_;
	$input =~ tr/ATGC/1234/;
	return ($input);
}
