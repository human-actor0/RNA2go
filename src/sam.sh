CIGAR_FORMAT='
Op BAM Description
M 0 alignment match (can be a sequence match or mismatch)
I 1 insertion to the reference
D 2 deletion from the reference
N 3 skipped region from the reference
S 4 soft clipping (clipped sequences present in SEQ)
H 5 hard clipping (clipped sequences NOT present in SEQ)
P 6 padding (silent deletion from padded reference)
= 7 sequence match
X 8 sequence mismatch

• H can only be present as the first and/or last operation.
• S may only have H operations between them and the ends of the CIGAR string.
• For mRNA-to-genome alignment, an N operation represents an intron. For other types of alignments,
the interpretation of N is not defined.
• Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
'

REF="
http://davetang.org/wiki/tiki-index.php?page=SAM
"

sam.variants(){
usage="$FUNCNAME <sam>
"
#1. calculate insertion and deletion points for each read
#2. collect coverage and other mutation for the indel sites (>10 reads maybe?).
#3. calculate significance of the indel by comparing mutation and control samples.
	if [ $# -lt 1 ];then echo "$usage"; return; fi
	cat $1 | perl -e 'use strict; my %res=();
	while(<STDIN>){ chomp; 
		if($_=~/^@/){ next;}
		my @a=split/\t/,$_;
		my $id=$a[0];
		my $flag=$a[1];
		my $chrom=$a[2];
		my $start=$a[3]-1;
		my $mapq=$a[4]; # -10log10 Pr( wrong )
		my $cigar=$a[5];
		my $seq=$a[9];
		my $len=0;
		#my $strand="+"; if ( $flag & 16 ){ $strand="-"; }
		my $strand="+"; if ( $flag & (0x10) ){ $strand="-"; }
		my $gpos=0; my $spos=0;
		while($cigar=~/(\d+)([MIDNSHPX=])/g){ 
			my ($n,$c)=($1,$2);
			if($c =~ /[MX=]/){
				if($c eq "X"){
					#print join("\t",( $chrom, $start+$gpos, $start+$gpos +$n, "X",substr($seq,$spos,$n),$strand)),"\n";
					for( my $i=$start+$gpos;  $i < $start+$gpos+$n; $i++){
						$res{ $chrom."\t".$i }{X}++;
					}
				}
				$gpos += $n;
				$spos += $n;
			}elsif($c =~ /[DN]/){
				#print join("\t",( $chrom, $start+$gpos, $start+$gpos +$n, "D",$n,$strand)),"\n";
				for( my $i=$start+$gpos;  $i < $start+$gpos+$n; $i++){
					$res{ $chrom."\t".$i }{D}++;
				}
				$gpos += $n;
			}elsif($c =~/I/){
				#print join("\t",( $chrom, $start+$gpos, $start+$gpos +$n, "I",substr($seq,$spos,$n),$strand)),"\n";
				$res{ $chrom."\t".($start+$gpos) }{I}++;
				$spos += $n;
			}elsif($c =~/S/){
				$spos += $n;	
			}else{ # P
			}
		}
	}
	print join("\t",("chrom","pos","X","D","I")),"\n";
	foreach my $k (keys %res){
		print $k;
		foreach my $t ( ("X","D","I")){
			my $v= defined $res{$k}{$t} ? $res{$k}{$t} : 0;
			print "\t$v";
		}
		print "\n";
	}
	'
}

sam.variants.test(){
## example from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/figure/F1/
echo \
"@HD	VN:1.0	SO:coordinate
@SQ	SN:chr1	LN:249250621
r1	163	chr1	7	30	8M2I4X1D3M	=	37	39	TTAGATAAAGGATACTG *
r2	0	chr1	9	30	3S6M1P1I4M	*	0	0	AAAAGATAAGGATA	*
r3	0	chr1	9	30	5H6X1S	*	0	0	AGCTAAa	*	NM:i:1
r4	0	chr1	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*
r3	16	chr1	29	30	6H5M	*	0	0	TAGGC	*	NM:i:0
r1	83	chr1	37	30	9M	=	7	-39	CAGCGCCAT	*"\
| sam.variants -
}


sam.bed12(){
usage="$FUNCNAME <sam> [-x]
"
# reference from https://samtools.github.io/hts-specs/SAMv1.pdf
if [ $# -lt 1 ]; then echo "$usage"; return; fi
	cat $1 | perl -ne 'chomp; my @a=split/\t/,$_; my $print_seq="'${2:-}'";
		if($_=~/^@/){ next;}
		my $id=$a[0];
		my $flag=$a[1];
		my $chrom=$a[2];
		next if $chrom eq "*";
		my $start=$a[3]-1;
		my $mapq=$a[4]; # -10log10 Pr( wrong )
		my $cigar=$a[5];
		my $seq=$a[9];
		my $len=0;
		my $strand="+"; if ( $flag & (0x10) ){ $strand="-"; }

		my $gseq=""; # genomic sequence 
		#\*|([0-9]+[MIDNSHPX=])+ 
		my @starts=(); 
		my @sizes=(); 
		my @seqs=();

		my $gpos=0; 
		my $spos=0; ## sequence offset
		my $prev_c="";
		while($cigar=~/(\d+)([MIDNSHPX=])/g){ 
			my ($x,$c)=($1,$2);
			if($c=~/[MX=]/){
				if($prev_c eq "D"){
					$sizes[$#sizes] += $x;
					$seqs[$#seqs] .= substr($seq,$spos,$x);
				}else{
					push @starts, $gpos; 
					push @sizes, $x;
					push @seqs, substr($seq,$spos,$x);
				}
				$spos += $x;
				$gpos += $x;
			}elsif($c=~/[DN]/){
				if( $c eq "D" ){
					$seqs[$#seqs] .= "-"x$x;
					$sizes[$#sizes] += $x;
				}
				$gpos += $x; 
			}elsif($c=~/[SI]/){
			## aware that soft/hard clipping does not affect genomic coordinates 
				$spos += $x;
			}elsif($c=~/P/){
				$seqs[$#seqs] .= "*"x$x;
				$sizes[$#sizes] += $x;
			}else{
				## P is not handled
			}
			$prev_c=$c;
		}
		my $end=$start+$starts[$#starts]+$sizes[$#sizes];
		my $str_sizes=join(",",@sizes);
		my $str_starts=join(",",@starts);
		my $str_seqs=join(",",@seqs);
		if($print_seq eq "-x"){
			$str_starts .= "\t".$str_seqs;
		}
		print join("\t",( 
			$chrom,$start,$end,$id,$mapq,$strand,
			$start,$end,"0,0,0",scalar @starts,
			$str_sizes,$str_starts
		)),"\n";
			
	'
}
sam.bed12.test(){
## example from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/figure/F1/
echo \
"@HD	VN:1.0	SO:coordinate
@SQ	SN:chr1	LN:249250621
r1	163	chr1	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG *
r2	0	chr1	9	30	3S6M1P1I4M	*	0	0	AAAAGATAAGGATA	*
r3	0	chr1	9	30	5H6M1S	*	0	0	AGCTAAa	*	NM:i:1
r4	0	chr1	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*
r3	16	chr1	29	30	6H5M	*	0	0	TAGGC	*	NM:i:0
r1	83	chr1	37	30	9M	=	7	-39	CAGCGCCAT	*"\
> tmp.inp
echo "## sam input";
cat tmp.inp
echo "## output:";
sam.bed12 tmp.inp -x
echo \
"chr1	6	22	TTAGATAAGATA*CTG	30	+
chr1	8	18	AGATAAGATA	30	+
chr1	8	14	AGCTAA	30	+
chr1	15	26	ATAGCT..............TCAGC	30	+
chr1	28	33	TAGGC	30	-
chr1	36	45	CAGCGCCAT	30	-" > tmp.exp

rm -rf tmp.exp tmp.obs tmp.inp
}
