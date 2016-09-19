
seq.dist(){
cat $1 | perl -e 'use strict;
	sub dist{
		my ($x, $y) = @_; my $n=0;
		for(my $i=0; $i<length($x);$i++){
			if( substr($x,$i,1) ne substr($y,$i,1)){ $n++; }
		}
		return $n;
	}
	my %S=();
	while(<STDIN>){ chomp;$_=~s/\s+/\t/g;  my ($id,$seq) = split/\t/,$_; $S{$id}=$seq; }
	foreach my $id1 (keys %S){
	foreach my $id2 (keys %S){ if( $id1 ne $id2){
		print $id1,"\t",$id2,"\t",$S{$id1},"\t",$S{$id2},"\t",dist($S{$id1},$S{$id2}),"\n";
	}}}
' | sort 
}
seq.dist.test(){
echo \
"1  GATCTG
2	GATCAG
3	GTAGCC
4	GTGGCC" | seq.dist -

}
seq.fasta(){
usage="$FUNCNAME <genome.fa> <bed>  [options]
 [options] : read \"bedtools getfasta\"
"
if [ $# -lt 2 ];then echo "$usage"; return; fi
#n1::one:0-3(+)
	bedtools getfasta -fi $1 -bed $2 ${@:3} -name -tab \
	| perl -ne 'chomp; my ($i,$s)=split/\t/,$_;
		$i=~/(\w+)::(\w+):(\d+)-(\d+)\(([\+|-])\)/;
		print $2,"\t",$3,"\t",$4,"\t",$1,"\t0\t",$5,"\t",$s,"\n";
	'
}

seq.read(){
usage="
USAGE: $FUNCNAME <fa> <bed> [options]
 [options]: 
	-s : reverse complements for the negative strand
	-i : ignore case 
	-uc : upper case

"
local option=${@:3};
if [ $# -lt 2 ]; then echo "$usage"; return; fi
	local tmpd=`hm util mktempd`;
	cat $2 > $tmpd/a;	
	cat $1 | perl -e 'use strict; 
		sub revc{
			my ($seq) = @_;
			$seq =~ tr/ACGTacgt/TGCAtgca/;
			return join("",reverse(split //,$seq));
		}

		my $file="'$tmpd/a'"; my $opt="'"$option"'";
		my %ref=(); my $chrom="";
		while(<STDIN>){ chomp;
			if($_=~/>(\S+)/){ $chrom=$1; next; }		
			$ref{$chrom} .= $_;
		}
		#foreach my $c (keys %ref){ print $c," ",substr($ref{$c},0,10),"\n"; }
		open(my $fh, "<", $file) or die "$!";
		while(<$fh>){chomp;
			my @a=split/\t/,$_;
			my $seq="NULL";
			if(defined $ref{$a[0]}){
				$seq=substr($ref{$a[0]},$a[1],$a[2]-$a[1]);
				if($opt =~ /s/ && $#a >=5 && $a[5] eq "-"){
					$seq=revc($seq);
				}
			}
			if($opt=~/lc/){ $seq = lc $seq; }elsif ($opt=~/uc/){ $seq = uc $seq; }
			print join("\t",@a),"\t",$seq,"\n";
		}
		close($fh);
		
	'
	rm -rf $tmpd;
}

seq.fasta.test(){
echo \
">one 
ATGCATGCATGCATGCATGCATGCATGCAT 
GCATGCATGCATGCATGCATGCATGCATGC 
ATGCAT 
>two another chromosome 
ATGCATGCATGCAT 
GCATGCATGCATGC" > tmp.fa
samtools faidx tmp.fa
echo \
"one	0	3	n1	0	+
two	0	5	n2	0	-" > tmp.bed
seq.fasta tmp.fa tmp.bed -s 
rm tmp.*
}
