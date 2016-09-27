
seq.kmers(){
## obtained from Cook malcolm
## http://comments.gmane.org/gmane.comp.lang.perl.bio.general/18242
perlcode="
        sub kmer{
                my ($prefix,$t,$l,$r) = @_;
                if($l<=0){ 
                        push @$r,$prefix;
                        return;
                }
                foreach my $nu ( "A", "C", "G", "T"){
                        $t->{$nu} = undef;
                        kmer($prefix.$nu,$t->{$nu},$l-1,$r);
                }       
        }       
        my %trie=();
        my @kmers=();
        kmer("",\%trie,3,\@kmers);
        print join("\n",@kmers),"\n";

"
k=$1;
s=$( printf "%${k}s" ); # a string with $k blanks
s=${s// /{A,T,G,C\}};   # substitute '{A,T,G,C}' for each of the k blanks
echo 'kmers using bash to expand:' $s > /dev/stderr
bash -c "echo  $s";     # let brace expanion of inferior bash compute the cross product

}

seq.mut(){
usage="$FUNCNAME <seq> <mutations>";
	echo $1 | perl -ne 'chomp; my @seq=split //,$_;
		my @b=("A","C","G","T");
		print $_;
		for(my $i=0; $i<=$#seq; $i++){
			foreach my $bi (@b){ if($bi ne  $seq[$i]){
				my @s=@seq;
				$s[$i]=$bi;
				print " ",join( "",@s);
			}}		
		}
		print "\n";
	'	
}


seq.mut.test(){
echo \
"AAAA CAAA GAAA TAAA ACAA AGAA ATAA AACA AAGA AATA AAAC AAAG AAAT" > exp
seq.mutate "AAAA" > obs
check obs exp
rm obs exp
}
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
seq.get(){
usage="$FUNCNAME <genome.fa> <bed>  [options]
 [options] : 
  -s : revCompment for the negative strand 
"
if [ $# -lt 2 ];then echo "$usage"; return; fi
#n1::one:0-3(+)
	bedtools getfasta -fi $1 -bed ${2/-stdin} ${@:3} -name  -tab \
	| perl -ne 'chomp; if( $_=~ /(.+)::(.+)\t([ACGT]+)/){
		print $1,"\t",$3,"\n";
		}
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

seq.get.test(){
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
seq.get tmp.fa tmp.bed -s 
rm tmp.*
}
