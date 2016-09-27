
bed.toy(){
usage="$FUNCNAME <intput.txt>"
if [ $# -lt 1 ];then echo "$usage"; return; fi
cat $1 | perl -ne 'use strict; chomp; 
	my @S=split//," ".$_." "; ## add pads for eacy calc.
	my @sizes=(); my @starts=();
	my $type="";
	my $start=0;
	for(my $i=1; $i< scalar @S; $i++){ 
		my $a=$S[$i-1];
		my $b=$S[$i];
		my $p=$i-1; ## original position before the padding
		if( ($a eq " " || $a eq "-") && $b =~ /[\w]/){
			$type=$b;
			if( $a eq " "){ $start=$p;}
			push @starts,$p-$start;
		}elsif($a =~ /[\w]/ && ($b eq " " || $b eq "-")){
			push @sizes, $p - $starts[$#starts] - $start;
		}
	}
	if( $type ne ""){
		my $strand = $type eq uc $type ? "+" : "-";
		my $end = $start + $starts[$#starts] + $sizes[$#sizes];
		print join("\t",(
			"chr1",
			$start, $end,
			$type, 0, $strand,
			$start, $end,
			"0,0,0",scalar @starts,
			join(",",@sizes),join(",",@starts)
		)),"\n";
	}
'
}

bed.toy.test(){
echo \
"01234567890123456789012345678901234567890123456789
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
 EEEE-------------EEEEEEEE------------EEEEEEEEEEEE
  RRR-------------RR
   rr-------------rrrrrr
                RRRR 
                                    RRRRR
" | bed.toy -
}

bed.exon(){
usage="
Exons of transcripts are merged into a gene when their boundaries are equal.
Exons are sorted by their 5 prime occurrence.
Suffix indexes are added (e.g., genename.exon#.sub ).
USAGE: $FUNCNAME <bed12>

"
if [ $# -ne 1 ];then echo "$usage"; return; fi
	cat $1 | perl -e 'use strict; my %res=();
	while(<STDIN>){ chomp;my @a=split/\t/,$_;
		my @sizes=split/,/,$a[10];	
		my @starts=split/,/,$a[11];	
		my $n=$a[9];
		for( my $i=0; $i<$n; $i++){
			my $id=$a[0]."\t".$a[3]."\t".$a[5];
			my $s= $a[1] + $starts[$i];
			my $e= $s + $sizes[$i];
			if( $a[5] eq "+"){
				$res{$id}{$s}{$e} ++; 
			}else{
				$res{$id}{-$e}{-$s} ++; 
			}
		}
	}
	foreach my $id (keys %res){
		my ($chr,$gene,$strand) = split /\t/,$id;
			my $i=0;
			foreach my $s (sort {$a<=>$b} keys %{$res{$id}}){ my $j=0;
				foreach my $e (sort {$a<=>$b} keys %{$res{$id}{$s}}){ my $n="$gene.E$i.$j";
					if($strand eq "+"){
						print $chr,"\t",$s,"\t",$e,"\t",$n,"\t0\t$strand\n";	
					}else{
						print $chr,"\t",-$e,"\t",-$s,"\t",$n,"\t0\t$strand\n";	
					}
					$j++;
				}
				$i++;
			}
		#if($strand eq "+"){ my $i=0;
		#	foreach my $s (sort {$a<=>$b} keys %{$res{$id}}){ my $j=0;
		#		foreach my $e (sort {$a<=>$b} keys %{$res{$id}{$s}}){ my $n="$gene.E$i.$j";
		#			print $chr,"\t",$s,"\t",$e,"\t",$n,"\t0\t$strand\n";	
		#			$j++;
		#		}
		#		$i++;
		#	}
		#}else{ my $i=0;
		#	foreach my $e (sort {$b<=>$a} keys %{$res{$id}}){ my $j=0;
		#		foreach my $s (sort {$b<=>$a} keys %{$res{$id}{$e}}){ my $n="$gene.E$i.$j";
		#			print $chr,"\t",$s,"\t",$e,"\t",$n,"\t0\t$strand\n";	
		#			$j++;
		#		}
		#		$i++;
		#	}
		#}
	}
		
	'
}
bed.exon.test(){
echo "
012345678901234567890123456789
 EEE---------EEE----EEEEEE
 EEE---------EEEE---EEEEEE
 eee---------eee----eeeeee
 eee---------eeee---eeeeee
" |  bed.toy - | bed.exon - 
}
bed.intron(){
usage="
USAGE: $FUNCNAME <bed12>
"
if [ $# -ne 1 ];then echo "$usage"; return; fi
	##  [es   ]s----e[    ee]
	cat $1 | awk -v OFS="\t" '$10 > 1{
		## take introns
		split($11,sizes,",");
		split($12,starts,",");
		for(i=1;i< $10;i++){
			## intron
			s = $2 + starts[i]+sizes[i];	
			e = $2 + starts[i+1];
			#ls = $2 + starts[i-1]; le = ls + sizes[i-1]; rs = $2 + starts[i]; re = rs + sizes[i];
			print $1,s,e,$4,$5,$6;
		}	
	}' 
}

bed.intron.test(){
echo "
012345678901234567890123456789
 EEE----------EEEE--EEEEEE
" |  bed.toy - | bed.intron - 
}

bed.rmchr(){
	cat $1 | perl -ne '$_=~s/^chr//;print $_;'
}
bed.enc(){
## encode bed6 in the name field 
	cat $1 | awk -v OFS="\t" -v j=${2:-"@"} '{ n=$1;for(i=2;i<= 6;i++){ n=n j $(i);} $4=n; }1' 
}

bed.split(){
usage="$FUNCNAMME <bed> <column> <output>"
if [ $# -ne 3 ];then echo "$usage"; return; fi
	if [ -z "${3##*\/*}" ];then
		mkdir -p ${3%/*}
	fi
        awk -v OFS="\t" -v C=$2 -v O=$3 '{
		fout=O"."$(C);
                print $0 >> fout;
        }' $1;
}

bed.split.test(){
echo \
"chr1	1	2
chr2	3	4" | bed.split - 1 tmp
rm -rf tmp
}

bed.3p(){
 	awk -v OFS="\t" '{ if($6=="-"){$3=$2+1;} $2=$3-1; print $0; }' $1
}

bed.5p(){
 	awk -v OFS="\t" '{if($6=="-"){$2=$3-1;}$3=$2+1; print $0; }' $1
}

bed.flank(){
usage="
FUNCT : extract flanking regions
USAGE : $FUNCNAME <bed> <left> <right> [<strand_opt>]
	<strand_opt> :  -s: strand specific
" 
if [ $# -lt 3 ];then echo "$usage"; return; fi
	awk -v OFS="\t" -v L=$2 -v R=$3 -v S=$4 '{ 
		if(S == "-s"  && $6 == "-"){
			$2=$2-R; $3=$3+L;
		}else{ $2=$2-L; $3=$3+R; } 
		if( $2 < 0){ $2=0;} if( $3 < 0){ $3=1;}
		print $0;
	}' $1;
}

bed.sum(){
usage="
$FUNCNAME <input.bed> [<input.bed>] > <out.bed>
"
if [ $# -lt 1 ];then echo "$usage"; return; fi
	perl -e 'use strict; my %res=();
	foreach my $f (@ARGV){
		my $fh; open($fh, $f) or die "$!";
		while(<$fh>){ chomp; my @a=split/\t/,$_;
			$res{ join("@",(@a[0..3],$a[5])) } += $a[4]; 
		}
		close($fh);
	}
	foreach my $k (keys %res){ 
		my ($c,$s,$e,$n,$st) = split /@/,$k;
		print $c,"\t",$s,"\t",$e,"\t",$n,"\t",$res{$k},"\t",$st,"\n";
	}
	' $@; 
}

bed.sum.test(){
echo \
"chr1	1	2	a	1	+
chr1	1	2	a	2	-
chr1	1	2	a	3	+" | bed.sum - > tmp.obs

echo \
"chr1	1	2	a	2	-
chr1	1	2	a	4	+" > tmp.exp
hm util check tmp.exp tmp.obs
rm tmp.obs tmp.exp
}
