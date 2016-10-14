
bed.sort(){
	sort -k1,1 -k2,3n
}
bed.toy(){
usage="$FUNCNAME <intput.txt>"
if [ $# -lt 1 ];then echo "$usage"; return; fi
cat $1 | perl -e 'use strict; 
	sub f{
		my ($s,$e,$g) = @_;
		my $n=scalar @$s;
		#print join(",",@$s)," ",join(",",@$e),"\n";
		my $s0=$s->[0];		
		my $e0=$e->[$n-1];		
		my $st="+"; if( lc $g eq $g){ $st="-";}
		print "chr1\t$s0\t$e0\t$g\t0\t$st\t$s0\t$e0\t0,0,0\t$n\t";
		print join(",",map { $e->[$_] - $s->[$_]} (0..($n-1))),"\t";
		print join(",",map { $_ - $s0 } @$s),"\n";

	}
	while(<STDIN>){chomp;  my $tmp=" ".$_." "; $tmp=~s/\d/ /g;
		my @S=split//,$tmp;
		my @sizes=(); my @starts=(); my @ends=();
		my $type="";
		my $start=0;
		for(my $i=1; $i< scalar @S; $i++){ 
			my ($a,$b)=($S[$i-1],$S[$i]);
			if($a eq " " && $b =~/\w/){
				push @starts,$i-1;
			}elsif($a =~/\w/ && $b eq " "){
				push @ends,$i-1;
				f(\@starts,\@ends,$a);
				@starts=();@ends=();
			}elsif($a eq "-" && $b =~/\w/){
				push @starts,$i-1;
			}elsif($a =~/\w/ && $b eq "-"){
				push @ends,$i-1;
			}elsif($a ne $b){
				push @ends,$i-1;
				f(\@starts,\@ends,$a);
				@starts=();@ends=();
				push @starts,$i-1;
			} 
		}
	}
'
}

bed.toy.test(){
echo \
"
01234567890123456789012345678901234567890123456789
     AABBBCCCC
  DDDDD-----DDDDD EEE-EEE
  e-----e--------e
" | bed.toy -
}

bed.exon(){
usage="
USAGE: $FUNCNAME <bed12> 
"
if [ $# -ne 1 ];then echo "$usage"; return; fi
	##  [es   ]s----e[    ee]
	cat $1 | awk -v OFS="\t" '{
		## take introns
		split($11,sizes,",");
		split($12,starts,",");
		for(i=1;i<= $10;i++){
			## exon
			s = $2 + starts[i];
			e = s + sizes[i];
			print $1,s,e,$4,$5,$6;
		}	
	}' 
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
