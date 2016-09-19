
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
USAGE : $FUNCNAME <bed> <left> <right> <strand_opt>
	<strand_opt> :  -s: strand specific
" 
if [ $# -ne 4 ];then echo "$usage"; return; fi
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
