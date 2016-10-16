
splicing.list_events(){
echo "
          _____a_____ 
        __b__     __c__
    ------[         ]--------
         d/         \e  
    \          f         /



## coverage :
          >>>>> (5) 
     >>>>>>>>(3) >>>>--------(4)
    ------[         ]--------
"
}
splicing.coverage(){
usage="
$FUNCNAME <flatexon.bed> <read.bed12> [options]
"
if [ $# -lt 2 ];then echo "$usage";return; fi
	intersectBed -a $1 -b $2 -wa -wb ${@:3} \
	| perl -e 'use strict;
	sub min{ my ($x,$y)=@_; return ($x > $y ? $y : $x);}
	sub max{ my ($x,$y)=@_; return ($x < $y ? $y : $x);}
	my %res=();
	while(<STDIN>){ chomp; my @a=split/\t/,$_;
		my @l=split/,/,$a[$#a-1];
		my @s=split/,/,$a[$#a];
		my $v=0;
		for(my $i=0; $i< $a[$#a-2];$i++){
			my $start=$a[7]+$s[$i];
			my $end=$start + $l[$i];
			my $ms=max($a[1],$start);
			my $me=min($a[2],$end);
			if( $me - $ms > 0){ $v += $me-$ms; }
		}
		if( $v > 0 ){
		my $k=join("\t",@a[0..5]);
		$res{$k} += $v;
		}
	}
	foreach my $k (keys %res){
		my @a=split /\t/,$k;
		$a[4]=$res{$k};
		print join("\t",@a),"\n";
	}
	'
}

splicing.coverage.test(){
echo "
012345678901234567890123456789
        AAABBBCCC
" | hm bed toy - | cut -f 1-6 > tmp.e
echo "
     RRRRRRRR
         RR---RR
" | hm bed toy - > tmp.r
splicing.coverage tmp.e tmp.r
rm tmp.*
}

splicing.flatexon(){
usage="
$FUNCNAME <gene.bed12> [options]
 [options]:
	-i : include introns
"
if [ $# -lt 1 ]; then echo "$usage";exit; fi
	cat $1 | hm bed exon - | perl -e 'use strict; my %res=();
	while(<STDIN>){chomp; my @a=split/\t/,$_;
		my $k=$a[0]."\t".$a[3]."\t".$a[5];
		$res{$k}{ $a[1] }{1} ++;
		$res{$k}{ $a[2]-1 }{2} ++;
	}
	foreach my $k (keys %res){
		my ($chrom,$gene,$strand)= split/\t/,$k;
		my @x= sort {$a<=>$b} keys %{$res{$k}};
		my %h= ();
		my %se=();
		foreach my $xi (@x){
			foreach my $xj (keys %{$res{$k}{$xi}}){
				$h{$xi} += $xj;
			}
		}
		for(my $i=0; $i<$#x; $i++){
			my ($s,$e)=( $x[$i], $x[$i+1]);
			my ($a,$b)=($h{$s},$h{$e});
			my $t="E";
			if( $a == 3){ $se{$s."\t".($s+1)."\t$gene\t$t"}=1; }
			if( $b == 3){ $se{$e."\t".($e+1)."\t$gene\t$t"}=1; }
			if( $a > 1){ $s++;}
			if( $b== 2){ $e++;}
			if( $a > 1  && $b !=2 ){ $t="I";}
			$se{$s."\t".$e."\t$gene\t$t"}=1;
		}
		foreach my $x (keys %se){
			print join("\t",($chrom,$x,$strand)),"\n";
		}
	}
	' | hm bed sort -
}
splicing.flatexon.test(){
echo \
"
0123456789012345678901234567890123456789
   EEEEEE-----EEEEEEE----EE-E  EEEEE
        EEEEEEE---------EEE-E  EE  
"> tmp.i 
cat tmp.i
cat tmp.i |  tail -n+3 | hm bed toy - \
| splicing.flatexon - -i

rm tmp.*
}
splicing.a2f_to_table(){
usage="
$FUNCNAME <cmd> <a2f> [<a2f>]
 <cmd> :
	a : count a per entry
	
"
if [ $# -lt 2 ];then echo "$usage"; return; fi
cmd='use strict;
	my $cmd="'$1'"; my $tmp="'${@:2}'"; 
	my $L=50;
	my @files=split/\s+/,$tmp;
	my %res=();
	my $i=0; my $nc=0;
	foreach my $f (@files){
		open(my $fh, "<", $f) or die  "$! $f";
		while(<$fh>){chomp; my @a=split/\t/,$_;
			my @s=split/,/,$a[4]; $a[4]=0; 
			my $k=join("\t",@a[0..5]);
			if( defined $res{$k}{$i}){
				for(my $j=0; $j<=$#s; $j++){
					$res{$k}{$i}[$j] += $s[$j];
				}
			}else{
				$res{ $k }{$i}=\@s;
				$nc=scalar @s;
			}
		}
		close($fh);
		$i++;
	}
	
	## print header
	my $n=$i;
	my $null = join(",", (0) x $nc);
	for(my $i=0; $i< $n; $i++){
		print "#f",$i,": ",$files[$i],"\n";
	}
	print "chrom\tstart\tend\tgene\tscore\tstrand";
	for(my $j=0; $j<$n; $j++){
		print "\tf$j.count";
	}
	print "\n";

	foreach my $k (keys %res){
		my @a=split/\t/,$k;
		print $k;
		for(my $i=0;$i<$n; $i++){
			my $v= defined $res{$k}{$i}[0] ? $res{$k}{$i}[0] : 0;
			print "\t$v";
		}
		print "\n";
	}
'
echo "$cmd" | perl 

}

splicing.a2f_to_table.test(){
echo \
"chr1	1	2	n1	1,2,3,4,5,6	+
chr1	3	5	n1	1,2,3,4,5,6	+
chr1	1	2	n2	1,2,3,4,5,6	+" > tmp.a 
echo \
"chr1	1	2	n1	1,2,3,4,5,6	+
chr1	3	4	n3	1,2,3,4,5,6	+" > tmp.b
splicing.a2f_to_table a tmp.a tmp.b
rm tmp.*
}

splicing.psi(){

usage="$FUNCNAME <read_len> <a2j> [<a2j> ..]"
if [ $# -lt 2 ];then echo "$usage";return; fi
cmd='use strict;
	my $L="'$1'";
	my $tmp="'${@:2}'"; 
	my @files=split/\s+/,$tmp;
	my %res=();
	my $i=0;
	my $nc=0;
	sub psi{
		my (@a)=@_;
		my $num=$a[0] + $a[3] + $a[4];
		my $den=$a[5];
		my $res="NaN";
		if( $num + $den > 0){ 
			$res=sprintf("%.4f",$num/($den + $num));
		}
		return( int($num+$den+0.5)."\t".$res) ;
	}
	foreach my $f (@files){
		open(my $fh, "<", $f) or die  "$! $f";
		while(<$fh>){chomp; my @a=split/\t/,$_;
			my @s=split/,/,$a[4]; $a[4]=0; 
			my $k=join("\t",@a[0..5]);
			if( defined $res{$k}{$i}){
				for(my $j=0; $j<=$#s; $j++){
					$res{$k}{$i}[$j] += $s[$j];
				}
			}else{
				$res{ $k }{$i}=\@s;
				$nc=scalar @s;
			}
			if($a[2]-$a[1] > $L){
				$res{$k}{$i}[0] *= $L/($a[2]-$a[1]-$L);
			}
		}
		close($fh);
		$i++;
	}
	my $n=$i;
	my $null = join(",", (0) x $nc);
	for(my $i=0; $i< $n; $i++){
		print "#f",$i,": ",$files[$i],"\n";
	}
	print "chrom\tstart\tend\tgene\tscore\tstrand";
	for(my $j=0; $j<$n; $j++){
		print "\tf$j.sup\tf$j.psi";
	}
	print "\n";
	foreach my $k (keys %res){
		my @a=split/\t/,$k;
		print $k;
		for(my $i=0;$i<$n; $i++){
			my $psi=psi(@{$res{$k}{$i}});
			print "\t$psi";

		}
		print "\n";
	}
'
echo "$cmd" | perl 
}
splicing.psi.test(){
echo \
"chr1	1	2	n1	1,2,3,4,5,6	+
chr1	3	5	n1	1,2,3,4,5,6	+
chr1	1	2	n2	1,2,3,4,5,6	+" > tmp.a 
echo \
"chr1	1	2	n1	1,2,3,4,5,6	+
chr1	3	4	n3	1,2,3,4,5,6	+" > tmp.b
echo \
"chr1	53699213	53699324	ENSG00000162385:E2.1	0,8,8,1,4,1	-
chr1	40030357	40030445	ENSG00000090621:E9.0	0,10,5,1,3,1	-
chr1	154940679	154940733	ENSG00000160691:E6.0	0,4,5,2,1,1	-
chr1	85724617	85724744	ENSG00000162642:E1.0	0,2,2,2,1,1	-
chr1	94654392	94654497	ENSG00000137962:E15.0	0,7,18,23,8,1	-" > tmp.c

splicing.psi 50 tmp.c
rm tmp.*
}
splicing.count_def(){
usage="
$FUNCNAME <exon.bed> <jc.bed> [options]
"
if [ $# -lt 2 ];then echo "$usage"; return; fi
	perl -e 'use strict; my $opt="'${@:3}'"; 
	my ($f1, $f2) = ("'$1'","'$2'");
	
	sub sw{
		my ($s)=@_;
		my $res=$s; $res=~tr/+-/-+/;
		return $res;
	}

	## read junctions
	my %J=();  my %JE=();
	open(my $fh,"<",$f2) or die "$!: $f2";
	while(<$fh>){chomp; my @a=split/\t/,$_;
		my ($l,$r)=($a[1],$a[2]-1);
		if(!( $opt=~/-S/) ){
			$J{$a[5]}{$a[0]}{$l}{$r} += $a[4];
			$J{$a[5]}{$a[0]}{$r}{$l} += $a[4];
		}
		if(!( $opt=~/-s/) ){
			$J{sw($a[5])}{$a[0]}{$l}{$r} += $a[4];
			$J{sw($a[5])}{$a[0]}{$r}{$l} += $a[4];
		}
	}
	close($fh);

	#foreach my $s (keys %J){
	#foreach my $c (keys %{$J{$s}}){
	#foreach my $l (keys %{$J{$s}{$c}}){
	#foreach my $r (keys %{$J{$s}{$c}{$l}}){
	#	print $s," ",$c," ",$l," ",$r," ",$J{$s}{$c}{$l}{$r},"\n";	
	#}}}}

	## read exons
	my %E=();
	open(my $fh,"<",$f1) or die "$!: $f1";
	while(<$fh>){chomp; my @a=split/\t/,$_;
		my $k=join("@",@a[0..5]);
		$E{$a[5]}{$a[0]}{l}{$a[1]}{$k}=$a[2]-1;
		$E{$a[5]}{$a[0]}{r}{$a[2]-1}{$k}=$a[1];
	}
	close($fh);

	##[    x]-----[l   r]----[y    ]
	my %E2=(); # exon-exon interaction
	foreach my $s (keys %E){ 
	foreach my $c (keys %{$E{$s}}){
	foreach my $l (keys %{$E{$s}{$c}{l}}){
	foreach my $e (keys %{$E{$s}{$c}{l}{$l}}){
		my $r=$E{$s}{$c}{l}{$l}{$e};
		my ($lc,$rc)=(0,0);
		my %xs=();
		my %ys=();
		if( defined $J{$s}{$c}{$l} ){ 
		foreach my $x (keys %{$J{$s}{$c}{$l}}){	
			if( defined $E{$s}{$c}{r}{$x} ){
				$lc+=$J{$s}{$c}{$x}{$l};
				$xs{$x}++;
			}
		}}
		if( defined $J{$s}{$c}{$r} ){ 
		foreach my $y (keys %{$J{$s}{$c}{$r}}){	
			if( defined $E{$s}{$c}{l}{$y} ){
				$rc+=$J{$s}{$c}{$r}{$y};
				$ys{$y}++;
			}
		}}
		my $f=0;
		foreach my $x (keys %xs){
		foreach my $y (keys %ys){
			$f += (defined $J{$s}{$c}{$x}{$y} ? $J{$s}{$c}{$x}{$y} : 0);
			
		}}
		my @b=split/@/,$e;
		$b[4] .= ",$lc,$rc,$f";
		print join ("\t",@b),"\n";
	}}}}

	'
}
splicing.count_def.test(){
echo \
"
0123456789012345678901234567890123456789
    JJJJ-------------JJJ
                           J--------JJJJJ
                          JJ--------JJJJJ
                          JJ-------JJJJJ
      jj----------------------------jj
                     AAAAAAA
      jj-------------jj
                          jj--------jj
    BBBB                             
                                    CCCC

" | hm bed toy - > tmp.i
tail -n+2 tmp.i | awk 'toupper($4)=="J"' | hm splicing count_jc - > tmp.j
tail -n+2 tmp.i | awk 'toupper($4)!="J"' | cut -f1-6 > tmp.e

echo "all";
splicing.count_def tmp.e tmp.j
echo "-s";
splicing.count_def tmp.e tmp.j -s
echo "-S";
splicing.count_def tmp.e tmp.j -S
rm tmp.*
}

splicing.exon(){
usage="
FUNCTION:
	Exons of transcripts are merged into a gene (4th column) when their boundaries are equal.
	Exons are sorted by their 5 prime occurrence.
	Suffix added (e.g., genename.exon#.sub ).
USAGE: 
	$FUNCNAME <bed12> 
"
if [ $# -lt 1 ];then echo "$usage"; return; fi
	hm bed exon $1 \
	| sort -u \
	| awk '{ print $1"@"$4"@"$6"\t"$2"\t"$3;}' \
	| hm bed sort \
	| mergeBed -i stdin -c 2,3 -o collapse,collapse \
	| perl -e 'use strict; 
	my %res=(); 
	while(<STDIN>){chomp; my @a=split/\t/,$_;
		my ($chrom,$name,$strand)=split /@/,$a[0];
		if($strand eq "+"){
			$res{$a[0]}{$a[1]}=$a[3]."\t".$a[4];
		}else{
			$res{$a[0]}{-$a[2]}=$a[3]."\t".$a[4];
		}
	}
	foreach my $k (keys %res){
		my $i=0;
		my ($chrom,$name,$strand)=split /@/,$k;
		foreach my $s (sort {$a<=>$b} keys %{$res{$k}}){
			my ($starts_str,$ends_str)=split/\t/,$res{$k}{$s};
			my @starts=split/,/,$starts_str;
			my @ends=split/,/,$ends_str;
			foreach my $j ( 0..$#starts){
				print $chrom,"\t",$starts[$j],"\t",$ends[$j],"\t",$name,":","E$i.$j\t0\t$strand\n";
			}
			$i++;
		}
	}
	
	'

}
splicing.exon.test(){
echo "
012345678901234567890123456789
 EEE---------EEE----EEEEEE
 EEE---------EEEE---EEEEEE
 EEE--------EEEE---EEEEEE
 eee---------eee----eeeeee
 eee---------eeee---eeeeee
" |  hm bed toy - | splicing.exon - tmp.out -s 
}
splicing.count_boundary(){
usage="
$FUNCNAME <bed> <bed> [options]
"
if [ $# -lt 2 ];then echo "$usage"; return; fi
	intersectBed -a ${1:-stdin} -b ${2:-stdin} -wa -wb ${@:3} \
	| awk -v OFS="@" '{
		l=0;r=0;
		if( $2 > $8 ){ l=1;}
		if( $3 < $9 ){ r=1;}
		print $1,$2,$3,$4,$5,$6"\t"l"\t"r;
	}' | hm stat sum -  | tr "@" "\t"
}

splicing.count_boundary.test(){
echo \
"01234567890123456789012345678901234567890123456789
        EEEEEE
   RRRRR
     RRRR
         RRRRR
             RRRRR
              RRRRR
         rrrrrr
" | hm bed toy - | awk '$4!=0' | cut -f1-6 > tmp.i
grep E tmp.i > tmp.e
grep -v E tmp.i > tmp.r
head tmp.*
	splicing.count_boundary tmp.e tmp.r -s
rm tmp.*
}
splicing.count_abc(){
usage="
FUNCTION: count number of contiguous reads (total, left crossing, right crossing)
USAGE: $FUNCNAME <target.bed> <read.bed> [options]
OPTIONS:
	-s : count reads on the same strand
	-S : count reads on the opposite starnd 
"
if [ $# -lt 2 ];then echo "$usage"; return; fi
	awk 'NF <= 6 || $10==1' $2 \
	| intersectBed -a ${1/-/stdin} -b stdin -wa -wb ${@:3} \
	| awk -v OFS="@" '{
		c=0;l=0;r=0;
		if( $2 > $8 ){ l=1;}
		if( $3 < $9 ){ r=1;}
		if( $2 <= $8 && $3 >= $9){ c=1;}
		print $1,$2,$3,$4,$5,$6"\t"c"\t"l"\t"r;
	}' | hm stat sum - | tr "@" "\t"  \
	| awk -v OFS="\t" '{
		print $1,$2,$3,$4,$7","$8","$9,$6;
	}'
}
splicing.count_abc.test(){
echo \
"01234567890123456789012345678901234567890123456789
      AAAAAAA
    RRRR
      RRRRRR
     RRRRRR
        RRRRRRR-RRRRR
        RRRR
             rr
                 rrrrrr
" | hm bed toy - | tail -n+2 > tmp.i
awk 'toupper($4)=="R"' tmp.i | cut -f1-6 > tmp.r
awk 'toupper($4)!="R"' tmp.i | cut -f1-6 > tmp.e
	splicing.count_abc tmp.e tmp.r -s
rm tmp.*
}

splicing.count_jc(){
usage="$FUNCNAME <bed12>" 
if [ $# -lt 1 ];then echo "$usage"; return; fi
	local tmp=${2:-"tmp"};
	awk -v OFS="\t" '{ $4="j"; $5=1;}1' $1 \
	| hm bed intron - \
	| hm bed flank - 1 1 \
	| hm bed sum - 
}


splicing.toy(){
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

splicing.toy.test(){
echo \
"01234567890123456789012345678901234567890123456789
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
 EEEE-------------EEEEEEEE------------EEEEEEEEEEEE
  RRR-------------RR
   rr-------------rrrrrr
                RRRR 
                                    RRRRR
" | splicing.toy -
}


splicing.table_exonpos(){
usage="
FUNC: calculate relative positions of exons from the minimum gene they belong to
$FUNCNAME <gene.bed12> <exon.bed> <output>
"
if [ $# -ne 3 ];then echo "$usage"; return; fi
overlap(){
	cat $1 | perl -e 'use strict; my %res=(); my $option="'$2'";
	while(<STDIN>){ chomp; my @a=split/\t/,$_;
		my $k=$a[0].":".$a[3].":".$a[5];
		if(defined $res{$k} ){
			if($option eq "max"){
				$res{ $k }{S}= $a[1] < $res{$k}{S} ? $a[1] : $res{$k}{S};  
				$res{ $k }{E}= $a[2] > $res{$k}{E} ? $a[2] : $res{$k}{E};  
			}else{
				$res{ $k }{S}= $a[1] > $res{$k}{S} ? $a[1] : $res{$k}{S};  
				$res{ $k }{E}= $a[2] < $res{$k}{E} ? $a[2] : $res{$k}{E};  
			}
		}else{
			$res{ $k }{S}=$a[1]; 
			$res{ $k }{E}=$a[2];
		}
	}
	foreach my $k (keys %res){
		my ($c,$g,$st)=split/:/,$k;
		print $c,"\t",$res{$k}{S},"\t",$res{$k}{E},"\t",$g,"\t0\t$st\n";
	}
	'
}
	hm util mkdir $3;
	cat $1 | cut -f1-6 > $3.gene;
	cat $2 > $3.exon;

	tail -n+2 $3.exon | overlap - max > $3.maxexon;
	intersectBed -a $3.gene -b $3.maxexon -F 1 -wa -wb -s  \
	| awk '$4==$(10)'\
	| cut -f 1-6 \
	| overlap - min > $3.mingene

	head -n 1 $3.exon | awk -v OFS="\t" '{ print $0,"relpos";}' > $3

	local nf=`awk '{print NF;}' $3.exon | head -n 1`;
	tail -n+2 $3.exon | intersectBed -a stdin -b $3.mingene -wa -wb -s \
	| perl -ne 'chomp; my @a=split/\t/,$_; my $M='$nf'; 
		my $s=$a[1]-$a[$M+1]; if($a[5] eq "-"){ $s=$a[$M+2]-$a[2];} 
		print join("\t",(@a[0..($M-1)],$s)),"\n";
	' >> $3;

}

splicing.table_exonpos.test(){
echo \
"
01234567890123456789012345678901234567890123456789
    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
       GGGGGGGGGGGGGGGG
         GGGGGGGGGGGGGGGGGGGGGG
         EEEEE
           EEEEEEEEEEEE
                     EEEEE
" | splicing.toy - > tmp.all
awk '$4 == "G"' tmp.all > tmp.gene
awk '$4 == "E"' tmp.all | cut -f1-6 > tmp.exon
	echo "==output==";
	splicing.table_exonpos tmp.gene tmp.exon
	rm tmp.*
}

splicing.psi_es(){
usage="
$FUNCNAME <transcript.bed12> <read.bed12> <output> [options]
     /             E              \
    ]----------[        ]----------[
               ---    --- ( I )
     \   I    /          \    I   /              
"
if [ $# -lt 3 ];then echo "$usage"; return; fi
	cat $1 | perl -ne 'chomp; my @a=split/\t/,$_;
		my @sizes=split/,/,$a[10];
		my @starts=split/,/,$a[11];
		for(my $i=1; $i < $a[9]-1; $i++ ){
			my $lexon_start = $a[1] + $starts[$i-1];
			my $lexon_end = $lexon_start + $sizes[$i-1];
			my $exon_start = $a[1] + $starts[$i];
			my $exon_end = $exon_start + $sizes[$i];
			my $rexon_start = $a[1] + $starts[$i+1];
			my $rexon_end = $rexon_start + $sizes[$i+1];
			my $content=join(",",($lexon_start,$lexon_end,$rexon_start,$rexon_end));
			my $content="null";
			print $a[0],"\t",$exon_start,"\t",$exon_end,"\t",$a[3],"\t",0,"\t",$a[5],"\n";
		}	
	' | sort -u > $3.exon
	intersectBed -a $3.exon -b ${2/-/stdin} ${@:4} -wa -wb \
	| perl -e 'use strict; 
	my %res=();
	my $ancor_sum=0;
	my $ancor_num=0;

	while(<STDIN>){ chomp; my @a=split/\t/,$_;
		my @l=split/,/,$a[16];
		my @s=split/,/,$a[17];
		my $id=join("\t",@a[0..5]);
		my ($les,$lee,$res,$ree) = split /,/,$a[4];
		my $max_l=$l[0];
		for(my $i=1; $i < $a[15]; $i++){
			$max_l=$l[$i] if $l[$i] > $max_l;
			### consider exact posisionts
			my $js= $a[7] + $s[$i-1]+$l[$i-1];
			my $je= $a[7] + $s[$i];
			#if ( $lee == $js  && $a[1] == $je || $a[2] == $js && $res == $je){
			if ( $a[1] == $je ){
				$res{$id}{I} ++;
			}elsif( $a[2] == $js ){
				$res{$id}{I} ++;
			#}elsif( $lee == $js && $res == $je ){
			}elsif( $a[1] > $js && $a[2] < $je ){
				$res{$id}{E} ++;
			} 
		}
		if( $a[15] > 1){
			$ancor_sum += $max_l;
			$ancor_num ++;
		}
		if( $a[15] == 1 && $a[1] <= $a[7] && $a[8] <= $a[2]){
			$res{$id}{I} ++;
		}

	}
	#print $ancor_sum/$ancor_num,"\n";
	## put here to normalize PSI;
	foreach my $k (keys %res){
		print $k;
		foreach my $i ( ("E","I")){
			my $v= defined $res{$k}{$i} ? $res{$k}{$i} : 0;
			print "\t$v";
		}
		print "\n";
	}
	' > $3
	rm $3.exon
}
splicing.psi_es.test(){
echo \
"
01234567890123456789012345678901234567890123456789
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
 EEEE-------------EEEEEEEE------------EEEEEEEEEEEE
 EEEE-------------EEEEEEEE-----------EEEEEEEEEEEE
  RRR-------------RR
                  RRR
                       RRRR
                         R-----------RRRRR
                         R------------RRRRR
    R---------------------------------RRRRR
    R--------------------------------RRRRR
       R---------------------------RRRRR
   rr-------------rrrrrr
                RRRR 
                         RRRRR
" | splicing.toy - > tmp.all 
cat tmp.all | awk '$4=="E"' > tmp.trans
cat tmp.all | awk 'toupper($4)=="R"' > tmp.read
splicing.psi_es tmp.trans tmp.read tmp.out -s
head tmp*
rm tmp.all tmp.trans tmp.read
}
splicing.relpos(){
usage="$FUNCNAME <exon.bed> <read.bed12> <u5p>,<d5p>,<u3p>,<d3p> <exp_type>
	--------[ exon          ]------------
	  | u5p | u3p |    | u3p| d3p}

 exp_type: netseq3p, rnaseq 
"
if [ $# -lt 1 ]; then echo "$usage"; return; fi
	if [ $4 == "netseq3p" ];then
	local tmpd=`hm util mktempd`;
	local p=( `echo $3 | tr "," " "` );
	cat $1 | hm bed enc | hm bed flank - ${p[0]} ${p[3]} -s > $tmpd/e
	awk '$10==1' $2 | hm bed 5p - \
	| intersectBed -a $tmpd/e -b stdin -wa -wb -S \
    | perl -e 'use strict;  my %res=(); my @b=('$3');
    	while(<STDIN>){chomp; my @a=split/\t/,$_;
			my $k=$a[3]; $k=~s/@/\t/g;
			my $l=$a[2]-$a[1];
			my $p=$a[7]-$a[1]; $p = $l - $p - 1 if $a[5] eq "-";
			if( $p <= $b[0] + $b[1]){
				$res{$k}{5}{$p-$b[0]}++;
			}elsif( $p >= $l - $b[2] -$b[3]-1){ 
				#print $p," ",$p-$l+$b[3]+1,"\n";
				$res{$k}{3}{$p-$l+$b[3]+1}++;
			}
       }
       foreach my $k (keys %res){
       		print $k;
			if( defined $res{$k}{5}){
                        	print "\t",join( ",", map {"$_:$res{$k}{5}{$_}"} keys %{$res{$k}{5}});
			}else{
				print "\tnull";
			}

			if( defined $res{$k}{3}){
                        	print "\t",join( ",", map {"$_:$res{$k}{3}{$_}"} keys %{$res{$k}{3}});
			}else{
				print "\tnull";
			}
			print "\n";
                }
        '
	rm -rf $tmpd;
	else
		echo "Error: $4 is not a known type!"; 
		echo "$usage";
	fi

}
splicing.relpos.test(){
echo \
"
012345678901234567890123456789
          EEEEEEEEEE
          eeeeeeeeee    
     rrr
      rrr
       rrr
        rrr
         rrr
          rrr
                    rrrrrr
                   rrrrrr
                  rrrrrr
                 rrrrrr
                rrrrrr
               rrrrrr
              rrrrrr
             rrrrrr
            rrrrrr
           rrrrrr
          rrrrrr
                   RRRRRR
      RRRRRR

" | splicing.toy - > tmp.all
cat tmp.all | awk 'toupper($4)=="E"' | cut -f 1-6 > tmp.e 
cat tmp.all | awk 'toupper( $4)=="R"' > tmp.r 
splicing.relpos tmp.e tmp.r 1,2,3,4 netseq3p 
rm tmp.*
}

splicing.relpos_to_table(){
usage="$FUNCNAME <relpos>"
if [ $# -lt 1 ];then echo "$usage"; return; fi
	cat $1 | perl -e 'use strict; 
	sub parse{
		my ($id, $i, $r, $c,$t)= @_;
		return if $i eq "null";
		foreach my $xy (split/,/,$i){
			my ($x,$y)=split/:/,$xy;
			$c->{$x} ++;
			$r->{$id}->{$t}->{$x} += $y;
		}
	}
	my %res=(); 
	my %col5=(); my %col3=();
	while(<STDIN>){ chomp; my @a=split/\t/,$_;
		parse(join("\t",@a[0..5]), $a[6], \%res, \%col5,5);
		parse(join("\t",@a[0..5]), $a[7], \%res, \%col3,3);
	}
	my @c5=sort {$a<=>$b} keys %col5;
	my @c3=sort {$a<=>$b} keys %col3;
	print join("\t",("chrom","start","end","name","score","strand")),"\t";
	print join("\t", map { "F.".$_ } @c5),"\t";
	print join("\t", map { "T.".$_ } @c3),"\n";
	foreach my $k ( keys %res ){
		print $k;
		foreach my $j (@c5){
			my $v= defined $res{$k}{5}{$j} ? $res{$k}{5}{$j} : 0;
			print "\t$v";
		}
		foreach my $j (@c3){
			my $v= defined $res{$k}{3}{$j} ? $res{$k}{3}{$j} : 0;
			print "\t$v";
		}
		print "\n";
	}
	' 
}

splicing.relpos_to_table.test(){
echo \
"chr1	100	200	n1	0	+	-4:1,-3:2,1:2,3:4	null
chr1	100	200	n1	0	+	-4:1,-3:2,1:2,3:4	-1:1,-2:2,10:3
chr1	200	400	n2	0	+	null	-4:1,-3:2,1:2,3:4" \
| splicing.relpos_to_table  -
}

splicing.table(){
usage="$FUNCNAME <ctr.ei>[,<ctr.ei..] <trt.ei>[,<trt.ei>] [options]
 [options]
"
if [ $# -ne 2 ];then echo "$usage"; return; fi
perl -e 'use strict;
	my $option="";
	my @fctr=map{"$_"} split /,/,"'$1'";
	my @ftrt=map{"$_"} split /,/,"'$2'";

	sub readf{
		my ($f,$r,$tag,$opt,$cols)=@_;
		open(my $fh,"<",$f) or die "$! : $f";
		while(<$fh>){chomp; my @a=split/\t/,$_;
			for(my $j=0; $j < $#a-5; $j++){
				my $tagj=$tag.".c".$j;
				$cols->{$tagj}=1;
				$r->{join("\t",@a[0..5])}{$tagj}=$a[$j+6];
			}
		}	
		close($fh);
	}
	my %res=();
	my $i=0;
	my %cols=();
	foreach my $f (@fctr){
		readf($f,\%res,"ctr".$i,$option,\%cols); $i++;
	}
	$i=0;
	foreach my $f (@ftrt){
		readf($f,\%res,"trt".$i,$option,\%cols); $i++;
	}
	print join("\t",("chrom","start","end","name","score","strand")),"\t";
	print join("\t",sort keys %cols),"\n";
	foreach my $k (keys %res){
		print $k;
		foreach my $c (sort keys %cols){
			my $v=0; $v=$res{$k}{$c} if defined $res{$k}{$c};
			print "\t",$v;
		}
		print "\n";
	}
'

}

splicing.table.test(){
	echo \
"chr1	100	200	s	0	+	1	2	
chr1	100	101	u	0	+	3	4
chr1	50	200	s	0	-	5	6
chr1	200	201	u	0	-	7	8
chr1	199	200	u	0	+	9	10" > tmp.a

splicing.table tmp.a tmp.a,tmp.a

	echo \
"chr1	100	200	s	0	+	1
chr1	100	101	u	0	+	3
chr1	50	200	s	0	-	5
chr1	200	201	u	0	-	7
chr1	199	200	u	0	+	9" > tmp.a
splicing.table tmp.a tmp.a,tmp.a

}

splicing.table2dexseq(){
## remove entries having one exon 
cat $1 | perl -e 'use strict; my $first=1; my $res=();
	while(<STDIN>){chomp; my @a=split/\t/,$_;	
		if($first){ 
			print join("\t",("id","gid",@a[6..$#a],)),"\n";
			$first=0;
		}elsif( $a[3]=~/(\w+)\.(E[\d\.]+)/){
			my $id=join ("@",@a[0..5]);
			my $gid=$1;
			my $vals=join ("\t",@a[6..$#a]);
			$res{$gid}{$id}=$vals;
			
		}else{	
			die "4th column needs gene.exon names e.g. \w+\.E[\d\.]+";
		}
	}
	foreach my $gid (keys %res){
		my @ids=keys %{$res{$gid}};	
		next if $#ids < 1;
		foreach my $id (@ids){
			print $id,"\t",$gid,"\t",$res{$gid}{$id},"\n";
		}
	}
'
}
splicing.table2dexseq.test(){
echo \
"chr	start	end	name	score	strand	ctr1	ctr2	trt1
chr1	100	200	G1.E1	0	+	1	2	3
chr1	100	200	G1.E2	0	+	1	2	3
chr1	100	200	G2.E1	0	+	1	2	3" | splicing.table2dexseq -
}

splicing.table_filter(){
usage="$FUNCNAME <table> <min> <non_zero>";
if [ $# -lt 2 ]; then echo "$usage";return ; fi
	cat $1 | perl -ne 'chomp; my @a=split /\t/,$_; 
		my $sum=0;foreach my $e (@a[1..$#a]){ $sum += $e;} 
		if($sum >= '$2'  || $a[0] eq "id"){ print $_,"\n";}
	' 
}

