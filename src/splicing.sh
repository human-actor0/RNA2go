
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

