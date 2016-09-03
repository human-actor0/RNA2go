
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

 exp_type: 3pnetseq, rnaseq 
"
if [ $# -lt 1 ]; then echo "$usage"; return; fi
	if [ $4 == "3pnetseq" ];then
	local tmpd=`hm util mktempd`;
	local p=( `echo $3 | tr "," " "` );
	cat $1 | hm bed enc | hm bed flank - ${p[0]} ${p[3]} -s > $tmpd/e
	awk '$10==1' $2 | hm bed 5p - \
	| intersectBed -a $tmpd/e -b stdin -wa -wb -S \
        | perl -e 'use strict;  my %res=(); my @b=('$3');
                while(<STDIN>){chomp; my @a=split/\t/,$_;
			my $k=join("\t",@a[0..5]);
			my $l=$a[2]-$a[1];
			my $p=$a[7]-$a[1]; $p = $l - $p if $a[5] eq "-";
			if( $p <= $b[0] + $b[1]){
				$res{$k}{5}{$p-$b[0]}++;
			}elsif( $p >= $l - $b[2] -$b[3]){ 
				print $p," ",$p-$l+$b[3]-1,"\n";
				$res{$k}{3}{$p-$l+$b[3]-1}++;
			}
                }
                foreach my $k (keys %res){
                        print $k,"\t";
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
     EEEEEEEEEEEEEEE
     eeeeeeeeeeeeeee    
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

" | splicing.toy - > tmp.all
cat tmp.all | awk 'toupper($4)=="E"' | cut -f 1-6 > tmp.e 
cat tmp.all | awk 'toupper( $4)=="R"' > tmp.r 
splicing.relpos tmp.e tmp.r 1,2,3,4 3pnetseq 
}
