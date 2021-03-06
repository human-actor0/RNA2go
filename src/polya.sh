
polya.pos(){
usage="
$FUNCNAME <input.bed> <type>  
 <type>: 3rev
	3rev : 3 prime end of the reverse complment  (default)

	          *
	    <------TTTTTTT      (read)
	===========AAAAAAAAAAAA (mRNA)
"
	cat $1 | perl -e 'use strict; my %res=();
		my $TYPE="'${2:-"3rev"}'";
		while(<STDIN>){ chomp; my @a=split /\t/,$_;
			if( $TYPE eq "3rev"){
				my ($p,$s)=($a[1],"-");
				if($a[5] eq "-"){ $p=$a[2]-1; $s="+"; }
				$res{$a[0]}{$s}{$p} ++;
			}
		}
		foreach my $c (keys %res){
		foreach my $s (keys %{$res{$c}}){
		foreach my $p (keys %{$res{$c}{$s}}){
			my $v=$res{$c}{$s}{$p};
			print join("\t",( $c, $p, $p+1,".",$v,$s)),"\n";
		}}}
	' 
if [ $# -lt 1 ];then echo "$usage"; return; fi
## find a position of polyA site from alignments
}

polya.pos.test(){
echo \
"chr1	100	200	r1	0	+
chr1	100	200	r2	0	+
chr1	100	200	r3	0	-" | polya.pos - 3rev
}


polya.filter(){ 
usage="
FUNCT: Filter inter-priming artefacts 
USAGE: $FUNCNAME <bed> <fasta>
REFER: heppard, S., Lawson, N.D. & Zhu, L.J., 2013.Bioinformatics (Oxford, England), 29(20), pp.2564–2571.
"
	if [ $# -ne 2 ]; then echo "$usage"; return; fi

	bed.flank $1 39 30 -s  | seq.read $2 - -s \
	| perl -ne 'chomp; my @a=split/\t/,$_; print join(";",@a[0..($#a-1)]),"\t",$a[$#a],"\n"' \
	| eval "$FILTER predict - $FILTER_M" \
	| perl -ne 'chomp; my ($bed,$seq,$pos,$score)=split/\t/,$_;
		my @a=split /;/,$bed;
		my $len=$a[2]-$a[1];
		my @p=split /,/,$pos;
		my @s=split /,/,$score;
		for(my $i=0; $i<=$#p;$i++){
			my $this_s=sprintf("%.3f",$s[$i]); 
			my $this_p=$p[$i];
			if ($#a > 4 && $a[5] eq "-"){
				$this_p=$len-$this_p-1;
			}
			print join("\t",($a[0],$a[1]+$this_p,$a[1]+$this_p+1,$a[3],$a[4],$a[5],$this_s)),"\n";
		}
	'
}

polya.filter.test(){
echo "hi";
}

