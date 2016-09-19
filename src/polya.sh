
polya.pos(){
usage="
$FUNCNAME <input.bed> > <out.bed> 
"
if [ $# -lt 1 ];then echo "$usage"; return; fi
## find a position of polyA site from alignments
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

