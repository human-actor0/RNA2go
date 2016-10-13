
hisat2.reptab(){
perl -e 'use strict;
	my %res=();
	my $id="";
	while(<STDIN>){chomp;
		if ($_=~/==> (\w+).(sam|bam)/){
			$id=$1;
		}elsif($_=~/(\d+) reads;/){
			$res{$id}{tot}=$1;
		}elsif($_=~/(\d+ \(.+\)) aligned exactly 1 time/){
			$res{$id}{uniq}=$1;
		}elsif($_=~/(\d+ \(.+\)) aligned >1 times/){
			$res{$id}{multi}=$1;
		}elsif($_=~/(.+) overall/){
			$res{$id}{mapped}=$1;
		}
	}
	my @C=("tot","uniq","multi","mapped");
	print join("\t",("id",@C)),"\n";
	foreach my $id (keys %res){
		print $id,"\t";
		print join("\t",map { $res{$id}{$_} } @C),"\n";
	}
'
}

hisat2.reptab.test(){
echo \
"==> hct116_wt_GTGGCC.sam.bsub.err.536735 <==
99596344 reads; of these:
  99596344 (100.00%) were unpaired; of these:
    11406734 (11.45%) aligned 0 times
    74452665 (74.75%) aligned exactly 1 time
    13736945 (13.79%) aligned >1 times
88.55% overall alignment rate
" | hisat2.reptab -


}
