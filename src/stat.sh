
stat.sum(){
	cat $1 | perl -e 'use strict; my %res=(); my $H='${2:-"0"}';
	while(<STDIN>){ chomp; my @a=split/\t/,$_;
		if($H > 0 ){ print $_,"\n"; $H--;next;}
		if(!defined $res{$a[0]}){
			$res{$a[0]} = ();
		}	
		for(my $i=0; $i< $#a; $i++){
			$res{ $a[0] }[$i] += $a[$i+1]; 
		}
	}
	foreach my $k (keys %res){
		print $k,"\t",join("\t",@{$res{$k}}),"\n";
	}'
}
stat.sum.test(){
echo \
"a	1	2
b	3	4
a	10	20
b	30	40" | stat.sum - 
}
