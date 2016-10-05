##http://www.dgp.toronto.edu/people/JamesStewart/378notes/22intervals/

## implement interval tree


interval.tree(){
cat $1 | perl -e 'use strict; 
	sub add{
		my ($i,$a,$b,$min, $max, $x)= @_;
		if( !defined $x){
			return;
		}
		if ( !defined $x->{s} ){ ## initiate this tree
			$x->{s}=$b;
			$x->{l}->{s}=$a;
			$x->{l}->{p}=$x;
			my %n=();
			$x->{r}= \%n; 
			return;
		}
		if ( $a <= $min && $max <= $b ){
			$x->{"i"}=$i;
		}else{
			if( $a < $x->{s} ){
				#print "add $a from $x->{s}\n";
				add($i,$a,$b, $min,$x->{s},$x->{l});
			}else{
				add($i,$a,$b, $x->{s},$max,$x->{r});
			}
			if( $b < $x->{s} ){
				add($i,$a,$b, $min,$x->{s},$x->{l});
			}else{
				add($i,$a,$b, $x->{s},$max,$x->{r});
			}
		}

	}
	sub trav{
		my ($T,$d) =@_;		
		if( ! defined $T){ return;}
		print " "x$d,$T->{s};
		print ".l:";
		trav($T->{l},$d++);
		print ".r;";
		trav($T->{r},$d++); 
		print "\n";
	}
	my %T=();
	my $i=0; 
	while(<STDIN>){ chomp; my ($a,$b)=split/\s+/,$_;
		add($i,$a,$b,"-inf","inf",\%T);
		$i++;
	}
	trav(\%T,0);

'
}

interval.tree.test(){
echo \
"10	20
18	22
15	25" \
| interval.tree -

}
