##http://www.dgp.toronto.edu/people/JamesStewart/378notes/22intervals/

## implement interval tree


interval.tree(){
cat $1 | perl -e 'use strict; 
	sub add{
		my ($i,$a,$b,$x)= @_;
		if ( !defined $x->{a} ){ ## insert this node 
			$x->{a}=$a;
			$x->{b}=$b;
			$x->{max}=$b;
			$x->{l}={};
			$x->{r}={};
			return $b;
		}else{
			my $max=0;
			if( $a < $x->{a} ){
				$max=add($i,$a,$b,$x->{l});
			}else{
				$max=add($i,$a,$b,$x->{r});
			}
			if($max > $x->{max}){
				$x->{max}=$max;
			}
			return $x->{max};
		}

	}
	sub trav{
		my ($T,$d) =@_;		
		if( ! defined $T->{a}){ return;}
		print "-"x$d,join(",",($T->{a},$T->{b},$T->{max})),"\n";
		trav($T->{l},$d+1);
		trav($T->{r},$d+1);
	}
	sub intersect{
		my ($a,$b,$x,$y)=@_;
		if( $b <= $x || $y <= $a){ return 0;}
		return 1;
	}
	sub search{
		my ($T,$a,$b,$d, $res) =@_;		
		if( ! defined $T->{a}){ return "";}
		if( intersect($T->{a},$T->{b},$a,$b) ){ 
			push @$res,$T->{a}."\t".$T->{b};
		}

		if(!defined $T->{l}->{max} || $T->{l}->{max} <= $a){
			search($T->{r},$a,$b,$d+1,$res);
		}else{
			search($T->{l},$a,$b,$d+1,$res);
		}
	}
	my %T=();
	my $i=0; 
	while(<STDIN>){ chomp; my ($a,$b)=split/\s+/,$_;
		add($i,$a,$b,\%T);
		$i++;
	}
	trav(\%T,0);
	my @res=();
	search(\%T,24,25,0,\@res);
	print "found--\n";
	print join ("\n",@res),"\n";

'
}

interval.tree.test(){
echo \
"10	20
15	25
18	22
14	30" \
| interval.tree -

}
