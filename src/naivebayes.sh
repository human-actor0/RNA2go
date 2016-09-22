
naivebase.predict(){
usage="$FUNCNAME <fea.txt> <model.txt>
"
if [ $# -ne 2 ];then echo "$usage"; return; fi
perl -e 'use strict; 
	my $f_fea="'$1'";
	my $f_mod="'$2'";
	my $PI= 3.141593; my $SQRT2PI= sqrt(2*$PI);
	my $EPS=0; 
	my $THRE=0.001; # min prob obtained from e1071.naiveBayes
	sub dnorm{
		my ($x,$mean,$sd) = @_;
		my $p = (1/($SQRT2PI*$sd)*exp(-($x-$mean)**2/(2*$sd**2)));
		return $p;
	}

	my %M=(); ## model
	open( my $fh, "<", $f_mod) or die "$f_fea";
	while(<$fh>){ chomp; my ($f,@tmp)=split /\t/,$_;
		foreach my $tmp1 (@tmp){
			my ($y,$tmp2)=split/:/,$tmp1;
			my @p=split/,/,$tmp2;
			$M{$f}{$y}=\@p;
			#print $f," ",$y," ",join(",",@{$M{$f}{$y}}),"\n";
		}
	}
	close($fh);

	my $n=0;
	my %ny=();
	foreach my $y (keys %{$M{t}}){
		$ny{$y}=$M{t}{$y};
		$n += $ny{$y};
	}

	open( my $fh, "<", $f_fea) or die "$f_fea";
	while(<$fh>){ chomp; my ($id,@a)=split /\t/,$_;
		my %logprob=();
		## prior p(y);
		foreach my $y (keys %ny){
			$logprob{$y} = log( $ny{$y} ) - log($n); 
		}
		foreach my $fv (@a){
			my ($f,$v) = split /:/,$fv;
			foreach my $y (keys %{$M{$f}}){
				my $prob=0;
				my @theta=@{$M{$f}{$y}};
				if($f=~/^n\./){
					#print $f," ",$y," ",$v," ",$theta[0]," ",$theta[1],"\n";
					$prob=dnorm( $v, $theta[0],$theta[1]);
				}elsif($f=~/^m\./){
					## theta^n
					$prob=$theta[0] ** $v;
				}elsif($f=~/^b\./){
					$prob=$theta[0] * $v + (1-$theta[0]) * (1-$v);
				}
				## set minimum prob
				if ($prob <= $EPS){ $prob=$THRE; }
				$logprob{$y} += log($prob);
			}
		}
		print join("\t", map{ "$_:$logprob{$_}" } sort keys %logprob),"\n";
	}
	close($fh);
'
}


naivebase.param(){
usage="$FUNCNAME <txt> [<alpha>]
 <alpha> : prior counts for smoothing ( default = 1 )
"
if [ $# -lt 1 ];then echo "$usage"; return; fi
	cat $1 | perl -e 'use strict; 
	my $ALPHA='${2:-1}'; ## smoothing
	my %Nyi=();   # sum of feature i in class y 
	my %N2yi=();  ## N^2_yi
	my %Ny=(); ## N_y
	my $N=0;   ## N
	my $n=0;  # number of instances
	my %ny=(); # number of instances per class


	while(<STDIN>){chomp; my ($y,@fea)=split/\t/,$_;
		$n++;
		$ny{$y}++;
		foreach my $tmp (@fea){
			if( $tmp=~/([nmb]\.\w+):(\d+)/){
				my ($f,$v)=($1,$2);
				$Nyi{$f}{$y} += $v;
				$N2yi{$f}{$y} += $v*$v;
				$Ny{$y} += $v;
				$N += $v;
			}
		}
	}

	sub mean_sd{
		my ($sum,$sumsq,$n) = @_; 
		my $esp = 0.001; ## minimum variance
		my $mean = $sum/$n;
		my $std = ($n > 1)? sqrt(($sumsq - $sum * $sum/$n)/($n-1)): $esp;
		if($std < $esp){ $std=$esp;}
		return ($mean,$std);
	}

	print "t\t",join("\t",map { "$_:$ny{$_}" } keys %ny),"\n";
	foreach my $f (keys %Nyi){
		print $f;
		foreach my $y (sort keys %ny){
			## calculate p(x_i|y)
			my $sum= defined $Nyi{$f}{$y} ? $Nyi{$f}{$y} : 0;
			my $sum2= defined $N2yi{$f}{$y} ? $N2yi{$f}{$y} : 0;
			my $p="NULL";
			my $n_y=$ny{$y};
			my $N_y=$Ny{$y};
			if($f=~/^n\./){
				my ($mu,$sig)=mean_sd($sum,$sum2,$n_y);
				$p="$mu,$sig";
			}elsif($f=~/^m\./){
				$p=($sum+$ALPHA)/( $N_y + 2*$ALPHA *$n);
			}elsif($f=~/^b\./){
				$p=($sum+$ALPHA)/( $N_y + 2*$ALPHA );
			}
			print "\t$y:$p";
		}
		print "\n";
	}
	'
}

naivebayes.test(){
echo \
"0	b.A:1	m.A:1	n.A:10
0	b.A:1	m.A:1	n.A:10
1	b.B:1	m.B:1	n.B:20" > tmp.fea 
naivebase.param tmp.fea  > tmp.param
naivebase.predict tmp.fea tmp.param

head tmp.*
rm tmp.*
}
