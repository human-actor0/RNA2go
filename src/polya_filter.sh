## data obtained from  https://raw.githubusercontent.com/Bioconductor-mirror/cleanUpdTSeq/release-3.3/inst/extdata/test.bed
#PeakName	prob False/oligodT internally primed	prob True	pred.class	UpstreamSeq	DownstreamSeq
data="6hpas-78439	0.999146181625342	0.000853818374657752	0	GGTCATTGTCCTGCAAAATGGACTACTTAACCGAACTGGA	GAAGTATAAGAAGTAAGTACATTAAAGCTAC
6hpas-78440	0.999999999999997	3.10583755675325e-15	0	TGGATTTAAATAACAAACAAGTTAAATAAAACGATTTGTA	AAAAAATAAAACAACTGAAGAAGAAAATGAA
6hpas-78441	4.44561326444396e-06	0.999995554386736	1	ATCTGCTTCAAAATGGATGCTCTGTTGAATCCTGAGCTCA	GGTAATCTTTCAAGTGCTGCTATTGAGCCAA
6hpas-78442	0.000478974691498084	0.999521025308502	1	AAATGCTTGCACATAATAAATGTAGGCTTAAAAGATTTCA	AAACGTTTGTGAGAGACGGATTTTACTTTGC
6hpas-22249	7.35421436011826e-07	0.999999264578564	1	TCTTCATCATGGTCATCTCGCACCAGAGAGTGTGCCAGGG	CAGGAAGTTTTACCTGTCTGTCATTATCGTC
6hpas-22250	0.93403928359027	0.0659607164097298	0	ACCCTGGTGAGGGTATAGAGCTGGTCCAGTGTGCCACGGC	AAAGAGGAAAACAGCATTGTTCCTCCTGGAT
6hpas-22251	0.999999999835442	1.645576640904e-10	0	TGATTTGTTTGTAACTGATTTTATCTTTTAATAAAAAAGA	AAAAAGAAAGTCAAGCCAAGAGGCAAATACA
6hpas-22252	0.974004754525814	0.0259952454741859	0	GGAGCGCGACCGCATCAACAAAATCTTGCAGGATTATCAG	AAGAAAAAGATGGTGAGTTATTATCATTCAA
6hpas-33204	0.99985217513107	0.000147824868930001	0	AGGGAAATAAATACAAAAGAATAAAAATATGATTCATTGT	AAGAAAAACACTTTAGCTACAAAAGTCCTTC
6hpas-33205	0.999999986000354	1.39996464368816e-08	0	ATTTAGTTGGGTATTATTTCAAATAAAGAGAGAGAGAGAC	ACAAAAACTACATCAAATTTGAGGACAAAAA
6hpas-33209	0.790767944303306	0.209232055696694	0	TCAAAGTTAATGTACATTAAAAATGAGTCAAAATGTTTAG	AATAAAAGAAGATTTGAATGATATATTCTTT
6hpas-33210	0.999872941841863	0.000127058158136887	0	TGAATGTATTTTCAAAGTTAATGTACATTAAAAATGAGTC	AAAATGTTTAGAATAAAAGAAGATTTGAATG
6hpas-33211	0.35084828810188	0.64915171189812	1	TTGGATAGTAAATTAATTATTTATAAAGTTTCTAGATTAC	ATAAAGAAAATAAATCTGTTATATCTGTATG
6hpas-33212	2.62255550185117e-10	0.999999999737744	1	TGATCTCCATATGATATCACCGTCCCTATTTAACTTAAAG	GTTTATCTTGTTTATAAGGGTGTGATAGAAA
6hpas-33213	0.999985145788329	1.48542116712713e-05	0	CCTCGATGATGCCGCCCGCAAAGCTGTCGCCGCCATTGCC	AAGAAATAAATGCAAATATTCATAATGCACA
6hpas-33214	0.999999914469604	8.55303959354549e-08	0	ATCGCCATTTTGCCCGTTCGTCATCGCATAAACCTGAGAC	AACCAAAAAAGGGCAAAGAGGCGGAGCTACA
6hpas-43855	1	4.8017325740641e-18	0	AAGGCCCAAACAGTAAAAAAAAATAAGACTGCTCTGCTTT	AAAAAAAAAAAAAAAAAAACCTTCAGTGGGA
6hpas-43856	0.297923854055364	0.702076145944636	1	AAGGTGTTTACATGTCTGTACTGCACTTCAATAATGTGAC	TAAAATAGGAATGCTCCAAATGGCTTCATTT
6hpas-43857	2.96212420355829e-05	0.999970378757964	1	TTACACAACGCTAATGGTTTTATTAAGGTGTTTACATGTC	TGTACTGCACTTCAATAATGTGACTAAAATA
6hpas-43858	4.22189144500321e-05	0.99995778108555	1	TTTATTTTAATAAATAAGCATTTTTAAAAGACTTCATATT	AATCAAACATTGTCTTGTCTATCATTGCCTA
6hpas-43859	1	5.0321956189001e-16	0	GAGAAGGAGAATGAGGAGAGTTTGAATCAAAATAATAATT	GAAAATAAAAAAAATAAAAAAAACTGGATGG
6hpas-43860	1	2.14246061067426e-17	0	GAGGAGAAGCAAGAGAAGGAGAATGAGGAGAGTTTGAATC	AAAATAATAATTGAAAATAAAAAAAATAAAA
6hpas-71514	0.999999999916475	8.35250718327141e-11	0	ATCTACAACCCCAAATCAGAAAAAGATTGGCACAGTATGG	AAAACACAAATAAAAAAGAAAGTGATTTACG
6hpas-71515	0.999995450433236	4.54956676366438e-06	0	TTTGTTACTTGAGACGCATCAAGATTTTATCTACAACCCC	AAATCAGAAAAAGATTGGCACAGTATGGAAA
6hpas-71516	5.21529000362569e-12	0.999999999994785	1	ATTTGGTCCGGATCAAGGGTAATAAATGACACATTGTTGC	ATTTTCTGCCGTCTTTGGGTCGTTTTCACAC
6hpas-71517	4.79012701408921e-06	0.999995209872986	1	GTTTTGAAATTGTGAGTATAAAGTAAATCTTTCAGTCATC	AGTGTTGAGTTTCATATACAGGAATCATGTA
6hpas-71518	0.993070695294844	0.00692930470515563	0	GCTTCACGGTTGCCCTCAGTGTGGGAAGAGCTTCACTTGG	AAAAAAACCCTTATTGAGCATATGAAGGTTC" 




program='
	use strict;
	my $PI = 3.141593;
	my $SQRT2PI= sqrt(2*$PI);
	sub mean_sd{
		my ($x, $n) = @_; ## $n: total to handle entiries with non-existing features
		my $esp = 0.001;
		my $sum = 0; my $sumsq = 0; 
		foreach my $v (@$x){
			$sum += $v; $sumsq += $v*$v;
		}
		my $mean = $sum/$n;
		my $std = ($n > 1)? sqrt(($sumsq - $sum * $sum/$n)/($n-1)): $esp;
		return ($mean,$std);
	}
	sub dnorm{
		my ($x,$mean,$sd) = @_;
		return (1/($SQRT2PI*$sd)*exp(-($x-$mean)**2/(2*$sd**2)));
	}
	#	print "1.6549~",dnorm(2.02,2,0.24),"\n";return;
	#my ($m,$std) = mean_sd([1,2,3,4,5]); print $m," ",$std,"\n";

	sub makeMotifFeature{
		my ($seq, $motif_len, $tag, $F) = @_;
		for(my $i=0;$i<length($seq) - $motif_len + 1;$i++){
			my $motif = substr($seq,$i,$motif_len);
			$F->{$tag.":".$motif} ++;
		}
	}
	sub makeMotifDistalFeature{
		my ($seq, $motif, $tag, $F) = @_;
		my $s = 0; my $n = 0;
		for(my $i=0;$i<length($seq) - length($motif) + 1;$i++){
			if(substr($seq,$i,length($motif)) eq $motif){
				$s += $i+1; # 1-base
				$n ++;
			}
		}
		$F->{$tag.":".$motif} = ($n >0)? $s/$n: length($seq);
	}
	sub makeFeature{
		my ($upseq,$dnseq) = @_;
		my %F = ();
		makeMotifFeature($upseq,6,"B",\%F);
		makeMotifFeature($dnseq,2,"N",\%F);
		makeMotifFeature($dnseq,1,"N",\%F);
		makeMotifDistalFeature($dnseq,"A","D",\%F);
		return \%F;
	}
	sub updateFeature{
		my ($F,$f,$op) = @_;
		print $f," ",$F->{$f}," ",$op,"\n";	
		my $type = [split /:/,$f]->[0];
		
		if($type eq "D" || $type eq "C"){ 
			if($op eq "inc"){ $F->{$f} ++;
			}elsif(defined $F->{$f}){
				$F->{$f} --;
				if ($F->{$f} == 0){ delete $F->{$f}; }
			}
		}elsif($type eq "B"){
			my ($ave,$n) = (split /:/,$F->{$f});
			if($op eq "inc"){
				$ave = $n/($n+1)*$ave + (30)/($n+1);
				$n++;
			}else{
				$ave =($n+1)/$n*$ave - (1)/$n;	 #  ave1 = n / n+1 ave0  + v/n+1 => ave0=  n+1/n (ave1 - v/n+1 )
				$n--;
			}
		}
	}
		
	sub printFeature{
		my ($f) = @_;
		foreach my $k (keys %$f){ 
			print $k,"\t",$f->{$k},"\n";
		}
	}
	sub collectFeatures{
		my ($y, $f, $F) = @_;
		foreach my $k (keys %$f){	
			if(!defined $F->{$y}->{$k}){
				$F->{$y}->{$k} =  [];
			}
			push @{$F->{$y}->{$k}}, $f->{$k};
		}
	}
	sub genMotifs{
		my ($motif,$m_len) = @_;
		return $motif if $m_len == 0;
		my @res= ();
		for my $nu ("A","C","G","T"){
			push @res,genMotifs($motif.$nu,$m_len-1);
		}
		return @res;
	} #print join("\n", genMotifs("",4)); exit(1);
	sub buildModel{
		my ($F, $NN) = @_;
		my %M = ();
		my $laplace = 1; ## e1071
		my $ny = scalar keys %$F;

		my $TOT=0;
		foreach my $y (keys %$F){
			my $N = $NN->{$y}; $TOT += $N;
			foreach my $k (genMotifs("",6)){
				$k = "B:".$k;
				my $v = $F->{$y}->{$k};
				my $n = (defined $v)? scalar @{$v}: 0;
				my $p=($n+$laplace)/($N + $laplace * $ny);	
				$M{$y}{$k} = $p;
				$M{$y}{"b:"} += log(1-$p);
			}
			
			foreach my $k (keys %{$F->{$y}}){
				my $type = [split /:/,$k]->[0];
				my $v = $F->{$y}->{$k};
				if($type eq "N" || $type eq "D"){ ## numeric or distal
					my ($mean,$sd) = mean_sd($v, $N); 
					$M{$y}{$k}= "$mean,$sd";
				}
			}
		}	
		foreach my $y (keys %$NN){ ## apriori
			$M{$y}{"T:"} = $NN->{$y}/$TOT;
		}
		return \%M;
	}
	sub writeModel{
	## class\ttype\tfea\tvalues
	## type::B : p(Y|X) = # entries containing X / total
	## type::N|D : p(Y|X) = dnorm(mu, sd)
		my ($fh,$M) = @_;
		foreach my $y (keys %{$M}){
		foreach my $k (keys %{$M->{$y}}){
			print {$fh} $y,"\t",$k,"\t",$M->{$y}->{$k},"\n";
		}}
	}
	sub readModel{
		my ($fh) = @_;
		my %M = ();
		while(<$fh>){ chomp;
			my ($y,$fea,$v) = split /\t/,$_;
			$M{$y}{$fea}=$v;
			#print $y," ",$fea,"\n";
		}
		return \%M;
	}
	sub s2h{
		my ($s,$h) = @_;
		my %h = ();
		foreach my $kv (@{[split /,/,$s]}){
			my ($k,$v) = split /:/,$kv;
			$h->{$k}=$v;
		}
	}
	sub predict{
		my ($M,$F) = @_;
		my $L=0;
		my $eps = 0.001;
		$L += log( $M->{"T:"}); ## apriori
		$L += $M->{"b:"}; ## sum of non-existence prob.
		foreach my $k (keys %$F){
			my $type = [split /:/,$k]->[0];
			my $v = $F->{$k};
			my $m = $M->{$k};
			next unless defined $m;
			if($type eq "B"){
				$L += log($m);
				$L -= log(1-$m); ## remove non-existance probs
			}elsif($type eq "N" || $type eq "D"){
				my ($mean,$sd) = split /,/,$m;
				my $p = dnorm($v,$mean,$sd);
				$L += log($p);
			}
		}
		return($L);
	}

	if("CMD" eq "predict"){
		open(my $fh, "<","MODEL");
		my $M = readModel($fh);
		close($fh);
		my $up_len = 40;
		my $dn_len = 30;

		my @classes = sort {$a<=>$b} keys %$M;
		while(<>){ chomp;
			my ($id,$seq) = split /\t/,$_;
			$seq=uc $seq;
			my @res_x=();
			my @res_y=();

			for(my $i=0; $i < length($seq) - $up_len - $dn_len + 1; $i++){
				my $up=substr($seq,$i,$up_len);
				my $dn=substr($seq,$i+$up_len,$dn_len);
				my $F = makeFeature($up,$dn);
				my %llh = ();
				my $denom = 0;
				foreach my $y (@classes){
					$llh{$y}=predict($M->{$y},$F);
					$denom += exp($llh{$y});
				}
				foreach my $y (@classes){
					#push @res, ($i+$up_len).":".$y.":".exp($llh{$y})/$denom;
					if($y == 1){
						push @res_x,($i+$up_len-1);
						push @res_y, exp($llh{$y})/$denom;
					}
				}
			}
			#print $id,"\t",$seq,"\t",join( ",",@res),"\n";
			print $id,"\t",$seq,"\t",join( ",",@res_x),"\t",join(",",@res_y),"\n";
		}
	}elsif("CMD" eq "train"){
		my %F = ();
		my %N = ();
		while(<>){ chomp;
			my ($y, $ups,$dns) = split /\t/,$_;
			my $f = makeFeature($ups,$dns);
			#print $_,"\n"; printFeature($f);
			collectFeatures($y,$f,\%F);
			$N{$y}++;
		}	
		my $M = buildModel(\%F,\%N);
		writeModel(*STDOUT,$M);
	}
'
MODEL=${BASH_SOURCE%/*}/polya_filter_nb.model
cmd=$program; 
cmd=${cmd//CMD/predict}; 
cmd=${cmd/MODEL/$MODEL};

echo "$data" | tail -n+2 | awk -v OFS="\t" '{ print $3,$5$6;}' \
| perl -e "$cmd"

