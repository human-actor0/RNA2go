
data="chr10 2965327 2965327 6hpas-22249   1       -     TCTTCATCATGGTCATCTCGCACCAGAGAGTGTGCCAGGG      CAGGAAGTTTTACCTGTCTGTCATTATCGT
chr10 2966558 2966558 6hpas-22250   1       -     ACCCTGGTGAGGGTATAGAGCTGGTCCAGTGTGCCACGGC      AAAGAGGAAAACAGCATTGTTCCTCCTGGA
chr10 2974251 2974251 6hpas-22251   2       -     TGATTTGTTTGTAACTGATTTTATCTTTTAATAAAAAAGA      AAAAAGAAAGTCAAGCCAAGAGGCAAATAC
chr10 2978441 2978441 6hpas-22252   1       -     GGAGCGCGACCGCATCAACAAAATCTTGCAGGATTATCAG      AAGAAAAAGATGGTGAGTTATTATCATTCA
chr11 16772291        16772291        6hpas-33204   1       -     AGGGAAATAAATACAAAAGAATAAAAATATGATTCATTGT      AAGAAAAACACTTTAGCTACAAAAGTCCTT
chr11 16777848        16777848        6hpas-33205   1       -     ATTTAGTTGGGTATTATTTCAAATAAAGAGAGAGAGAGAC      ACAAAAACTACATCAAATTTGAGGACAAAA
chr11 17122845        17122845        6hpas-33209   1       -     TCAAAGTTAATGTACATTAAAAATGAGTCAAAATGTTTAG      AATAAAAGAAGATTTGAATGATATATTCTT
chr11 17122856        17122856        6hpas-33210   2       -     TGAATGTATTTTCAAAGTTAATGTACATTAAAAATGAGTC      AAAATGTTTAGAATAAAAGAAGATTTGAAT
chr11 17123062        17123062        6hpas-33211   1       -     TTGGATAGTAAATTAATTATTTATAAAGTTTCTAGATTAC      ATAAAGAAAATAAATCTGTTATATCTGTAT
chr11 17123194        17123194        6hpas-33212   1       -     TGATCTCCATATGATATCACCGTCCCTATTTAACTTAAAG      GTTTATCTTGTTTATAAGGGTGTGATAGAA
chr11 17123754        17123754        6hpas-33213   3       -     CCTCGATGATGCCGCCCGCAAAGCTGTCGCCGCCATTGCC      AAGAAATAAATGCAAATATTCATAATGCAC
chr11 17204740        17204740        6hpas-33214   1       -     ATCGCCATTTTGCCCGTTCGTCATCGCATAAACCTGAGAC      AACCAAAAAAGGGCAAAGAGGCGGAGCTAC
chr12 25855374        25855374        6hpas-43855   1       -     AAGGCCCAAACAGTAAAAAAAAATAAGACTGCTCTGCTTT      AAAAAAAAAAAAAAAAAAACCTTCAGTGGG
chr12 26155099        26155099        6hpas-43856   2       -     AAGGTGTTTACATGTCTGTACTGCACTTCAATAATGTGAC      TAAAATAGGAATGCTCCAAATGGCTTCATT
chr12 26155123        26155123        6hpas-43857   2       -     TTACACAACGCTAATGGTTTTATTAAGGTGTTTACATGTC      TGTACTGCACTTCAATAATGTGACTAAAAT
chr12 26170452        26170452        6hpas-43858   1       -     TTTATTTTAATAAATAAGCATTTTTAAAAGACTTCATATT      AATCAAACATTGTCTTGTCTATCATTGCCT
chr12 26295950        26295950        6hpas-43859   3       -     GAGAAGGAGAATGAGGAGAGTTTGAATCAAAATAATAATT      GAAAATAAAAAAAATAAAAAAAACTGGATG
chr12 26295962        26295962        6hpas-43860   5       -     GAGGAGAAGCAAGAGAAGGAGAATGAGGAGAGTTTGAATC      AAAATAATAATTGAAAATAAAAAAAATAAA
chr15 733981  733981  6hpas-71514   10      -     ATCTACAACCCCAAATCAGAAAAAGATTGGCACAGTATGG      AAAACACAAATAAAAAAGAAAGTGATTTAC
chr15 734009  734009  6hpas-71515   2       -     TTTGTTACTTGAGACGCATCAAGATTTTATCTACAACCCC      AAATCAGAAAAAGATTGGCACAGTATGGAA
chr15 735146  735146  6hpas-71516   13      -     ATTTGGTCCGGATCAAGGGTAATAAATGACACATTGTTGC      ATTTTCTGCCGTCTTTGGGTCGTTTTCACA
chr15 735378  735378  6hpas-71517   5       -     GTTTTGAAATTGTGAGTATAAAGTAAATCTTTCAGTCATC      AGTGTTGAGTTTCATATACAGGAATCATGT
chr15 735597  735597  6hpas-71518   2       -     GCTTCACGGTTGCCCTCAGTGTGGGAAGAGCTTCACTTGG      AAAAAAACCCTTATTGAGCATATGAAGGTT
chr16 18304846        18304846        6hpas-78439   10      +     GGTCATTGTCCTGCAAAATGGACTACTTAACCGAACTGGA      GAAGTATAAGAAGTAAGTACATTAAAGCTA
chr16 18312684        18312684        6hpas-78440   1       +     TGGATTTAAATAACAAACAAGTTAAATAAAACGATTTGTA      AAAAAATAAAACAACTGAAGAAGAAAATGA
chr16 18316016        18316016        6hpas-78441   1       +     ATCTGCTTCAAAATGGATGCTCTGTTGAATCCTGAGCTCA      GGTAATCTTTCAAGTGCTGCTATTGAGCCA
chr16 18316389        18316389        6hpas-78442   1       +     AAATGCTTGCACATAATAAATGTAGGCTTAAAAGATTTCA      AAACGTTTGTGAGAGACGGATTTTACTTTG"

polya_filter.build_feature(){
cat $1 | perl -e 'use strict;
	while(<STDIN>){chomp; $_=~s/\s+/\t/g;
		my ($chrom,$start,$end,$name,$score,$strand,$upseq,$dnseq) = split /\t/,$_;
		## n.N.Downstream
		my %D=();
		my %d=();
		## m: multinomial, n:normal, b:binomial
		for(my $i=0; $i<length($dnseq); $i++){
			my $nu=substr($dnseq,$i,1);
			$D{"m:".$nu} ++;
			$d{"n:".$nu}{sum} += $i;
			$d{"n:".$nu}{num} ++;
		}
		foreach my $k (keys %d){
			$D{$k} = $d{$k}{sum}/$d{$k}{num};
		}
		for(my $i=0; $i<length($upseq) - 5; $i++){
			my $n6=substr($upseq,$i,6);
			$D{"b:".$n6} = 1;
		}
		print $_,"\t";
		print join("\t",( map{ "$_:$D{$_}"} keys %D)),"\n";
	}
'
}

polya_filter.build_feature.test(){
	echo "$data" | polya_filter.build_feature -
}
cleanUpdTSeq.build_feature(){
'
             library(cleanUpdTSeq)
             testFile = system.file("extdata", "test.bed", package="cleanUpdTSeq")
             testSet = read.table(testFile, sep = "\t", header = TRUE)
             
             #convert the test set to GRanges with upstream and downstream sequence information
             peaks = BED2GRangesSeq(testSet,upstream.seq.ind = 7, downstream.seq.ind = 8, withSeq=TRUE)
             #build the feature vector for the test set with sequence information 
             testSet.NaiveBayes = buildFeatureVector(peaks,BSgenomeName = Drerio, upstream = 40,
              downstream = 30, wordSize = 6, alphabet=c("ACGT"),
              sampleType = "unknown",replaceNAdistance = 30, 
             method = "NaiveBayes", ZeroBasedIndex = 1, fetchSeq = FALSE)

write.table(testSet.NaiveBayes$data,file="testSet.feature",sep="\t",quote=F)
             data(data.NaiveBayes)
             predictTestSet(data.NaiveBayes$Negative, data.NaiveBayes$Positive, testSet.NaiveBayes,
             outputFile = "test-predNaiveBayes.tsv", assignmentCutoff = 0.5)
'

}
