
extract_seq(){
cmd='
	## go to the data directory of the R package
	load("data.NaiveBayes.rda")
	d=data.NaiveBayes$Positive;
	d1=data.frame(id=1:nrow(d),upseq=d$upstream.seq,dnseq=d$downstream.seq);
	write.table(file="cleanUpdTSeq.pdata",d1,row.names=F,col.names=T,quote=F,sep="\t")

	d=data.NaiveBayes$Negative;
	d1=data.frame(id=1:nrow(d),upseq=d$upstream.seq,dnseq=d$downstream.seq);
	write.table(file="cleanUpdTSeq.ndata",d1,row.names=F,col.names=T,quote=F,sep="\t")
' 
echo "$cmd" | R --no-save 
}
