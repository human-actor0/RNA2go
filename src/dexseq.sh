
dexseq.gen(){
usage="$FUNCNAME <out>";
if [ $# -lt 1 ];then echo "$usage"; return; fi
OUT=$1;
cmd='
countData <- matrix( rpois(100, 100), nrow=25 )
condition=rep( c("ctr", "trt"), each=2 ) 
sampleData <- data.frame(condition=condition)
design <- formula( ~ sample + exon + condition:exon )
groupID <- rep(
paste0("gene", 1:10),
each= 100 )
featureID <- rep(
paste0("exon", 1:100),
times= 10 )

colnames(countData) = condition;
res=data.frame(id=featureID, gid=groupID, countData);
write.table(res,file="'$OUT'",col.names=T,row.names=F,quote=F,sep="\t");
'
echo "$cmd" | R --no-save 
}
dexseq.run(){
usage="$FUNCNAME <table> <ctr> <trt> <out>"
if [ $# -lt 4 ];then echo "$usage"; return; fi
	cmd='
	library(DEXSeq);
	tt=read.table("'$1'",header=T);
	ctr.j=grep("'$2'",colnames(tt));
	trt.j=grep("'$3'",colnames(tt));
	countData=cbind( tt[,ctr.j],tt[,trt.j]);
	sampleData <- data.frame( condition=c(
		rep( "ctr", length(ctr.j)), 
		rep( "trt", length(trt.j)) ));
	design <- formula( ~ sample + exon + condition:exon )
	groupID <- tt$gid
	featureID <- tt$id
	dxd=DEXSeqDataSet( countData, sampleData, design, featureID, groupID )
	dxd = estimateSizeFactors( dxd )
	dxd = estimateDispersions( dxd )
	dxd = testForDEU( dxd )
	dxd = estimateExonFoldChanges( dxd )
	dxr <- DEXSeqResults( dxd )
	write.table(dxr,file="'$4'",col.names=T,row.names=F,quote=F,sep="\t");
	'
	echo "$cmd" | R --no-save 
}

dexseq.run.test(){
	dexseq.gen tmp.i
	dexseq.run tmp.i ctr trt tmp.o
	head tmp.o
	rm tmp.*
}
