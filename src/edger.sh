edger.rep(){
usage(){ echo "
$FUNCNAME <table> <ctr_col> <trt_col> <out> 
 [option]
	nooptions 
"
}
	if [ $# -lt 4 ];then usage;return; fi
	cmd='
	library(edgeR)
	CTR="^'$2'"; TRT="^'$3'"; 
	tt=read.table("'$1'",header=T,check.names=F);
	cn=colnames(tt);
	m=tt[,c(grep(CTR,cn), grep(TRT,cn))];

	m.cn=colnames(m);
	group=rep(1,length(m.cn)); 
	group[ grep(CTR,m.cn)]=1;
	group[ grep(TRT,m.cn)]=2;
	coef=2;
	
	y=DGEList(counts=m,group=factor(group));
	y=calcNormFactors(y);

	group.this=factor(group);
	H1 <- model.matrix(~ group.this )
	y = estimateGLMCommonDisp(y,H1);
	y = estimateGLMTrendedDisp(y,H1);
	y = estimateGLMTagwiseDisp(y,H1);

	fit=glmFit(y$counts, H1, y$tagwise.dispersion,offset=0,prior.count=0)
	llh=glmLRT(fit,coef=coef)
	#nn=paste("logFC",TRT,"/",CTR,sep="")
	tt[["logCPM"]]=llh$table$logCPM;
	tt[["logFC"]]=llh$table$logFC;
	tt[["pval"]]=llh$table$PValue;
	##res=data.frame(tt, nn=llh$table$logFC, pval=llh$table$PValue)
	## chrom start end logFC pval
	write.table(tt,file="'$4'", col.names=T,quote=F,row.names=F,sep="\t");
	' 
	echo "$cmd" | R --no-save
}
edger.rep.test(){
echo \
"x@trt@ctr	trt1	trt2	ctr1	ctr2
a	1	2	10	20
b	4	0	50	0" > tmp.in
edger.rep tmp.in ctr trt tmp.out
head tmp.out
rm tmp.*


}
edger.interact(){
usage(){ echo "
$FUNCNAME <table> <ctr_prefix> <trt_prefix> <count1_suffix> <cout2_suffix> <output> [option]
"
}
	if [ $# -lt 6 ];then usage;return; fi
	cmd='
	library(edgeR);
	CTR="^'$2'"; TRT="^'$3'"; 
	C1="'$4'$"; C2="'$5'$";
	tt=read.table("'$1'",header=T,check.names=F);
	cn=colnames(tt);
	m=tt[,c(grep(CTR,cn), grep(TRT,cn))];

	m.cn=colnames(m);
	group=rep(1,length(m.cn)); 
	group[ grep(CTR,m.cn)]=1;
	group[ grep(TRT,m.cn)]=2;
	event=rep(1,length(m.cn)); event[ grep(C2,m.cn)]=2;
	
	y=DGEList(counts=m,group=factor(group));
	y=calcNormFactors(y);

	event.this=factor(event);
	group.this=factor(group);
	H1 <- model.matrix(~ event.this + group.this + event.this:group.this, data=y$samples )
	H0 <- model.matrix(~ event.this + group.this )
	coef <- (ncol(H0)+1):ncol(H1)
	#y=estimateCommonDisp(y)
	#y=estimateTagwiseDisp(y, trend="movingave")
	y = estimateGLMCommonDisp(y,H1);
	y = estimateGLMTrendedDisp(y,H1);
	y = estimateGLMTagwiseDisp(y,H1);

	fit=glmFit(y$counts, H1, y$tagwise.dispersion,offset=0,prior.count=0)
	llh=glmLRT(fit,coef=coef)
	nn=paste("logFC",TRT,"/",CTR,sep="")
	tt[[nn]]=llh$table$logFC;
	tt[["pval"]]=llh$table$PValue;
	##res=data.frame(tt, nn=llh$table$logFC, pval=llh$table$PValue)
	## chrom start end logFC pval
	write.table(tt,file="'$6'", col.names=T,quote=F,row.names=F,sep="\t");
	' 
	echo "$cmd" > $6.cmd
	R --no-save -f $6.cmd 
}
edger.interact.test(){
echo \
"x@trt@ctr	trt1.c1	trt1.c2	trc.c1	trt2.c1	trt2.c2	ctr1.c1	ctr1.c2
a	1	10	NA	11	10	11	10
b	2	20	NA	22	20	22	20
c	3	30	NA	33	30	33	30" | edger.interact - ctr trt c2 c1 out

}


edger.norep(){
usage="USAGE $FUNCNAME <nput> <BCV> <output> 
	<BCV>: 	0.4: human data, 
		0.1: genetically identical, 
		0.01: technical replicates (edgeRUsersGuide)
"
if [ $# -lt 3 ];then echo "$usage";return; fi
	cmd='O="'$3'";
		library(edgeR);
		tt=read.table("'$1'",header=T);
		bcv = '$2';
		y = DGEList(counts=tt[,2:3],group=c(1,2))
		y = calcNormFactors(y);
		et = exactTest(y, dispersion=bcv^2)	
		write.table(file=O,cbind(tt,et$table),quote=F,sep="\t",row.names=F,col.names=T);
	'
	echo "$cmd" | R --no-save
}
edger.norep.test(){
echo \
"id	A	B
a	2	3
b	4	19" > tmp.in
edger.norep tmp.in 0.2 tmp.out
cat tmp.out
rm tmp.*
}
