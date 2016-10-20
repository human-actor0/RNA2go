
	library(edgeR);
	CTR="^ctr"; TRT="^trt"; 
	C1="c2$"; C2="c1$";
	tt=read.table("stdin",header=T,check.names=F);
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
	#nn=paste("logFC",TRT,"/",CTR,sep="")
	#tt[[nn]]=llh$table$logFC;
	tt[["logCPM"]]=llh$table$logCPM;
	tt[["logFC"]]=llh$table$logFC;
	tt[["pval"]]=llh$table$PValue;
	##res=data.frame(tt, nn=llh$table$logFC, pval=llh$table$PValue)
	## chrom start end logFC pval
	write.table(tt,file="out", col.names=T,quote=F,row.names=F,sep="\t");
	
