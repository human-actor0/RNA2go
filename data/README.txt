
MODEL=cleanUpdTSeq.model.rda
LAPLACE=1
TP_DATA=$HMHOME/data/cleanUpdTSeq.pdata
TN_DATA=$HMHOME/data/cleanUpdTSeq.ndata


#hm polyafilter updnseq2fea $TP_DATA > tmp.pos
#hm polyafilter updnseq2fea $TN_DATA > tmp.neg 
hm polyafilter train tmp.pos tmp.neg $MODEL $LAPLACE
#hm polyafilter predict pos.fea mod.rda pos.fea.pred
#hm polyafilter predict neg.fea mod.rda neg.fea.pred


