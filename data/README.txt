
MODEL=cleanUpdTSeq.model.rda
LAPLACE=1
TP_DATA=$HMHOME/data/cleanUpdTSeq.pdata
TN_DATA=$HMHOME/data/cleanUpdTSeq.ndata


#hm polyafilter updnseq2fea $TP_DATA > tmp.pos
#hm polyafilter updnseq2fea $TN_DATA > tmp.neg 
#hm polyafilter train tmp.pos tmp.neg $MODEL $LAPLACE
cmd="
hm polyafilter predict tmp.neg $MODEL tmp.neg.pred 
hm polyafilter predict tmp.pos $MODEL tmp.pos.pred 
"; 
#echo "$cmd" | bsub
rm tmp.*


