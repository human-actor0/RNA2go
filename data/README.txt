
MODEL=cleanUpdTSeq.model.rda
LAPLACE=1
TP_DATA=$HMHOME/data/cleanUpdTSeq.pdata
TN_DATA=$HMHOME/data/cleanUpdTSeq.ndata


#hm polyafilter build_feature $TP_DATA  > pos.fea
#hm polyafilter build_feature $TN_DATA > neg.fea
#hm polyafilter train pos.fea neg.fea $MODEL $LAPLACE
#hm polyafilter predict pos.fea mod.rda pos.fea.pred
#hm polyafilter predict neg.fea mod.rda neg.fea.pred


