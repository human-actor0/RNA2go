
cleanUpdTSeq.build_feature(){
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

}
