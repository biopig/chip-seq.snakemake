#########################################################################
# the sqlite can build follow the example [植物构建txdb对象](https://www.jianshu.com/writer#/notebooks/26834060/notes/37325933)
#
#
#########################################################################

library(ChIPseeker)
library(GenomicFeatures)

Args <- commandArgs(T)

#maize_txdb <- loadDb("../maize_v4_2018_11.sqlite")
maize_txdb <- loadDb(Args[1])
#file_peak <-readPeakFile("../../7_meme_analysis/get_fasta/AG_1.GPS_peaks.bed")
file_peak <- readPeakFile(Args[2])

plot1<-paste(Args[3],"/","result_summary.pdf",sep="")
plot2<-paste(Args[3],"/","result_result.pdf",sep="")
file_result<-paste(Args[3],"/","result_result.txt",sep="")

file_annotation<-annotatePeak(file_peak,tssRegion = c(-3000, 3000), TxDb = maize_txdb)
#ZmARF10_annotation<-annotatePeak(ZmARF10_peak,tssRegion = c(-3000, 3000), TxDb = maize_txdb, addFlankGeneInfo=TRUE, flankDistance=5000)

pdf(plot1, width=12 , height=12)
covplot(file_peak,chrs=c("1","2","3","4","5","6","7","8","9","10"))
plotAvgProf2(file_peak, TxDb=maize_txdb, upstream=3000, downstream=3000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
#peakHeatmap(file_peak, TxDb=maize_txdb, upstream=1000, downstream=1000, color="red")
plotAnnoPie(file_annotation)
plotAnnoBar(file_annotation)
dev.off()

pdf(plot2, width=12 , height=12)
ChIPseeker::upsetplot(file_annotation, vennpie=TRUE)
plotDistToTSS(file_annotation,title="Distribution of transcription factor-binding loci\nrelative to TSS")
dev.off()

#write.table(as.data.frame(file_annotation),file="ZmARF10_annotaton.txt",append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE)
write.table(as.data.frame(file_annotation),file=file_result, append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE)
q(save = "no")
