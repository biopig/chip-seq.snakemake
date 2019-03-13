#使用genomicfeatures来构建txdb对象，首先就是安装对应的包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicFeatures", version = "3.8")
BiocManager::install("biomaRt", version = "3.8")
library("GenomicFeatures")
library("biomaRt")

#列出可用的数据库的
listMarts()
listMarts(host = "http://plants.ensembl.org")

#查找数据库中所有的植物的物种
mart<-useMart(biomart = "plants_mart",host = "http://plants.ensembl.org")
datasets <- listDatasets(mart)
datasets$dataset

#下载数据
maize_txdb<-makeTxDbFromBiomart(biomart = "plants_mart",dataset = "zmays_eg_gene",host = "http://plants.ensembl.org")

#保存数据
saveDb(maize_txdb, file="maize_v4_2018_11.sqlite")

#载入
maize_txdb <- loadDb("maize_v4_2018_11.sqlite")
