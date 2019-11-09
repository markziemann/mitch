#install.packages("../..",repos = NULL, type="source")
library(mitch)

DATAURL=paste("https://raw.githubusercontent.com/",
"markziemann/mitch_paper/master/data-raw/",sep="")
RNAURL=paste(DATAURL,"rna.tsv",sep="")
rna <- read.table(RNAURL, header=TRUE)
K9AURL=paste(DATAURL,"k9a.tsv",sep="")
k9a <- read.table(K9AURL, header = TRUE)
K36AURL=paste(DATAURL,"k36a.tsv",sep="")
k36a <- read.table(K36AURL, header = TRUE)
myList <- list(rna = rna, k9a = k9a, k36a = k36a)

myImportedData <- mitch_import(myList, "edger")
myImportedData <- head(myImportedData, 1000)

rna <- rna[rownames(rna) %in% rownames(myImportedData), ]
k9a <- k9a[rownames(k9a) %in% rownames(myImportedData), ]
k36a <- k36a[rownames(k36a) %in% rownames(myImportedData), ]

myList <- list(rna = rna, k9a = k9a, k36a = k36a)

REACTOMEURL=paste(DATAURL,"ReactomePathways.gmt",sep="")
download.file(REACTOMEURL,destfile="ReactomePathways.gmt")
genesetsExample <- head(gmt_import(REACTOMEURL), 200)
system("head -200 ReactomePathways.gmt>../../inst/extdata/sample_genesets.gmt")

resExample<-mitch_calc(myImportedData,genesetsExample,
resrows=5,priority="significance")

save(genesetsExample,file="../../data/genesetsExample.RData",compress="xz")
save(k36a, file = "../../data/k36a.RData", compress = "xz")
save(k9a, file = "../../data/k9a.RData", compress = "xz")
save(rna, file = "../../data/rna.RData", compress = "xz")
save(myList, file = "../../data/myList.RData", compress = "xz")
save(myImportedData, file = "../../data/myImportedData.RData",compress="xz")
save(resExample, file = "../../data/resExample.RData", compress = "xz")

write.table(rna, file = "../../inst/extdata/rna.tsv",
sep = "\t", quote = FALSE)
write.table(k9a, file = "../../inst/extdata/k9a.tsv",
sep = "\t", quote = FALSE)
write.table(k36a, file = "../../inst/extdata/k36a.tsv",
sep = "\t", quote = FALSE)

message("1d")
myImportedData <- mitch_import(rna, "edger" )
resExample <- mitch_calc(myImportedData, genesetsExample,
resrows = 5, priority = "significance")
nrow(resExample$enrichment_result)

message("2d")
myList <- list(rna = rna, k9a = k9a)
myImportedData <- mitch_import(myList, "edger" )
resExample <- mitch_calc(myImportedData, genesetsExample,
resrows = 5, priority = "significance")
nrow(resExample$enrichment_result)

message("3d")
myList <- list(rna = rna, k9a = k9a, k36a = k36a)
myImportedData <- mitch_import(myList, "edger" )
resExample <- mitch_calc(myImportedData, genesetsExample,
resrows = 5, priority = "significance")
nrow(resExample$enrichment_result)

unlink("ReactomePathways.gmt")

