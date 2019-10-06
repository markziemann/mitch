install.packages("..", repos = NULL, type = "source")
library(mitch)

rna <- read.table("rna.tsv", header = TRUE)
rna <- rna[, 1:6]

k9a <- read.table("k9a.tsv", header = TRUE)
k9a <- k9a[, 1:6]

k36a <- read.table("k36a.tsv", header = TRUE)
k36a <- k36a[, 1:6]

myList <- list(rna = rna, k9a = k9a, k36a = k36a)

myImportedData <- mitch_import(myList, "edger", geneIDcol = "Name")
myImportedData <- head(myImportedData, 1000)

rna <- rna[rna$Name %in% rownames(myImportedData), ]
k9a <- k9a[k9a$Name %in% rownames(myImportedData), ]
k36a <- k36a[k36a$Name %in% rownames(myImportedData), ]

myList <- list(rna = rna, k9a = k9a, k36a = k36a)

genesetsExample <- head(gmt_import("ReactomePathways.gmt"), 200)

resExample <- mitch_calc(myImportedData, genesetsExample, resrows = 5, priority = "significance")

save(genesetsExample, file = "../data/genesetsExample.RData", compress = "xz")
save(k36a, file = "../data/k36a.RData", compress = "xz")
save(k9a, file = "../data/k9a.RData", compress = "xz")
save(rna, file = "../data/rna.RData", compress = "xz")
save(myList, file = "../data/myList.RData", compress = "xz")
save(myImportedData, file = "../data/myImportedData.RData", compress = "xz")
save(resExample, file = "../data/resExample.RData", compress = "xz")

write.table(rna, file = "../inst/extdata/rna.tsv", sep = "\t", quote = FALSE)
write.table(k9a, file = "../inst/extdata/k9a.tsv", sep = "\t", quote = FALSE)
write.table(k36a, file = "../inst/extdata/k36a.tsv", sep = "\t", quote = FALSE)

system2("head -200 ReactomePathways.gmt > ../inst/extdata/sample_genesets.gmt")

message("1d")
myImportedData <- mitch_import(rna, "edger", geneIDcol = "Name")
resExample <- mitch_calc(myImportedData, genesetsExample, resrows = 5, priority = "significance")
nrow(resExample$enrichment_result)

message("2d")
myList <- list(rna = rna, k9a = k9a)
myImportedData <- mitch_import(myList, "edger", geneIDcol = "Name")
resExample <- mitch_calc(myImportedData, genesetsExample, resrows = 5, priority = "significance")
nrow(resExample$enrichment_result)

message("3d")
myList <- list(rna = rna, k9a = k9a, k36a = k36a)
myImportedData <- mitch_import(myList, "edger", geneIDcol = "Name")
resExample <- mitch_calc(myImportedData, genesetsExample, resrows = 5, priority = "significance")
head(resExample$enrichment_result)

