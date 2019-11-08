library(mitch)

rna <- read.table(system.file("extdata/rna.tsv", package = "mitch"), header = TRUE)
k9a <- read.table(system.file("extdata/k9a.tsv", package = "mitch"), header = TRUE)
k36a <- read.table(system.file("extdata/k36a.tsv", package = "mitch"), header = TRUE)
x <- list(rna = rna, k9a = k9a, k36a = k36a)

y <- mitch_import(x, "edger", geneID="Name")

genesets <- gmt_import(system.file("extdata/sample_genesets.gmt", package = "mitch"))

res <- mitch_calc(y, genesets, resrows = 5, priority = "effect", cores = 2)

mitch_plots(res, outfile = "3dcharts.pdf")

unlink("3dreport.html")

mitch_report(res, "3dreport.html")

unlink("3dreport.html")

unlink("3dcharts.pdf")

detach("package:mitch", unload = TRUE)
