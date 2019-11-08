library(mitch)

rna <- read.table(system.file("extdata/rna.tsv", package = "mitch"), header = TRUE)
k9a <- read.table(system.file("extdata/k9a.tsv", package = "mitch"), header = TRUE)
x <- list(rna = rna, k9a = k9a)

y <- mitch_import(x, "edger", geneID="Name")

genesets <- gmt_import(system.file("extdata/sample_genesets.gmt", package = "mitch"))

res <- mitch_calc(y, genesets, resrows = 5, priority = "effect", cores = 2)

mitch_plots(res, outfile = "2dcharts.pdf")

unlink("2dreport.html")

mitch_report(res, "2dreport.html")

unlink("2dreport.html")

unlink("2dcharts.pdf")

detach("package:mitch", unload = TRUE)

