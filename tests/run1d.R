library(mitch)

rna <- read.table(system.file("extdata/rna.tsv", package = "mitch"), header = TRUE)

y <- mitch_import(rna, "edger", geneIDcol = "Name")

genesets <- gmt_import(system.file("extdata/sample_genesets.gmt", package = "mitch"))

res <- mitch_calc(y, genesets, resrows = 5, priority = "effect", cores = 2)

mitch_plots(res, outfile = "1dcharts.pdf", cores = 2)

unlink("1dreport.html")

mitch_report(res, "1dreport.html")

detach("package:mitch", unload = TRUE)
