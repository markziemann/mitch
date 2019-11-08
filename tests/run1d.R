library(mitch)

rna <- read.table(system.file("extdata/rna.tsv", package = "mitch"), header = TRUE)

y <- mitch_import(rna, "edger", geneID="Name")

genesets <- gmt_import(system.file("extdata/sample_genesets.gmt", package = "mitch"))

res <- mitch_calc(y, genesets, resrows = 5, priority = "effect", cores = 2)

mitch_plots(res, outfile = "1dcharts.pdf")

unlink("1dreport.html")

mitch_report(res, "1dreport.html")

unlink("1dreport.html")

unlink("1dcharts.pdf")

detach("package:mitch", unload = TRUE)
