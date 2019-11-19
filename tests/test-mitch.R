library("mitch")
library("testthat")

test_that("multiplication works", {
    expect_equal(2 * 2, 4)
})

# 1d
data(rna,genesetsExample)
y<-mitch_import(rna,DEtype="edgeR")
res<-mitch_calc(y,genesetsExample,cores=2)
mitch_plots(res,outfile="1d.pdf") 
if (file.exists("1d.html")) { unlink("1d.html") } 
mitch_report(res,"1d.html")

test_that("1d works", {
    expect_equal(  length(which(res$enrichment_result$p.adjustANOVA<0.1)) ,1)
    expect_true(file.info("1d.pdf")$size>10000)
    expect_true(file.info("1d.html")$size>4000000)
})

unlink("1d.html")
unlink("1d.pdf")


# 2d
data(rna,k9a,genesetsExample)
x<-list("rna"=rna,"k9a"=k9a)
y<-mitch_import(x,DEtype="edgeR")
res<-mitch_calc(y,genesetsExample,cores=2)
mitch_plots(res,outfile="2d.pdf")
if (file.exists("2d.html")) { unlink("2d.html") }
mitch_report(res,"2d.html")

test_that("2d works", {
    expect_equal(  length(which(res$enrichment_result$p.adjustMANOVA<0.1)) ,1)
    expect_true(file.info("2d.pdf")$size>100000)
    expect_true(file.info("2d.html")$size>6000000)
})

unlink("2d.html")
unlink("2d.pdf")

# 3d
data(rna,k9a,k36a,genesetsExample)
x<-list("rna"=rna,"k9a"=k9a,"k36a"=k36a)
y<-mitch_import(x,DEtype="edgeR")
res<-mitch_calc(y,genesetsExample,cores=2)
mitch_plots(res,outfile="3d.pdf")
if (file.exists("3d.html")) { unlink("3d.html") }
mitch_report(res,"3d.html")

test_that("3d works", {
    expect_equal(  length(which(res$enrichment_result$p.adjustMANOVA<0.1)) ,1)
    expect_true(file.info("3d.pdf")$size>100000)
    expect_true(file.info("3d.html")$size>6000000)
})

unlink("3d.html")
unlink("3d.pdf")


