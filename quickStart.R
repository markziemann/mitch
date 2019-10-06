# Might need to install devtools
install.packages("devtools")
library("devtools")

# Install mitch and taucharts from github
devtools::install_github("markziemann/Mitch")
devtools::install_github("hrbrmstr/taucharts")
library("mitch")

# Download a 2 column DGE rank file from an RNA-seq experiment
# Genes can be ranked by directional p-value, logFC or other metric
download.file("https://raw.githubusercontent.com/markziemann/Mitch/master/ex/rna_LGHGvHGHGV_jn.rnk",destfile="rna_LGHGvHGHGV_jn.rnk")

# Import the rank file
x<-as.matrix(read.table("rna_LGHGvHGHGV_jn.rnk",header=TRUE))

# Get Reactome gene sets
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip",destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets<-gmt_import("ReactomePathways.gmt")

# Run the analysis
res<-mitch_calc(x,genesets,resrows=25,priority="effect")
head(res$enrichment_result)
str(res$analysis_metrics)

# Generate a report in html format
mitch_report(res,"myreport.html")

# Generate high res plots
mitch_plots(res,outfile="mycharts.pdf")
