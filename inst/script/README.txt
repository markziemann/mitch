In the inst/extdata/ folder there are some "raw" data files which are used to
demonstrate mitch functions in the examples. The instructions used to create
these data are in the inst/script/example_datagen.R script.

These are the files:
.
├── k36a.tsv 
├── k9a.tsv
├── rna.tsv
└── sample_genesets.gmt

These three tables are edgeR differential expression/binding tables comparing
samples before and after 60 minutes stimulation. The 'rna.tsv' file is for 
gene expression measured by RNA-seq. The 'k9a.tsv' file represents histone 3 
lysine 9 acetylation responses measured with ChIP-seq. The 'k36a.tsv' file
represents histone 3 lysine 36 acetylation responses measured with ChIP-seq. 
These files represent only 1000 genes to keep file sizes small for the 
purposes of the examples. The full files are to be found in another repo:
https://github.com/markziemann/mitch_paper/tree/master/data-raw
