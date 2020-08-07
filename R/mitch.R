#' mitch: An R package for multi-dimensional pathway enrichment analysis
#'
#' mitch is an R package for multi-dimensional enrichment analysis. At it's 
#' heart, it uses a rank-MANOVA based statistical approach to detect sets of 
#' genes that exhibit enrichment in the multidimensional space as compared to 
#' the background. mitch is useful for pathway analysis of profiling studies 
#' with two to or more contrasts, or in studies with multiple omics profiling, 
#' for example proteomic, transcriptomic, epigenomic analysis of the same 
#' samples. mitch is perfectly suited for pathway level differential analysis 
#' of scRNA-seq data.
#'
#' A typical mitch workflow consists of:
#' 1) Import gene sets with gmt_import()
#' 2) Import profiling data with mitch_import()
#' 3) Calculate enrichments with mitch_calc()
#' 4) And generate plots and reports with mitch_plots() and mitch_report()
#'
#' More documentation on the github page https://github.com/markziemann/mitch
#' or with ?<function>, eg: ?mitch_import
#'
#' @docType package
#' @name mitch
#' @examples
#' # Example workflow
#' # Import some gene sets
#' genesetsExample<-gmt_import(system.file('extdata/sample_genesets.gmt', 
#' package = 'mitch'))
#' # Load some edgeR tables (rna, k9a, k36a). 
#' data(rna,k9a,k36a)
#' # Create a list of differential profiles
#' myList<-list('rna'=rna,'k9a'=k9a,'k36a'=k36a)
#' # Import as edgeR table 
#' myImportedData<-mitch_import(myList,DEtype='edger')
#' # Calculate enrichment using MANOVA
#' resExample<-mitch_calc(myImportedData,genesetsExample,priority='effect',
#' resrows=5,cores=2)
#' # Generate some high res plots in PDF format
#' mitch_plots(resExample,outfile='outres.pdf')
#' #' Generate a report of the analysis in HTML format
#' mitch_report(resExample,'outres.html')
NULL

#' @import utils
utils::globalVariables(c("p.adjustMANOVA", "effect", "p.adjustANOVA", "Var2",
    "value", "..density..","dummy_x","dummy_y"))


mapGeneIds <- function(y, z) {
    if (!is.null(attributes(y)$geneTable)) {
        gt <- attributes(y)$geneTable
        col1 <- length(which(z$geneidentifiers %in% gt[, 1]))
        col2 <- length(which(z$geneidentifiers %in% gt[, 2]))
        
        if (col1 + col2 < (nrow(y)/2)) {
            stop("Error it looks as if the Gene IDs in the profile don't match
            the geneTable")
        }
        
        if (col1 > col2) {
            colnames(gt) = c("geneidentifiers", "GeneSymbol")
            z <- merge(gt, z, by = "geneidentifiers")
            z$geneidentifiers = NULL
        } else {
            colnames(gt) = c("GeneSymbol", "geneidentifiers")
            z <- merge(gt, z, by = "geneidentifiers")
            z$geneidentifiers = NULL
        }
        z <- aggregate(. ~ GeneSymbol, z, function(x) {
            sum(as.numeric(as.character(x)))
        })

        colnames(z) = c("geneidentifiers", "y")
    }
    z
}


edger_score <- function(y , geneIDcol = geneIDcol ) {
    
    NCOL = ncol(y)
    if (NCOL < 2) {
        stop("Error: there are <2 columns in the input, 'PValue' and 'logFC' 
        are required ")
    }
    
    PCOL = length(which(names(y) == "PValue"))
    if (PCOL > 1) {
        stop("Error, there is more than 1 column named 'PValue' in the input")
    }
    if (PCOL < 1) {
        stop("Error, there is no column named 'PValue' in the input")
    }
    
    FCCOL = length(which(names(y) == "logFC"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'logFC' in the input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'logFC' in the input")
    }
    
    s <- sign(y$logFC) * -log10(y$PValue)
    
    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


deseq2_score <- function(y , geneIDcol = geneIDcol ) {

    ZCOL = length(which(names(y) == "stat"))
    if (ZCOL > 1) {
        stop("Error, there is more than 1 column named 'stat' in the input")
    }
    if (ZCOL < 1) {
        stop("Error, there is no column named 'stat' in the input")
    }

    s <- y$stat

    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


absseq_score <- function(y, geneIDcol = geneIDcol ) {
    
    NCOL = ncol(y)
    if (NCOL < 2) {
        stop("Error: there are <2 columns in the input, 'pvalue' and 
        'foldChange' are required ")
    }
    
    PCOL = length(which(names(y) == "pvalue"))
    if (PCOL > 1) {
        stop("Error, there is more than 1 column named 'pvalue' in the input")
    }
    if (PCOL < 1) {
        stop("Error, there is no column named 'pvalue' in the input")
    }
    
    FCCOL = length(which(names(y) == "foldChange"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'foldChange' in
        the input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'foldChange' in the input")
    }
    
    s <- sign(y$foldChange) * -log10(y$pvalue)
    
    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


sleuth_score <- function(y , geneIDcol = geneIDcol ) {
    
    NCOL = ncol(y)
    if (NCOL < 2) {
        stop("Error: there are <2 columns in the input, 'pval' and 'b'
        are required ")
    }
    
    PCOL = length(which(names(y) == "pval"))
    if (PCOL > 1) {
        stop("Error, there is more than 1 column named 'pval' in the input")
    }
    if (PCOL < 1) {
        stop("Error, there is no column named 'pval' in the input")
    }
    
    FCCOL = length(which(names(y) == "b"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'b' in the input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'b' in the input")
    }
    
    s <- sign(y$b) * -log10(y$pval)
    
    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


topconfect_score <- function(y , geneIDcol = geneIDcol ) {
    
    FCCOL = length(which(names(y) == "confect"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'confect' in the input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'confect' in the input")
    }
    
    # better to get the sign of fold change from the effect column
    FCCOL = length(which(names(y) == "effect"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'effect' in the input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'effect' in the input")
    }
    
    # there is a problem with topconfects having some NA values
    yy <- y[!is.na(y$effect), ]
    
    pos <- subset(yy, effect > 0)
    neg <- subset(yy, effect < 0)
    pos$mitchrank <- rev(seq(from = 1, to = nrow(pos)))
    neg$mitchrank <- rev(seq(from = -1, to = -nrow(neg)))
    yy <- rbind(pos, neg)
    s <- yy$mitchrank
    
    if (!is.null(attributes(y)$geneIDcol)) {
        g <- yy[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(yy)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


ballgown_score <- function(y , geneIDcol = geneIDcol ) {

    NCOL = ncol(y)
    if (NCOL < 2) {
        stop("Error: there are <2 columns in the input, 'pval' and 'fc' are
        required ")
    }

    PCOL = length(which(names(y) == "pval"))
    if (PCOL > 1) {
        stop("Error, there is more than 1 column named 'pval' in the input")
    }
    if (PCOL < 1) {
        stop("Error, there is no column named 'pval' in the input")
    }

    FCCOL = length(which(names(y) == "fc"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'fc' in the input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'fc' in the input")
    }

    s <- sign(log2(y$fc)) * -log10(y$pval)

    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


noiseq_score <- function(y , geneIDcol = geneIDcol ) {

    ZCOL = length(which(names(y) == "ranking"))
    if (ZCOL > 1) {
        stop("Error, there is more than 1 column named 'ranking' in the input")
    }
    if (ZCOL < 1) {
        stop("Error, there is no column named 'ranking' in the input")
    }

    s <- y$ranking

    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


tcc_score <- function(y , geneIDcol = geneIDcol ) {

    NCOL = ncol(y)
    if (NCOL < 2) {
        stop("Error: there are <2 columns in the input, 'p.value' and 'm.value'
        are required ")
    }

    PCOL = length(which(names(y) == "p.value"))
    if (PCOL > 1) {
        stop("Error, there is more than 1 column named 'p.value' in the input")
    }
    if (PCOL < 1) {
        stop("Error, there is no column named 'p.value' in the input")
    }

    FCCOL = length(which(names(y) == "m.value"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'm.value' in the input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'm.value' in the input")
    }

    s <- sign(y$m.value) * -log10(y$p.value)

    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


deds_score <- function(y , geneIDcol = geneIDcol ) {

    ZCOL = length(which(names(y) == "t"))
    if (ZCOL > 1) {
        stop("Error, there is more than 1 column named 't' in the input")
    }
    if (ZCOL < 1) {
        stop("Error, there is no column named 't' in the input")
    }

    s <- y$t

    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


cuffdiff_score <- function(y , geneIDcol = geneIDcol ) {

    ZCOL = length(which(names(y) == "test_stat"))
    if (ZCOL > 1) {
        stop("Error, there is more than 1 column named 'test_stat' in the
        input")
    }
    if (ZCOL < 1) {
        stop("Error, there is no column named 'test_stat' in the input")
    }

    s <- y$test_stat

    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


seurat_score <- function(y , geneIDcol = geneIDcol ) {
    
    NCOL = ncol(y)
    if (NCOL < 2) {
        stop("Error: there are <2 columns in the input, 'p_val' and 'avg_logFC'
        are required ")
    }
    
    PCOL = length(which(names(y) == "p_val"))
    if (PCOL > 1) {
        stop("Error, there is more than 1 column named 'p_val' in the input")
    }
    if (PCOL < 1) {
        stop("Error, there is no column named 'p_val' in the input")
    }
    
    FCCOL = length(which(names(y) == "avg_logFC"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'avg_logFC' in the
        input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'avg_logFC' in the input")
    }
    
    s <- sign(y$avg_logFC) * -log10(y$p_val)
    
    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    
    z$y[is.infinite(z$y) & z$y < 0] <- min(z$y[!is.infinite(z$y)]) - 0.01
    z$y[is.infinite(z$y) & z$y > 0] <- max(z$y[!is.infinite(z$y)]) + 0.01
    z
}


muscat_score <- function(y , geneIDcol = geneIDcol ) {
    
    NCOL = ncol(y)
    if (NCOL < 2) {
        stop("Error: there are <2 columns in the input, 'p_val' and 'logFC'
        are required ")
    }
    
    PCOL = length(which(names(y) == "p_val"))
    if (PCOL > 1) {
        stop("Error, there is more than 1 column named 'p_val' in the input")
    }
    if (PCOL < 1) {
        stop("Error, there is no column named 'p_val' in the input")
    }
    
    FCCOL = length(which(names(y) == "logFC"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'logFC' in the input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'logFC' in the input")
    }
    
    s <- sign(y$logFC) * -log10(y$p_val)
    
    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}

scde_score <- function(y , geneIDcol = geneIDcol ) {

    ZCOL = length(which(names(y) == "Z"))
    if (ZCOL > 1) {
        stop("Error, there is more than 1 column named 'Z' in the input")
    }
    if (ZCOL < 1) {
        stop("Error, there is no column named 'Z' in the input")
    }

    s <- y$Z

    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


mast_score <- function(y , geneIDcol = geneIDcol ) {

    NCOL = ncol(y)
    if (NCOL < 2) {
        stop("Error: there are <2 columns in the input, 'Pr(>Chisq)' and 'coef'
        are required ")
    }

    PCOL = length(which(names(y) == "Pr(>Chisq)"))
    if (PCOL > 1) {
        stop("Error, there is more than 1 column named 'Pr(>Chisq)' in the
        input")
    }
    if (PCOL < 1) {
        stop("Error, there is no column named 'Pr(>Chisq)' in the input")
    }

    FCCOL = length(which(names(y) == "coef"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'coef' in the input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'coef' in the input")
    }

    # because data.table object  doesn't seem to work
    y<-as.data.frame(y)

    s <- sign(y$coef) * -log10(y[,"Pr(>Chisq)"])

    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


desingle_score <- function(y , geneIDcol = geneIDcol ) {

    NCOL = ncol(y)
    if (NCOL < 2) {
        stop("Error: there are <2 columns in the input, 'pvalue' and 
        'foldChange' are required ")
    }

    PCOL = length(which(names(y) == "pvalue"))
    if (PCOL > 1) {
        stop("Error, there is more than 1 column named 'pvalue' in the input")
    }
    if (PCOL < 1) {
        stop("Error, there is no column named 'pvalue' in the input")
    }

    FCCOL = length(which(names(y) == "foldChange"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'foldChange' in the
        input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'foldChange' in the input")
    }

    s <- sign(log2(y$foldChange)) * -log10(y$pvalue)

    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


dmrcate_score <- function(y , geneIDcol = geneIDcol ) {

    NCOL = ncol(y)
    if (NCOL < 2) {
        stop("Error: there are <2 columns in the input, 'Stouffer' and 
        'meanbetafc' are required ")
    }

    PCOL = length(which(names(y) == "Stouffer"))
    if (PCOL > 1) {
        stop("Error, there is more than 1 column named 'Stouffer' in the
        input")
    }
    if (PCOL < 1) {
        stop("Error, there is no column named 'Stouffer' in the input")
    }

    FCCOL = length(which(names(y) == "meanbetafc"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'meanbetafc' in the
        input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'meanbetafc' in the input")
    }

    y<-as.data.frame(y)
    s <- sign(y$meanbetafc) * -log10(y$Stouffer)

    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


dep_score <- function(y , geneIDcol = geneIDcol ) {

    NCOL = ncol(y)
    if (NCOL < 2) {
        stop("Error: there are <2 columns in the input, '*p.val' and '*ratio' 
        are required ")
    }

    PCOL = length(grep("p.val",names(y)))
    if (PCOL > 1) {
        message("Note, using the leftmost column with '*p.val' in the name")
    }
    if (PCOL < 1) {
        stop("Error, there is no column named '*p.val' in the input")
    }

    FCCOL = length(grep("ratio",names(y)))
    if (FCCOL > 1) {
        message("Note, using the leftmost column with 'ratio' in the name")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'ratio' in the input")
    }

    PCOL = grep("p.val",names(y))[1]

    FCCOL = grep("ratio",names(y))[1]

    s <- sign(y[,FCCOL]) * -log10(y[,PCOL])

    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


msmstests_score <- function(y , geneIDcol = geneIDcol ) {

    NCOL = ncol(y)
    if (NCOL < 2) {
        stop("Error: there are <2 columns in the input, 'p.value' and 'LogFC'
        are required ")
    }

    PCOL = length(which(names(y) == "p.value"))
    if (PCOL > 1) {
        stop("Error, there is more than 1 column named 'p.value' in the input")
    }
    if (PCOL < 1) {
        stop("Error, there is no column named 'p.value' in the input")
    }

    FCCOL = length(which(names(y) == "LogFC"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'LogFC' in the input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'LogFC' in the input")
    }

    s <- sign(y$LogFC) * -log10(y$p.value)

    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


plgem_score <- function(y , geneIDcol = geneIDcol ) {

    LEN = length(y)
    if (LEN < 2) {
        stop("Error: there are <2 items in the input list, '$p.value' and
        '$PLGEM.STN' are required")
    }

    PCOL = length(which(names(y) == "p.value"))
    if (PCOL > 1) {
        stop("Error, there is more than 1 item in the list named 'p.value'")
    }
    if (PCOL < 1) {
        stop("Error, there is no item in the list named 'p.value'")
    }

    FCCOL = length(which(names(y) == "PLGEM.STN"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'PLGEM.STN' in the
        input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'PLGEM.STN' in the input")
    }

    s <- sign(y$PLGEM.STN) * -log10(y$p.value)
    g <- rownames(y$PLGEM.STN)
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


sdams_score <- function(y , geneIDcol = geneIDcol ) {

    LEN = length(y)
    if (LEN < 2) {
        stop("Error: there are <2 items in the input list, '$pv_2part' and
        '$beta' are required")
    }

    PCOL = length(which(names(y) == "pv_2part"))
    if (PCOL > 1) {
        stop("Error, there is more than 1 item in the list named 'pv_2part'")
    }
    if (PCOL < 1) {
        stop("Error, there is no item in the list named 'pv_2part'")
    }

    FCCOL = length(which(names(y) == "beta"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'beta' in the input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'beta' in the input")
    }

    s <- sign(y$beta) * -log10(y$pv_2part)
    g <- y$feat.names
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


diffbind_score <- function(y , geneIDcol = geneIDcol ) {

    NCOL = ncol(y)
    if (NCOL < 2) {
        stop("Error: there are <2 columns in the input, 'p.value' and 'Fold'
        are required ")
    }

    PCOL = length(which(names(y) == "p.value"))
    if (PCOL > 1) {
        stop("Error, there is more than 1 column named 'p.value' in the input")
    }
    if (PCOL < 1) {
        stop("Error, there is no column named 'p.value' in the input")
    }

    FCCOL = length(which(names(y) == "Fold"))
    if (FCCOL > 1) {
        stop("Error, there is more than 1 column named 'Fold' in the input")
    }
    if (FCCOL < 1) {
        stop("Error, there is no column named 'Fold' in the input")
    }

    y <- as.data.frame(y)
    s <- sign(y$Fold) * -log10(y$p.value)

    if (!is.null(attributes(y)$geneIDcol)) {
        g <- y[, attributes(y)$geneIDcol]
    } else {
        g <- rownames(y)
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    z <- mapGeneIds(y, z)
    z
}


preranked_score <- function(y, joinType , geneIDcol = geneIDcol ) {

    if (!is.null(attributes(y)$geneIDcol)) {
        NCOL = ncol(y)
        if (NCOL > 2) {
            stop("Error: there are >2 columns in the input. Your files need to 
            have only the gene ID and rank stat ")
        }
        g <- y[, attributes(y)$geneIDcol]
        s <- y[ , !(names(y) %in% geneIDcol)]
    } else {
        NCOL = ncol(y)
        if (NCOL > 1) {
            stop("Error: there are >1 columns in the input. Should only contain
            the gene ID and rank stat ")
        }
        g <- rownames(y)
        s <- y[,1]
    }
    z <- data.frame(g, s, stringsAsFactors = FALSE)
    colnames(z) <- c("geneidentifiers", "y")
    if ( is.null(joinType) ) {
        joinType <- "inner" 
    }
    if ( joinType == "inner" ) {
        z <- na.omit(z)
    }
    z <- mapGeneIds(y, z)
    z
}


#' mitch_import
#'
#' This function imports differential omics data from common differential tools
#' like edgeR, limma and DESeq2. It calculates a summarised differential
#' expression metric by multiplying the sign of the log fold change by the 
#' -log10 of the p-value. If this behaviour is not desired, mitch_import can be
#' bypassed in favour of another scoring metric.
#' @param x a list of differential expression tables
#' @param DEtype the program that generated the differential expression table
#' Valid options are 'edgeR', 'DESeq2', 'limma', 'ABSSeq', 'Sleuth', 'Seurat',
#' 'topConfects', 'muscat', 'Swish', 'scDE', 'MAST', 'DEsingle', 'ballgown',
#' 'NOIseq', 'TCC', 'DEDS', 'cuffdiff', 'fishpond', 'missMethyl', 'DMRcate',
#' 'DEP', 'msmsTests', 'plgem', 'SDAMS', 'DEqMS', 'DiffBind' and 'prescored'. 
#' Where 'prescored' is a dataframe containing the test statistic and gene ID 
#' (either in rowname or separate column) and nothing else. 'preranked' is an 
#' alias for 'prescored'.
#' @param geneIDcol the column containing gene names. If gene names are 
#' @param joinType the type of join to perform, either 'inner' or 'full'.
#' By default, joins are 'inner' except for Seurat and muscat where full is 
#' used. specified as row names, then geneIDcol=NULL.
#' @param geneTable a 2 column table mapping gene identifiers in the profile to
#' gene identifiers in the gene sets. 
#' @return a multi-column table compatible with mitch_calc analysis.
#' @keywords import mitch
#' @export
#' @examples
#' # first step is to create a list of differential profiles
#' data(rna,k9a,k36a)
#' x<-list('rna'=rna,'k9a'=k9a,'k36a'=k36a)
#' # import as edgeR table 
#' imported<-mitch_import(x,DEtype='edger')
#' @importFrom plyr join_all
mitch_import <- function(x, DEtype, geneIDcol = NULL, geneTable = NULL, 
    joinType = NULL) {
    
    if (is.data.frame(x)) {
        message("The input is a single dataframe; one contrast only. Converting
        it to a list for you.")
        NAME = deparse(substitute(x))
        x <- list(x = x)
    }
    
    if (!is.list(x)) {
        stop("Error: Input (x) must be a LIST of dataframes.")
    }
    
    if (is.null(names(x))) {
        stop("Error: Input (x) must be a NAMED list of dataframes.")
    }

    if ( length(x)>69 ) {
        stop("Error: mitch is currently limited to 69 dimensions or fewer.")
    }
    
    if (!is.null(geneTable) && !is.data.frame(geneTable)) {
        stop("Error: the geneTable needs to be a dataframe.")
    }
    
    if (!is.null(geneTable) && (ncol(geneTable) < 2 || ncol(geneTable) > 2)) {
        stop("Error: the geneTable needs to be a dataframe of 2 columns.")
    }
    
    # the geneIDcol should be an attribute added to each list item
    for (i in seq_len(length(x))) {
        if (!is.null(geneIDcol)) {
            LEN = length(which(names(x[[i]]) %in% geneIDcol))
            if (LEN < 1) {
                stop("Error: the specified geneIDcol doesn't seem to exist")
            }
            if (LEN > 1) {
                stop("Error: there are multiple matches for the specified
                geneIDcol")
            }
            attributes(x[[i]])$geneIDcol <- which(names(x[[i]]) %in% geneIDcol)
        } else {
            attributes(x[[i]])$geneIDcol <- NULL
        }
        if (!is.null(geneTable)) {
            if (!is.data.frame(geneTable)) {
                stop("Error: geneTable is not a data frame.")
            }
            if (ncol(geneTable) != 2) {
                stop("Error: geneTable must be a 2 column dataframe.")
            }
            attributes(x[[i]])$geneTable <- geneTable
        }
    }
    
    DEtype = tolower(DEtype)
    
    validDEtype = c("edger", "deseq2", "limma", "absseq", "sleuth", "seurat",
        "topconfects", "muscat", "swish", "scde", "mast", "desingle",
        "ballgown", "noiseq", "tcc", "deds", "cuffdiff", "preranked", 
        "prescored","fishpond", "missmethyl", "dmrcate", "dep", "msmstests",
        "plgem", "sdams", "deqms", "diffbind")

    if (DEtype == "edger") {
        xx <- lapply(x, edger_score)
    } else if (DEtype == "deseq2" || DEtype == "swish" || 
    DEtype == "fishpond" ) {
        xx <- lapply(x, deseq2_score)
    } else if (DEtype == "absseq") {
        xx <- lapply(x, absseq_score)
    } else if (DEtype == "sleuth") {
        xx <- lapply(x, sleuth_score)
    } else if (DEtype == "seurat") {
        xx <- lapply(x, seurat_score)
    } else if (DEtype == "topconfects") {
        xx <- lapply(x, topconfect_score)
    } else if (DEtype == "muscat") {
        xx <- lapply(x, muscat_score)
    } else if (DEtype == "scde") {
        xx <- lapply(x, scde_score)
    } else if (DEtype == "mast" ) {
        xx <- lapply(x, mast_score)
    } else if (DEtype == "desingle" ) {
        xx <- lapply(x, desingle_score)
    } else if (DEtype == "ballgown" ) {
        xx <- lapply(x, ballgown_score)
    } else if (DEtype == "noiseq" ) {
        xx <- lapply(x, noiseq_score)
    } else if (DEtype == "tcc" ) {
        xx <- lapply(x, tcc_score)
    } else if (DEtype == "deds" || DEtype == "missmethyl" || 
    DEtype == "limma" || DEtype == "deqms" ) {
        xx <- lapply(x, deds_score)
    } else if (DEtype == "cuffdiff" ) {
        xx <- lapply(x, cuffdiff_score)
    } else if (DEtype == "dmrcate" ) {
        xx <- lapply(x, dmrcate_score)
    } else if (DEtype == "dep" ) {
        xx <- lapply(x, dep_score)
    } else if (DEtype == "msmstests" ) {
        xx <- lapply(x, msmstests_score)
    } else if (DEtype == "plgem" ) {
        xx <- lapply(x, plgem_score)
    } else if (DEtype == "sdams" ) {
        xx <- lapply(x, sdams_score)
    } else if (DEtype == "diffbind" ) {
        xx <- lapply(x, diffbind_score)
    } else if (DEtype == "preranked" || DEtype == "prescored") {
        xx <- lapply(x, preranked_score, joinType = joinType)
    } else {
        stop(paste("Specified DEtype does not match one of the following:",
        validDEtype))
    }
    
    # give the colums a unique name otherwise join_all will fail
    for (i in seq_len(length(xx))) {
        colnames(xx[[i]]) <- c("geneidentifiers", paste("y", i, sep = ""))
    }
    
    if (is.null(joinType)) {
        if (DEtype == "seurat" || DEtype == "muscat" || DEtype == "scde" ||
        DEtype == "mast" || DEtype == "desingle" ) {
            xxx <- join_all(xx, by = "geneidentifiers", type = "full")
        } else {
            xxx <- join_all(xx, by = "geneidentifiers", type = "inner")
        }
    } else {
        xxx <- join_all(xx, by = "geneidentifiers", type = joinType)
    }
    rownames(xxx) <- xxx$geneidentifiers
    xxx$geneidentifiers = NULL
    colnames(xxx) <- names(x)
    
    STARTSWITHNUM = length(grep("^[0-9]", colnames(xxx)))
    if (STARTSWITHNUM > 0) {
        stop("Error: it looks like one or more column names starts with a
        number. This is incompatible with downstream analysis. Please modify")
    }
    
    MEAN_N_GENES_IN = mean(unlist(lapply(x, nrow)))
    N_GENES_OUT = nrow(xxx)
    PROP = signif(N_GENES_OUT/MEAN_N_GENES_IN, 3)
    message(paste("Note: Mean no. genes in input =", MEAN_N_GENES_IN))
    message(paste("Note: no. genes in output =", N_GENES_OUT))
    if (PROP < 0.05) {
        warning("Warning: less than half of the input genes are also in the
        output")
    } else {
        message(paste("Note: estimated proportion of input genes in output =",
        PROP))
    }
    return(xxx)
}

#' gmt_import
#'
#' This function imports GMT files into a list of character vectors for mitch 
#' analysis. GMT files are a commonly used
#' format for lists of genes used in pathway enrichment analysis. GMT files can
#' be obtained from Reactome, MSigDB, etc.
#' @param gmtfile a gmt file.
#' @return a list of gene sets.
#' @keywords import genesets
#' @export
#' @examples
#' # Import some gene sets
#' genesetsExample<-gmt_import(system.file('extdata/sample_genesets.gmt', 
#' package = 'mitch'))
#' @import utils
gmt_import <- function(gmtfile) {
    genesetLines <- strsplit(readLines(gmtfile), "\t")
    genesets <- lapply(genesetLines, utils::tail, -2)
    names(genesets) <- unlist(lapply(genesetLines, head, 1))
    attributes(genesets)$originfile <- gmtfile
    genesets
}

MANOVA <- function(x, genesets, minsetsize = 10, cores = detectCores() - 1,
    priority = NULL) {
    
    STARTSWITHNUM = length(grep("^[0-9]", colnames(x)))
    if (STARTSWITHNUM > 0) {
        stop("Error: it looks like one or more column names starts with a
        number. This is incompatible with downstream analysis. Please modify")
    }
    
    sets <- names(genesets)
    
    if (is.null(priority)) {
        priority = "significance"
    }
    
    if (priority !="significance" && priority !="effect" && priority !="SD"){
        stop("Error: Parameter 'priority' must be either 'significance'(the 
        default), 'effect' or 'SD'.")
    }
    
    hypotenuse <- function(x) {
        sqrt(sum(unlist(lapply(x, function(x) { x^2 })),na.rm = TRUE ))
    }
    
    # calculate the hypotenuse for downstream use
    HYPOT = hypotenuse(apply(x, 2, length))
    
    res <- pbmclapply(sets, function(set) {
        
        inset <- rownames(x) %in% as.character(unlist(genesets[set]))
        
        NROW = nrow(x)
        
        if (length(which(inset)) >= minsetsize) {
            fit <- manova(x ~ inset)
            sumMANOVA <- summary.manova(fit)
            sumAOV <- summary.aov(fit)
            pMANOVA <- sumMANOVA$stats[1, "Pr(>F)"]
            raov <- lapply(sumAOV, function(zz) {
                zz[1, "Pr(>F)"]
            })
            raov <- unlist(raov)
            names(raov) <- gsub("^ Response ", "p.", names(raov))
            # S coordinates
            NOTINSET <- colMeans(x[!inset, ], na.rm=TRUE)
            scord <- (2 * (colMeans(x[inset, ], na.rm=TRUE) - NOTINSET))/NROW
            names(scord) <- paste0("s-", names(scord))
            # calculate the hypotenuse length of s scores
            s.dist <- hypotenuse(scord)
            names(s.dist) = "s.dist"
            mysd = sd(scord)
            names(mysd) = "SD"
            
            return(data.frame(set, setSize = sum(inset), pMANOVA, t(scord), 
                t(raov),t(s.dist), t(mysd), stringsAsFactors = FALSE))
        }
    }, mc.cores = cores)
    
    fres <- ldply(res, data.frame)
    
    if (nrow(fres) < 1) {
        
        message("Warning: No results found. Check that the gene names in the
            profile match the gene sets and consider loosening the minsetsize
            parameter.")
        
    } else {
        fres$p.adjustMANOVA <- p.adjust(fres$pMANOVA, "fdr")
        
        # prioritisation
        if (priority == "significance") {
            fres <- fres[order(fres$pMANOVA), ]
            message("Note: When prioritising by significance (ie: small 
            p-values), large effect sizes might be missed.")
        }
        if (priority == "effect") {
            fres <- fres[order(-fres$s.dist), ]
            message("Note: Enrichments with large effect sizes may not be 
            statistically significant.")
        }
        if (priority == "SD") {
            fres <- fres[order(-fres$SD), ]
            fres <- subset(fres, p.adjustMANOVA <= 0.05)
            message("Note: Prioritisation by SD after selecting sets with 
            p.adjustMANOVA<=0.05.")
        }
        attributes(fres)$priority <- priority
        return(fres)
    }
}


ANOVA <- function(x, genesets, minsetsize = 10, cores = detectCores() - 1, 
priority = NULL) {
    
    STARTSWITHNUM = length(grep("^[0-9]", colnames(x)))
    if (STARTSWITHNUM > 0) {
        stop("Error: it looks like one or more column names starts with a
        number. This is incompatible with downstream analysis. Please modify")
    }
    
    sets <- names(genesets)
    
    if (is.null(priority)) {
        priority = "significance"
    }
    
    if (priority != "significance" && priority != "effect") {
        stop("Error: Parameter 'priority' must be either 'significance'
        (default) or 'effect'.")
    }

    x<-x[which(!is.na(x)),,drop=FALSE]
    
    res <- pbmclapply(sets, function(set) {
        resample <- function(x, set) {
            sss <- x[which(rownames(x) %in% 
                as.character(unlist(genesets[set]))),]
            mysample <- sample(sss, length(sss), replace = TRUE)
            mean(mysample)
        }
        
        inset <- rownames(x) %in% as.character(unlist(genesets[set]))
        
        NROW = nrow(x)
        
        if (length(which(inset)) >= minsetsize) {
            fit <- aov(x[, 1] ~ inset)
            pANOVA <- summary(fit)[[1]][, 5][1]
            NOTINSET <- mean(x[!inset, ])
            s.dist <- (2 * (mean(x[inset, ]) - NOTINSET))/NROW
            gres <- data.frame(set, setSize = sum(inset), pANOVA, s.dist,
            stringsAsFactors = FALSE)
            gres
        }
    }, mc.cores = cores)
    
    fres <- ldply(res, data.frame)
    
    if (nrow(fres) < 1) {
        
        message("Warning: No results found. Check that the gene names in the
        profile match the gene sets and consider loosening the minsetsize
        parameter.")
        
    } else {
        fres$p.adjustANOVA <- p.adjust(fres$pANOVA, "fdr")
        
        # prioritisation
        if (priority == "significance") {
            fres <- fres[order(fres$pANOVA), ]
            message("Note: When prioritising by significance (ie: small
            p-values), large effect sizes might be missed.")
        }
        if (priority == "effect") {
            fres <- fres[order(-abs(fres$s.dist)), ]
            message("Note: Enrichments with large effect sizes may not be
            statistically significant.")
        }
        attributes(fres)$priority <- priority
        return(fres)
    }
}

mitch_metrics_calc<-function(x, genesets, enrichment_result, minsetsize = 10) {
    
    if (!is.null(enrichment_result)) {
        
        num_genesets = length(genesets)
        included_genesets <- nrow(enrichment_result)
        geneset_counts <- as.data.frame(as.vector(unlist(lapply(genesets, 
        function(set) {
            length(which(as.vector(unlist(set)) %in% rownames(x)))
        }))))
        rownames(geneset_counts) <- names(genesets)
        colnames(geneset_counts) = "count"
        genesets_excluded=
            names(genesets)[which(geneset_counts$count<minsetsize)]
        genesets_included=
            names(genesets)[which(geneset_counts$count>=minsetsize)]
        num_genesets_excluded = length(genesets_excluded)
        num_genesets_included = length(genesets_included)
        num_genes_in_genesets = length(unique(as.vector(unlist(genesets))))
        num_genes_in_profile = length(unique(rownames(x)))
        duplicated_genes_present = length(rownames(x)) > num_genes_in_profile
        num_profile_genes_in_sets=length(which(rownames(x) %in% 
        as.vector(unlist(genesets))))
        num_profile_genes_not_in_sets = num_genes_in_profile - 
        num_profile_genes_in_sets
        num_sets_significant = 
            nrow(enrichment_result[which(enrichment_result$p.adjustMANOVA < 
            0.05), ])
        profile_pearson_correl = cor(x, method = "p")[2, 1]
        profile_spearman_correl = cor(x, method = "s")[2, 1]
        
        # genes in each quadrant
        g1 = length(which(x[, 1] > 0 & x[, 2] > 0))
        g2 = length(which(x[, 1] > 0 & x[, 2] < 0))
        g3 = length(which(x[, 1] < 0 & x[, 2] < 0))
        g2 = length(which(x[, 1] < 0 & x[, 2] > 0))
        
        # genesets in each quadrant
        ns1 = nrow(subset(enrichment_result, p.adjustMANOVA < 0.05 &
            enrichment_result[,4] > 0 & enrichment_result[, 5] > 0))
        ns2 = nrow(subset(enrichment_result, p.adjustMANOVA < 0.05 & 
            enrichment_result[,4] > 0 & enrichment_result[, 5] < 0))
        ns3 = nrow(subset(enrichment_result, p.adjustMANOVA < 0.05 &
            enrichment_result[,4] < 0 & enrichment_result[, 5] < 0))
        ns4 = nrow(subset(enrichment_result, p.adjustMANOVA < 0.05 &
            enrichment_result[,4] < 0 & enrichment_result[, 5] > 0))
        num_sets_significant_by_quadrant = paste(ns1, ns2, ns3, ns4, sep = ",")
        
        dat <- list(num_genesets = num_genesets, 
            num_genes_in_profile = num_genes_in_profile, 
            duplicated_genes_present = duplicated_genes_present, 
            num_profile_genes_in_sets = num_profile_genes_in_sets, 
            num_profile_genes_not_in_sets = num_profile_genes_not_in_sets,
            num_genesets_excluded = num_genesets_excluded, 
            num_genesets_included = num_genesets_included, 
            num_genes_in_genesets = num_genes_in_genesets, 
            genesets_excluded = genesets_excluded,
            genesets_included = genesets_included, 
            profile_pearson_correl = profile_pearson_correl,
            profile_spearman_correl = profile_spearman_correl, 
            num_sets_significant = num_sets_significant,
            num_sets_significant_by_quadrant=num_sets_significant_by_quadrant,
            geneset_counts = geneset_counts)
        dat
    }
}

mitch_metrics_calc1d <- function(x, genesets, anova_result, minsetsize = 10) {
    
    if (!is.null(anova_result)) {
        
        num_genesets = length(genesets)
        included_genesets <- nrow(anova_result)
        geneset_counts <- as.data.frame(as.vector(unlist(lapply(genesets,
            function(set) {
            length(which(as.vector(unlist(set)) %in% rownames(x)))
        }))))
        rownames(geneset_counts) <- names(genesets)
        colnames(geneset_counts) = "count"
        genesets_excluded=
            names(genesets)[which(geneset_counts$count<minsetsize)]
        genesets_included=
            names(genesets)[which(geneset_counts$count >= minsetsize)]
        num_genesets_excluded = length(genesets_excluded)
        num_genesets_included = length(genesets_included)
        num_genes_in_genesets = length(unique(as.vector(unlist(genesets))))
        num_genes_in_profile = length(unique(rownames(x)))
        duplicated_genes_present = length(rownames(x)) > num_genes_in_profile
        num_profile_genes_in_sets = length(which(rownames(x) %in% 
            as.vector(unlist(genesets))))
        num_profile_genes_not_in_sets=num_genes_in_profile-num_profile_genes_in_sets
        num_sets_significant=
            nrow(anova_result[which(anova_result$p.adjustANOVA<0.05),])
        
        # genes up and down
        g1 = length(which(x[, 1] > 0))
        g2 = length(which(x[, 1] < 0))
        
        # genesets in each quadrant
        num_sets_up = nrow(subset(anova_result, p.adjustANOVA < 0.05 &
            anova_result[,4] > 0))
        num_sets_dn = nrow(subset(anova_result, p.adjustANOVA < 0.05 &
            anova_result[,4] < 0))
        
        dat <- list(num_genesets = num_genesets,
            num_genes_in_profile = num_genes_in_profile, 
            duplicated_genes_present = duplicated_genes_present,
            num_profile_genes_in_sets = num_profile_genes_in_sets, 
            num_profile_genes_not_in_sets = num_profile_genes_not_in_sets,
            num_genesets_excluded = num_genesets_excluded, 
            num_genesets_included = num_genesets_included,
            num_genes_in_genesets = num_genes_in_genesets, 
            genesets_excluded = genesets_excluded,
            genesets_included = genesets_included, 
            num_sets_significant = num_sets_significant,
            num_sets_up = num_sets_up, 
            num_sets_dn = num_sets_dn,
            geneset_counts = geneset_counts)
        dat
    }
}


mitch_rank <- function(x) {
    
    for (i in seq_len(ncol(x))) {
        LEN = length(x[, i])
        UNIQLEN = length(unique(x[, i]))
        if (UNIQLEN/LEN < 0.4) {
            warning("Warning: >60% of genes have the same score. This isn't
            optimal for rank based enrichment analysis.")
        }
    }
    
    rank_adj <- function(x) {
        xx <- rank(x,na.last = "keep")
        num_neg = length(which(x < 0))
        num_zero = length(which(x == 0))
        num_adj = num_neg + (num_zero/2)
        adj <- xx - num_adj
        adj
    }
    adj <- apply(x, 2, rank_adj)
    adj
}


detailed_sets <- function(res, resrows = 50) {
    # collect ranked genelist of each genest
    genesets <- res$input_genesets
    ss <- res$ranked_profile
    mykeys <- as.character(res$enrichment_result[seq_len(resrows), 1])
    dat <- vector(mode = "list", length = resrows)
    names(dat) <- mykeys
    
    for (i in seq_len(resrows)) {
        sss <- ss[which(rownames(ss) %in% genesets[[which(names(genesets) %in%
        as.character(res$enrichment_result[i, 
            1]))]]), ]
        dat[[i]] <- sss
    }
    dat
}


get_os <- function(){
    sysinf <- Sys.info()
    if (!is.null(sysinf)){
        os <- sysinf['sysname']
        if (os == 'Darwin')
            os <- "osx"
    } else { ## mystery machine
        os <- .Platform$OS.type
        if (grepl("^darwin", R.version$os))
            os <- "osx"
        if (grepl("linux-gnu", R.version$os))
            os <- "linux"
    }
    tolower(os)
}


#' mitch_calc
#'
#' This function performs multivariate gene set enrichment analysis. 
#' @param x a multicolumn numerical table with each column containing 
#' differential expression scores for a contrast.
#' Rownames must match genesets.
#' @param genesets lists of genes imported by the gmt_imprt function or
#' similar.
#' @param minsetsize the minimum number of genes required in a set for it to be
#' included in the statistical analysis.
#' Default is 10.
#' @param cores the number of parallel threads for computation. Defaults to the 
#' number of cores present minus 1.
#' @param resrows an integer value representing the number of top genesets for
#' which a detailed report is to be 
#' generated. Default is 50.
#' @param priority the prioritisation metric used to selecting top gene sets. 
#' Valid options are 'significance', 
#' 'effect' and 'SD'. 
#' @return mitch res object with the following parts:
#' $input_profile: the supplied input differential profile
#' $input_genesets: the supplied input gene sets
#' $ranked_profile: the differential profile after ranking
#' $enrichment_result: the table of MANOVA/ANOVA enrichment results for each
#' gene set
#' $analysis_metrics:  several metrics that are important to the interpretation
#' of the results
#' $detailed_sets: a list of dataframes containing ranks of members of
#' prioritised gene sets.
#' @keywords mitch calc calculate manova 
#' @import parallel
#' @importFrom pbmcapply pbmclapply
#' @import stats
#' @importFrom plyr ldply
#' @export
#' @examples
#' # Example using mitch to calculate multivariate enrichments and
#' # prioritise based on effect size 
#' data(myImportedData,genesetsExample)
#' resExample<-mitch_calc(myImportedData,genesetsExample,priority='effect',
#' minsetsize=5,cores=2)
mitch_calc <- function(x, genesets, minsetsize = 10, cores = detectCores() - 1,
    resrows = 50, priority = NULL) {
    
    colnames(x) <- sub("-", "_", colnames(x))
    input_profile <- x
    input_genesets <- genesets
    ranked_profile <- mitch_rank(input_profile)
    if (get_os() == "windows") { cores=1 }
    
    if ( ncol(x)>69 ) {
        stop("Error: mitch is currently limited to 69 dimensions or fewer.")
    }

    if (ncol(x) > 1) {
        enrichment_result <- MANOVA(ranked_profile, genesets, 
            minsetsize = minsetsize, cores = cores, priority = priority)
        
        if (!is.null(enrichment_result)) {
            mitch_metrics <- mitch_metrics_calc(x, genesets, enrichment_result)
            dat <- list(input_profile = input_profile,
                input_genesets = input_genesets,
                ranked_profile = ranked_profile,
                enrichment_result = enrichment_result, 
                analysis_metrics = mitch_metrics)
            
            if (nrow(enrichment_result) < resrows) {
                resrows <- nrow(enrichment_result)
            }
            dat$detailed_sets <- detailed_sets(dat, resrows)
            attr(dat, "profile_dimensions") <- colnames(dat$input_profile)
            dat
        }
    } else if (ncol(x) == 1) {
        
        enrichment_result <- ANOVA(ranked_profile, genesets,
            minsetsize = minsetsize, cores = cores, priority = priority)
        
        if (!is.null(enrichment_result)) {
            mitch_metrics <- mitch_metrics_calc1d(x, genesets,
            enrichment_result)
            dat <- list(input_profile = input_profile,
                input_genesets = input_genesets, 
                ranked_profile = ranked_profile,
                enrichment_result = enrichment_result, 
                analysis_metrics = mitch_metrics)
            
            if (nrow(enrichment_result) < resrows) {
                resrows <- nrow(enrichment_result)
            }
            dat$detailed_sets <- detailed_sets(dat, resrows)
            attr(dat, "profile_dimensions") <- colnames(dat$input_profile)
            dat
        }
    }
}

plot1d_profile_dist <- function(res) {
    par(mfrow = c(2, 1))
    hist(res$input_profile[, 1], breaks = 50,
        main = "Distribution of DE scores", xlab = paste("DE score for ",
        colnames(res$input_profile)))
    plot(res$input_profile, xlab = paste("DE score for ",
        colnames(res$input_profile)), 
        pch = "|", frame.plot = FALSE)
    UPS = length(which(res$input_profile > 0))
    DNS = length(which(res$input_profile < 0))
    TOTAL = nrow(res$input_profile)
    mtext(paste(TOTAL, "genes in total,", UPS, "trending up-regulated,", 
    DNS, "trending down-regulated"))
    pl <- recordPlot()
    pl
}

plot_geneset_hist <- function(res) {
    par(mfrow = c(3, 1))
    geneset_counts <- res$analysis_metrics$geneset_counts
    boxplot(geneset_counts$count, horizontal = TRUE, frame = FALSE, 
        main = "Gene set size", 
        xlab = "number of member genes included in profile")
    hist(geneset_counts$count, 100, xlab = "geneset size",
        main = "Histogram of geneset size")
    hist(geneset_counts$count, 100, xlim = c(0, 500), xlab = "geneset size",
        main = "Trimmed histogram of geneset size")
    pl <- recordPlot()
    pl
}

plot1d_volcano <- function(res) {
    par(mfrow = c(1, 1))
    sig <- subset(res$enrichment_result, p.adjustANOVA <= 0.05)
    plot(res$enrichment_result$s.dist, -log10(res$enrichment_result$pANOVA),
    xlab = "s score", 
        ylab = "-log10(p-value)",
        main = "volcano plot of gene set enrichments",
        pch = 19, cex = 0.8)
    points(sig$s.dist, -log10(sig$pANOVA), pch = 19, cex = 0.85, col = "red")
    TOTAL = nrow(res$enrichment_result)
    SIG = nrow(sig)
    UP = length(which(sig$s.dist > 0))
    DN = length(which(sig$s.dist < 0))
    SUBHEADER = paste(TOTAL, "gene sets in total,", UP, "upregulated and ",
        DN, "downregulated (FDR<=0.05)")
    mtext(SUBHEADER)
    pl <- recordPlot()
    pl
}

plot1d_detailed <- function(res, i) {
    par(mfrow = c(3, 1))
    ss <- res$ranked_profile
    sss <- res$detailed_sets[[i]]
    set <- names(res$detailed_sets[i])
    size <- length(sss)
    
    beeswarm(sss, vertical = FALSE, cex = 0.75, xlim = c(min(ss), max(ss)),
        col = "darkgray", pch = 19, main = set, cex.main = 1.5,
        xlab = paste("ranked DE score in:", colnames(ss)))
    mtext("beeswarm plot", cex = 0.8)
    
    hist(sss, xlim = c(min(ss), max(ss)), breaks = 15, col = "darkgray", 
        main = NULL, border = "black", 
        xlab = paste("ranked DE score in:", colnames(ss)))
    mtext("histogram", cex = 0.8)
    
    plot(sss, rep(1, length(sss)), type = "n", xlim = c(min(ss), max(ss)), 
        frame = FALSE, axes = FALSE, ylab = "", 
        xlab = paste("ranked DE score in:", colnames(ss)))
    rug(sss, ticksize = 0.9)
    axis(1)
    mtext("rugplot", cex = 0.8)
    pl <- recordPlot()
    pl
}

plot2d_profile_dist <- function(res) {
    plot(res$input_profile, pch = 19, 
        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2), 
        main = "Scatterplot of all genes")
    abline(v = 0, h = 0, lty = 2, lwd = 2, col = "blue")
    pl <- recordPlot()
    pl
}


plot2d_profile_density <- function(res) {
    palette <- colorRampPalette(c("white", "yellow", "orange", "red", 
        "darkred", "black"))
    
    ss <- res$ranked_profile
    xmin = min(ss[, 1])
    xmax = max(ss[, 1])
    ymin = min(ss[, 2])
    ymax = max(ss[, 2])
    
    k <- MASS::kde2d(ss[, 1], ss[, 2])
    X_AXIS = paste("Rank in contrast", colnames(ss)[1])
    Y_AXIS = paste("Rank in contrast", colnames(ss)[2])
    
    filled.contour(k, xlim = c(xmin, xmax), ylim = c(ymin, ymax), 
        color.palette = palette, plot.title = {
            abline(v = 0, h = 0, lty = 2, lwd = 2, col = "blue")
            title(main = "Rank-rank plot of all genes", xlab = X_AXIS, 
            ylab = Y_AXIS)
        })
    pl <- recordPlot()
    pl
}

plot2d_gene_quadrant_barchart <- function(res) {
    uu = length(which(res$input_profile[, 1] > 0 & res$input_profile[, 2] > 0))
    ud = length(which(res$input_profile[, 1] > 0 & res$input_profile[, 2] < 0))
    dd = length(which(res$input_profile[, 1] < 0 & res$input_profile[, 2] < 0))
    du = length(which(res$input_profile[, 1] < 0 & res$input_profile[, 2] > 0))
    a <- as.data.frame(c(uu, ud, dd, du))
    rownames(a) = c("top-right", "bottom-right", "bottom-left", "top-left")
    colnames(a) = "a"
    barplot(a$a, names.arg = rownames(a), 
    main = "number of genes in each quadrant")
    pl <- recordPlot()
    pl
}

plot2d_set_quadrant_barchart <- function(res) {
    par(mfrow = c(1, 1))
    a <- res$analysis_metrics[14]
    a <- as.data.frame(as.numeric(unlist(strsplit(as.character(a), ","))), 
    stringsAsFactors = FALSE)
    rownames(a) = c("top-right", "bottom-right", "bottom-left", "top-left")
    colnames(a) = "a"
    barplot(a$a, names.arg = rownames(a), main = "number of genesets FDR<0.05")
    pl <- recordPlot()
    pl
}


plot2d_set_scatter <- function(res) {
    sig <- subset(res$enrichment_result, p.adjustMANOVA < 0.05)
    plot(res$enrichment_result[, 4:5], pch = 19, col = rgb(red = 0, green = 0, 
        blue = 0, alpha = 0.2),
        main = "Scatterplot of all gene sets; FDR<0.05 in red")
    abline(v = 0, h = 0, lty = 2, lwd = 2, col = "blue")
    points(sig[, 4:5], pch = 19, col = rgb(red = 1, green = 0, blue = 0,
    alpha = 0.5))
    pl <- recordPlot()
    pl
}

plot2d_set_scatter_top <- function(res) {
    resrows = length(res$detailed_sets)
    top <- head(res$enrichment_result, resrows)
    plot(res$enrichment_result[, 4:5], pch = 19, col = rgb(red = 0, green = 0,
        blue = 0, alpha = 0.2),
        main = paste("Scatterplot of all gene sets; top", resrows, 
        "in red"))
    abline(v = 0, h = 0, lty = 2, lwd = 2, col = "blue")
    points(top[, 4:5], pch = 19, col = rgb(red = 1, green = 0, blue = 0,
    alpha = 0.5))
    pl <- recordPlot()
    pl
}

plot2d_heatmap <- function(res) {
    d = ncol(res$input_profile)
    resrows = length(res$detailed_sets)
    pl = NULL
    if (resrows > 2) {
        hmapx <- head(res$enrichment_result[, 4:(4 + d - 1)], resrows)
        rownames(hmapx) <- head(res$enrichment_result$set, resrows)
        colnames(hmapx) <- gsub("^s.", "", colnames(hmapx))
        my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
        heatmap.2(as.matrix(hmapx), scale = "none", margins = c(10, 25),
            cexRow = 0.8, trace = "none", cexCol = 0.8, col = my_palette)
        pl <- recordPlot()
    }
    pl
}

plot_effect_vs_significance <- function(res) {
    par(mfrow = c(1, 1))
    plot(res$enrichment_result$s.dist,
        -log(res$enrichment_result$p.adjustMANOVA), 
        xlab = "s.dist (effect size)", 
        ylab = "-log(p.adjustMANOVA) (significance)", 
        pch = 19, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2),
        main = "effect size versus statistical significance")
    pl <- recordPlot()
    pl
}


plot2d_detailed_density <- function(res, i) {
    palette <- colorRampPalette(c("white", "yellow", "orange", "red",
        "darkred", "black"))
    ss <- res$ranked_profile
    xmin = min(ss[, 1])
    xmax = max(ss[, 1])
    ymin = min(ss[, 2])
    ymax = max(ss[, 2])
    ll <- res$enrichment_result[i, ]
    size <- ll$setSize
    sss <- res$detailed_sets[[i]]
    X_AXIS = paste("Rank in contrast", colnames(ss)[1])
    Y_AXIS = paste("Rank in contrast", colnames(ss)[2])
    par(mar = c(5, 4, 4, 2))
    k <- MASS::kde2d(sss[, 1], sss[, 2])
    filled.contour(k, color.palette = palette, xlim = c(xmin, xmax),
        ylim = c(ymin, ymax), plot.title = {
        abline(v = 0, h = 0, lty = 2, lwd = 2, col = "blue")
        title(main = ll$set, xlab = X_AXIS, ylab = Y_AXIS)
    })
    pl <- recordPlot()
    pl
}


plot2d_detailed_scatter <- function(res, i) {
    ss <- res$ranked_profile
    xmin = min(ss[, 1])
    xmax = max(ss[, 1])
    ymin = min(ss[, 2])
    ymax = max(ss[, 2])
    sss <- res$detailed_sets[[i]]
    X_AXIS = paste("Rank in contrast", colnames(ss)[1])
    Y_AXIS = paste("Rank in contrast", colnames(ss)[2])
    ll <- res$enrichment_result[i, ]
    plot(sss, pch = 19, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2),
        main = ll$set, xlim = c(xmin, xmax), ylim = c(ymin, ymax), 
        xlab = X_AXIS, ylab = Y_AXIS)
    abline(v = 0, h = 0, lty = 2, lwd = 2, col = "blue")
    pl <- recordPlot()
    pl
}


plot2d_detailed_violin <- function(res, i) {
    pl <- list()
    
    ss <- res$ranked_profile
    xmin = min(ss[, 1])
    xmax = max(ss[, 1])
    ymin = min(ss[, 2])
    ymax = max(ss[, 2])
    ll <- res$enrichment_result[i, ]
    size <- ll$setSize
    sss <- res$detailed_sets[[i]]
    X_AXIS = paste("Rank in contrast", colnames(ss)[1])
    Y_AXIS = paste("Rank in contrast", colnames(ss)[2])
    ss_long <- melt(ss)
    sss_long <- melt(sss)
    
    p <- ggplot(ss_long, aes(Var2, value)) + geom_violin(data = ss_long,
        fill = "grey", colour = "grey") + geom_boxplot(data = ss_long, 
        width = 0.9, fill = "grey", outlier.shape = NA, 
        coef = 0) + geom_violin(data = sss_long, fill = "black", 
        colour = "black") + geom_boxplot(data = sss_long, width = 0.1, 
        outlier.shape = NA) + labs(y = "Position in rank", title = ll[, 1])
    
    print(p + theme_bw() + theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 15), 
        plot.title = element_text(size = 20)))
    pl <- recordPlot()
    pl
}


ggpairs_points <- function(res) {
    ggpairs_points_plot <- function(data, mapping, ...) {
        p <- ggplot(data = data, mapping = mapping) + 
            geom_point(alpha = 0.05) + 
            geom_vline(xintercept = 0, linetype = "dashed") +
            geom_hline(yintercept = 0, 
            linetype = "dashed")
    }
    
    p <- ggpairs(as.data.frame(res$input_profile),
        title = "Scatterplot of all genes", 
        lower = list(continuous = ggpairs_points_plot))
    print(p + theme_bw())
}


ggpairs_points_subset <- function(res) {
    d <- ncol(res$ranked_profile)
    ggpairs_points_plot <- function(data, mapping, ...) {
        p <- ggplot(data = data, mapping = mapping) +
            geom_point(alpha = 0.05) + 
            geom_vline(xintercept = 0, linetype = "dashed") +
            geom_hline(yintercept = 0, 
            linetype = "dashed")
    }
    enrichment_result_clipped <- res$enrichment_result[, 4:(3 + d)]
    colnames(enrichment_result_clipped) <- colnames(res$input_profile)
    p <- ggpairs(enrichment_result_clipped,
        title = "Scatterplot of all genessets; FDR<0.05 in red", 
        lower = list(continuous = ggpairs_points_plot))
    print(p + theme_bw())
}


ggpairs_contour <- function(res) {
    palette <- colorRampPalette(c("white", "yellow", "orange", "red",
        "darkred", "black"))
    ggpairs_func <- function(data, mapping, ...) {
        p <- ggplot(data = data, mapping = mapping) +
            stat_density2d(aes(fill = ..density..), 
            geom = "tile", contour = FALSE) +
            geom_vline(xintercept = 0, linetype = "dashed") + 
            geom_hline(yintercept = 0, linetype = "dashed") +
            scale_fill_gradientn(colours = palette(25))
        p
    }
    ss <- res$ranked_profile
    p <- ggpairs(as.data.frame(ss),
        title = "Contour plot of all genes after ranking", 
        lower = list(continuous = ggpairs_func),
            diag = list(continuous = wrap("barDiag", 
            binwidth = nrow(ss)/100)))
    print(p + theme_bw())
}


colname_substitute <- function(res) {
    d <- ncol(res$ranked_profile)
    if (d > 5) {
        mydims <- data.frame(attributes(res)$profile_dimensions)
        colnames(mydims) <- "dimensions"
        colnames(res$input_profile) <- paste("d",
            seq_len(ncol(res$input_profile)), sep = "")
        colnames(res$ranked_profile) <- paste("d",
            seq_len(ncol(res$ranked_profile)), sep = "")
    }
    res
}


gene_sector_table <- function(res) {
    mytheme <- gridExtra::ttheme_default(
        core = list(fg_params = list(cex = 0.5)), 
        colhead = list(fg_params = list(cex = 0.7)),
        rowhead = list(fg_params = list(cex = 0.7)))
    
    d <- ncol(res$ranked_profile)
    ss <- res$ranked_profile
    sig <- sign(ss)
    if (d < 6) {
        sig <- sign(ss)
        sector_count <- aggregate(seq(from = 1, to = nrow(sig)) ~ ., 
            sig, FUN = length)
        colnames(sector_count)[ncol(sector_count)] <- 
            "Number of genes in each sector"
        grid.newpage()
        grid.table(sector_count, theme = mytheme)
    }
}


geneset_sector_table <- function(res) {
    d <- ncol(res$ranked_profile)
    mytheme <- gridExtra::ttheme_default(
        core = list(fg_params = list(cex = 0.5)), 
        colhead = list(fg_params = list(cex = 0.7)),
        rowhead = list(fg_params = list(cex = 0.7)))
    
    sig <- 
        sign(res$enrichment_result[which(res$enrichment_result$p.adjustMANOVA< 
        0.05), 4:(4 + d - 1)])
    if (d < 6) {
        if (nrow(sig) > 0) {
            sector_count <- aggregate(seq(from = 1, to = nrow(sig)) ~ ., 
                sig, FUN = length)
            colnames(sector_count)[ncol(sector_count)] <-
                "Number of gene sets in each sector"
            grid.newpage()
            grid.table(sector_count, theme = mytheme)
        }
    }
}

heatmapx <- function(res) {
    d <- ncol(res$ranked_profile)
    resrows = length(res$detailed_sets)
    hmapx <- head(res$enrichment_result[, 4:(4 + d - 1)], resrows)
    rownames(hmapx) <- head(res$enrichment_result$set, resrows)
    colnames(hmapx) <- gsub("^s.", "", colnames(hmapx))
    my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
    heatmap.2(as.matrix(hmapx), scale = "none", margins = c(10, 25),
        cexRow = 0.8, trace = "none", cexCol = 0.8, col = my_palette)
    pl <- recordPlot()
    pl
}

plot3d_detailed_density <- function(res, i) {
    palette <- colorRampPalette(c("white", "yellow", "orange", "red",
        "darkred", "black"))
    d <- ncol(res$ranked_profile)
    ss <- res$ranked_profile
    ll <- res$enrichment_result[i, ]
    size <- ll$setSize
    sss <- res$detailed_sets[[i]]
    
    if (d > 5) {
        colnames(sss) <- paste("d", seq_len(ncol(res$input_profile)), sep = "")
    }
    
    ggpairs_contour_limit_range <- function(data, mapping, ...) {
        p <- ggplot(data = data, mapping = mapping) +
            stat_density2d(aes(fill = ..density..), 
            geom = "tile", contour = FALSE) +
            geom_vline(xintercept = 0, linetype = "dashed") + 
            geom_hline(yintercept = 0, linetype = "dashed") +
            scale_fill_gradientn(colours = palette(25)) + 
            scale_x_continuous(limits = 
                range(min(ss[, gsub("~", "", as.character(mapping[1]))]), 
                max(ss[, gsub("~", "", 
                as.character(mapping[1]))]))) +
            scale_y_continuous(limits = range(min(ss[, 
            gsub("~", "", as.character(mapping[2]))]),
            max(ss[, gsub("~", "", as.character(mapping[2]))])))
        p
    }
    
    p <- ggpairs(as.data.frame(sss), title = ll[, 1],
        lower = list(continuous = ggpairs_contour_limit_range), 
        diag = list(continuous = wrap("barDiag", binwidth = nrow(ss)/10)))
    print(p + theme_bw())
}


plot3d_detailed_points <- function(res, i) {
    d <- ncol(res$ranked_profile)
    ss <- res$ranked_profile
    ll <- res$enrichment_result[i, ]
    size <- ll$setSize
    sss <- res$detailed_sets[[i]]
    
    if (d > 5) {
        colnames(sss) <- paste("d", seq_len(ncol(res$input_profile)), sep = "")
    }
    
    ggpairs_points_limit_range <- function(data, mapping, ...) {
        p <- ggplot(data = data, mapping = mapping) +
            geom_point(alpha = 0.1) +
            geom_vline(xintercept = 0, linetype = "dashed") +
            geom_hline(yintercept = 0, linetype = "dashed") + 
            scale_x_continuous(limits =
                range(min(ss[, gsub("~", "", as.character(mapping[1]))]), 
                max(ss[, gsub("~", "", as.character(mapping[1]))]))) +
            scale_y_continuous(limits = range(min(ss[, 
            gsub("~", "", as.character(mapping[2]))]), 
            max(ss[, gsub("~", "", as.character(mapping[2]))])))
        p
    }
    
    p <- ggpairs(as.data.frame(sss), title = ll[, 1],
        lower = list(continuous = ggpairs_points_limit_range), 
        diag = list(continuous = wrap("barDiag", binwidth = nrow(ss)/10)))
    print(p + theme_bw())
}


plot3d_detailed_violin <- function(res, i) {
    d <- ncol(res$ranked_profile)
    ss <- res$ranked_profile
    ll <- res$enrichment_result[i, ]
    sss <- res$detailed_sets[[i]]
    
    if (d > 5) {
        colnames(sss) <- paste("d", seq_len(ncol(res$input_profile)), sep = "")
    }
    ss_long <- melt(ss)
    sss_long <- melt(sss)
    p <- ggplot(ss_long, aes(Var2, value)) +
        geom_violin(data = ss_long, fill = "grey", colour = "grey") +
        geom_boxplot(data = ss_long, width = 0.9, fill = "grey", 
        outlier.shape = NA, coef = 0) +
        geom_violin(data = sss_long, fill = "black", colour = "black") +
        geom_boxplot(data = sss_long, width = 0.1, outlier.shape = NA) + 
        labs(y = "Position in rank", title = ll[, 1])
    
    print(p + theme_bw() +
        theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        plot.title = element_text(size = 20)))
}


#' mitch_plots
#'
#' This function generates several plots of multivariate gene set enrichment in
#' high resolution PDF format.
#' The number of detailed sets to generate is dictated by the resrows set in
#' the mitch_calc command.
#' @param res a mitch results object.
#' @param outfile the destination file for the plots in PDF format. should
#' contain 'pdf' suffix. Defaults to 
#' 'Rplots.pdf'
#' @return generates a PDF file containing enrichment plots.
#' @keywords mitch plot plots pdf 
#' @export
#' @examples
#' data(resExample)
#' mitch_plots(resExample,outfile='outres.pdf')
#' @import grDevices
#' @import graphics
#' @import GGally
#' @import grid
#' @import gridExtra
#' @importFrom beeswarm beeswarm
#' @importFrom gplots heatmap.2
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom MASS kde2d
mitch_plots <- function(res, outfile = "Rplots.pdf") {
    
    resrows = length(res$detailed_sets)
    d = ncol(res$ranked_profile)
    pdf(outfile)
    
    if ( d>20 ) {
        stop("Error: mitch plotting features are impractical for over 20
        dimensions.")
    }

    if (d == 1) {
        
        plot1d_profile_dist(res)
        plot_geneset_hist(res)
        plot1d_volcano(res)
        lapply(seq_len(resrows), function(i) {
            plot1d_detailed(res, i)
        })
        
    } else if (d == 2) {
        
        plot2d_profile_dist(res)
        plot2d_profile_density(res)
        plot2d_gene_quadrant_barchart(res)
        plot_geneset_hist(res)
        plot2d_set_quadrant_barchart(res)
        plot2d_set_scatter(res)
        plot2d_set_scatter_top(res)
        plot2d_heatmap(res)
        plot_effect_vs_significance(res)
        
        lapply(seq_len(resrows), function(i) {
            plot2d_detailed_density(res, i)
            plot2d_detailed_scatter(res, i)
            plot2d_detailed_violin(res, i)
        })
    } else if (d > 2) {
        
        res <- colname_substitute(res)
        ggpairs_points(res)
        ggpairs_contour(res)
        gene_sector_table(res)
        plot_geneset_hist(res)
        geneset_sector_table(res)
        ggpairs_points_subset(res)
        heatmapx(res)
        plot_effect_vs_significance(res)
        
        lapply(seq_len(resrows), function(i) {
            plot3d_detailed_density(res, i)
            plot3d_detailed_points(res, i)
            plot3d_detailed_violin(res, i)
        })
    }
    dev.off()
}


#' mitch_report
#'
#' This function generates an R markdown based html report containing tables
#' and several plots of mitch results 
#' The plots are in png format, so are not as high in resolution as compared to
#' the PDF generated by mitch_plots 
#' function. The number of detailed sets to generate is dictated by the resrows
#' set in the mitch_calc command.
#' @param res a mitch results object.
#' @param outfile the destination file for the html report. should contain
#' 'html' suffix. Defaults to 
#' 'Rplots.pdf'
#' @param overwrite should overwrite the report file if it already exists?
#' @return generates a HTML file containing enrichment plots.
#' @keywords mitch report html markdown knitr
#' @export
#' @examples
#' data(resExample)
#' mitch_report(resExample,'outres2.html')
#' @import knitr
#' @importFrom rmarkdown render
#' @import echarts4r
mitch_report <- function(res, outfile , overwrite=FALSE) {

    df <- data.frame(dummy_x = seq(20), dummy_y = rnorm(20, 10, 3))
    trash<-df %>% 
        e_charts(dummy_x) %>% 
        e_scatter(dummy_y, symbol_size = 10)

    DIRNAME <- normalizePath(dirname(outfile))
    HTMLNAME <- paste( basename(outfile), ".html", sep = "")
    HTMLNAME <- gsub(".html.html$", ".html", HTMLNAME)
    
    if (file.exists(HTMLNAME)) {
        if (overwrite==FALSE) {
            stop("Error: the output HTML file aready exists.")
        } else {
            message("Note: overwriting existing report")
        }
    }

    if (!file.exists(DIRNAME)) {
        stop("Error: the output folder does not exist.")
    }
    
    rmd_tmpdir <- tempdir()
    rmd_tmpfile <- paste(rmd_tmpdir, "/mitch.Rmd", sep = "")
    html_tmp <- paste(paste(rmd_tmpdir, "/mitch_report.html", sep = ""))
    
    DATANAME <- gsub(".html$", ".RData", HTMLNAME)
    DATANAME <- paste(rmd_tmpdir, "/", DATANAME, sep = "")
    save.image(DATANAME)
    MYMESSAGE = paste("Dataset saved as \"", DATANAME, "\".")
    message(MYMESSAGE)
    
    knitrenv <- new.env()
    assign("DATANAME", DATANAME, knitrenv)
    assign("res", res, knitrenv)
    
    rmd = system.file("mitch.Rmd", package = "mitch")
    rmarkdown::render(rmd, intermediates_dir = "." , output_file = html_tmp)
    file.copy(html_tmp, outfile, overwrite=overwrite)
}



#' Reactome gene sets
#'
#' Genesets from Reactome database suitable for enrichment analysis.
#' Acquired August 2019. The structure of this data is a named list of
#' vectors, containing human gene names (character strings). This is 
#' a sample of 200 gene sets from the approximately 2000 present in the
#' full dataset.
#' @docType data
#' @usage data(genesetsExample)
#' @format A list of gene sets
#' @keywords datasets
#' @references Fabregat et al. (2017) BMC Bioinformatics volume 18, 
#' Article number: 142, https://www.ncbi.nlm.nih.gov/pubmed/28249561
#' @source Reactome website: https://reactome.org/
#' @examples
#' data(genesetsExample)
"genesetsExample"

#' H3K36ac profile
#' 
#' Example edgeR result of differential ChIP-seq H3K36ac.
#' This is a dataframe which contains columns for log fold change, log counts
#' per million, p-value and FDR adjusted p-value. These columns consist of 
#' numerical values. The row names represent human gene names. This is a sample
#' of 1000 gene of an original dataset that contains measurements of ~30000 
#' genes.
#' @docType data
#' @usage data(k36a)
#' @format data frame
#' @keywords datasets
#' @examples
#' data(k36a)
"k36a"

#' H3K9ac profile
#' 
#' Example edgeR result of differential ChIP-seq H3K9ac.
#' This is a dataframe which contains columns for log fold change, log counts
#' per million, p-value and FDR adjusted p-value. These columns consist of 
#' numerical values. The row names represent human gene names. This is a sample
#' of 1000 gene of an original dataset that contains measurements of ~30000 
#' genes.
#' @docType data
#' @usage data(k9a)
#' @format data frame     
#' @keywords datasets
#' @examples
#' data(k9a)
"k9a"

#' RNA profile
#' 
#' Example edgeR result of differential RNA expression.
#' This is a dataframe which contains columns for log fold change, log counts
#' per million, p-value and FDR adjusted p-value. These columns consist of 
#' numerical values. The row names represent human gene names. This is a sample
#' of 1000 gene of an original dataset that contains measurements of ~15000 
#' genes.
#' @docType data
#' @usage data(rna)
#' @format data frame     
#' @keywords datasets
#' @examples
#' data(rna)
"rna"

#' myList: A list of three edgeR results 
#' 
#' Example edgeR results of differential RNA, H3K9ac and H3K36ac profiling.
#' The structure of this data is a list of three dataframes. Each data frame
#' is 1000 lines only.
#' @docType data
#' @usage data(myList)
#' @format data frame     
#' @keywords datasets
#' @examples
#' data(myList)
"myList"

#' myImportedData: Example imported profiles
#' 
#' Example of three edgeR profiles imported using mitch.
#' The structure of this data is a dataframe where each column represents
#' one of the following profiling datasets after scoring: RNA, H3K9ac and
#' H3K36ac. Each row represents one gene and this dataset contains just 1000
#' rows to keep the example dataset small.
#' @docType data
#' @usage data(myImportedData)
#' @format data frame     
#' @keywords datasets
#' @examples
#' data(myImportedData)
"myImportedData"

#' resExample: Example mitch result
#' 
#' Example of mitch results. Enrichment of the Reactome gene sets in the RNA,
#' H3K9ac and H3K36ac datasets. The structure of this data set is a list where
#' the 1st element is "input_profile" that has been imported (data frame), 2nd
#' element is the "input_genesets" (names list of gene names [character 
#' vectors]), 3rd is "ranked_profile" which is the input profiling data after
#' ranking (data frame), 4th is "enrichment_result" which is a data frame which
#' provides enrichment information on each gene set in the profiling data
#' including s scores, p-values and FDR adjusted p-values. 5th is 
#' "analysis_metrics" (list). The 6th slot is "detailed_sets" which is a list
#' of 5 matrices which details the enrichment of members of selected gene sets.
#' @docType data
#' @usage data(resExample)
#' @format list of mixed data types
#' @keywords datasets
#' @examples
#' data(resExample)
"resExample"
