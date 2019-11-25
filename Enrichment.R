library(R6)
library(tidyverse)

# Extra information:
# For GAGE, the type of input IDs is important
# KEGG sets can be setup using gage::kegg.gsets(species='ath')
# GO sets can be setup using gage::go.gsets(species='Arabidopsis')

Enrichment <- R6Class(
    public = list(
        
        perform_gage_gse = function(mat, annot_rows, cond_cols, levels, annot_sets, quiet=FALSE, same_dir=TRUE) {
        
            library(gage)
            
            if (!is.matrix(mat)) {
                stop("Input data must be a matrix, found: ", typeof(mat))
            }

            ref.idx <- which(cond_cols == levels[1])
            cond.idx <- which(cond_cols == levels[2])

            message("Using indices: ", paste(ref.idx, collapse=" "), " and ", paste(cond.idx, collapse=" "))

            if (length(ref.idx) == 0 || length(cond.idx) == 0) {
                stop("At least one of provided condition levels had no matches in condition column")
            }

            rownames(mat) <- annot_rows
            
            if (!quiet) message("Calculating annotation sets enrichment...")
            gage_out <- gage::gage(
                mat, 
                gsets=annot_sets, 
                ref=ref.idx, 
                samp=cond.idx, 
                compare="unpaired",
                same.dir=same_dir
            )
            if (!quiet) message("Done!")
            gage_out
        },
        gse_summary = function(gage_mat, pthres=0.05, qthres=0.1, title=NULL) {
            psig_count <- gage_mat %>% data.frame() %>% filter(p.val < pthres) %>% nrow()
            qsig_count <- gage_mat %>% data.frame() %>% filter(q.val < qthres) %>% nrow()
            tot_count <- nrow(gage_mat)
            perc_psig <- round(100 * psig_count / tot_count, 3)
            message(
                "p < ", pthres, ": ", psig_count,
                ", perc psig: ", perc_psig,
                "%, q < ", qthres, ": ", qsig_count,
                ", tot count: ", tot_count
            ) 
        },
        convert_to_entrez = function(db, ids, id_type) {
        
            # Db could be org.At.tair.db
            # select(org.At.tair.db, keys=ids, columns=c(), keytype=c())
            select(org.At.tair.db, columns='ENTREZID', keytype='TAIR')
        },
        stat_hists = function(gsea_mat, title="No title") {
       
            plts <- list()

            ps <- gsea_mat %>% data.frame() %>% dplyr::select(c(p=p.val, q=q.val))
            plts[["p"]] <- ggplot(ps, aes(p)) + 
                geom_histogram(bins=100, na.rm=TRUE) + 
                ggtitle(paste0(title, " (p)")) + 
                xlim(0, 1)

            plts[["q"]] <- ggplot(ps, aes(q)) + 
                geom_histogram(bins=100, na.rm=TRUE) + 
                ggtitle(paste0(title, " (q)")) + 
                xlim(0, 1)
            plts
        },
        heatmap = function() {},
        network_graph = function() {}
    ),
    private = list()
)

enrich <- Enrichment$new()
er <- enrich

message("Loaded Enrichment module as 'enrich' and 'er'")
