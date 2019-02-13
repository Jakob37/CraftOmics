library(R6)
library(grid)
library(gridExtra)
library(tidyverse)
library(Biostrings)

SimulateBatchUtils <- R6Class(
    public = list(
        
        roc_obj = NULL,
        ev_obj = NULL,
        stv_obj = NULL,
        
        initialize = function(roc_obj = NULL, ev_obj=NULL, stv_obj=NULL) {
            self$roc_obj <- roc_obj
            self$ev_obj <- ev_obj
            self$stv_obj <- stv_obj
        },
        
        simulate_batches = function(base_fasta, spikein_fasta, out_dir, base_count, spike_count, 
                        base_intensity_low, base_intensity_high, feature_spread_frac, batch_diff, batch_spread, sample_names, 
                        batch_impacts=c(0,1,0,1,0,1), spike_levels=c(1,1,1,1,1,1), design_fp=NULL, raw_data_fp=NULL, eval_plots=TRUE) {
    
            if (length(sample_names) != length(batch_impacts) || length(sample_names) != length(spike_levels)) {
                stop(
                    "batch_impacts needs to be same length as sample_names and spike_levels. Found: ",
                    "\nsample_names: ", paste(sample_names, collapse=", "),
                    "\nbatch_impacts: ", paste(batch_impacts, collapse=", "),
                    "\nspike_levels: ", paste(spike_levels, collapse=", ")
                )
            }
            
            design_df <- data.frame(sample=sample_names, group=spike_levels, batch=batch_impacts)
            if (!is.null(design_fp)) {
                write_tsv(design_df, design_fp)
                message("Writing design to ", design_fp)
            }
            else {
                message("No design_fp option provided, no design written")
            }
            
            base_fa <- Biostrings::readAAStringSet(base_fasta)
            spikein_fa <- Biostrings::readAAStringSet(spikein_fasta)
            
            base_select_fa <- sample(base_fa, base_count)
            spike_select_fa <- sample(spikein_fa, spike_count)
            
            name_base_base <- vapply(names(base_select_fa), function(name) { unlist(strsplit(name, " "))[1] } , "")
            name_base_spike <- vapply(names(spike_select_fa), function(name) { unlist(strsplit(name, " "))[1] } , "")
            
            # Background features
            background_df <- data.frame(matrix(nrow=0, ncol=length(sample_names)))
            for (feature_i in seq_len(base_count)) {
                
                feature_base <- runif(1, base_intensity_low, base_intensity_high)
                sample_ints <- rnorm(length(sample_names), feature_base, feature_base * feature_spread_frac)
                background_df <- rbind(background_df, sample_ints)
            }
            colnames(background_df) <- sample_names
            
            # Spike-in features
            spike_df <- data.frame(matrix(nrow=0, ncol=length(sample_names)))
            for (feature_i in seq_len(spike_count)) {
                
                feature_base <- runif(1, base_intensity_low, base_intensity_high)
                spike_ints <- rnorm(length(sample_names), feature_base, feature_base * feature_spread_frac) * spike_levels
                spike_df <- rbind(spike_df, spike_ints)
            }
            colnames(spike_df) <- sample_names
            raw_df <- rbind(background_df, spike_df)
            
            # Apply batch effect
            for (sample_i in seq_len(length(sample_names))) {
                raw_df[, sample_i] <- raw_df[, sample_i] + rnorm(nrow(raw_df), batch_diff, batch_spread) * batch_impacts[sample_i]
            }
            
            # Write final proteins
            all_base <- c(name_base_base, name_base_spike)
            all_fa <- c(base_select_fa, spike_select_fa)
            for (sample_name in sample_names) {
                sample_vals <- raw_df[, sample_name]
                names(all_fa) <- paste0(all_base, " [# intensity=", sample_vals, " #]")
                sample_path <- paste0(out_dir, "/", sample_name, ".fa")
                message("Writing ", length(all_fa), " entries with total intensity ", sum(sample_vals), " to ", sample_path)
                Biostrings::writeXStringSet(all_fa, filepath = sample_path)
            }
            
            raw_df <- cbind(annot=as.character(c(name_base_base, name_base_spike)), raw_df)
            write_tsv(raw_df, path=raw_data_fp)
            message("Writing raw results to ", raw_data_fp)
        },
        
        eval_plots = function(ddf, rdf, annot_col_name, spike_pattern=NULL, 
                              sampleCol="sample", groupCol="group", batchCol="batch", 
                              sig_thres=0.1, 
                              pval_name="P.Value", adjp_name="adj.P.Val", expr_name="AveExpr", log_name="logFC") {

            sdf <- rdf %>% select(as.character(ddf[[sampleCol]]))
            annot_col <- rdf %>% select(annot_col_name) %>% unlist()
            
            rdf$pvals <- rdf[[pval_name]]
            rdf$adjp <- rdf[[adjp_name]]
            rdf$means <- rdf[[expr_name]]
            rdf$fold <- rdf[[log_name]]
            rdf$is_sig <- rdf$adjp < sig_thres
                        
            if (!is.null(spike_pattern)) {
                rdf$is_spike <- grepl(spike_pattern, rdf[[annot_col_name]])
                color_col <- "is_spike"
            }
            else {
                color_col <- "is_sig"
            }

            ma_plt <- private$ma(
                rdf, x_lims=NULL, y_lims=NULL, sig_thres=0.2, 
                title="Expression pattern", sig_col_name="is_sig", 
                fold_col_name="fold", expr_col_name="means", na.rm=FALSE, color_col_name=color_col)
            
            vulc_plt <- private$ma(
                rdf, x_lims=NULL, y_lims=NULL, sig_thres=0.2, 
                title="Expression pattern", sig_col_name="is_sig", 
                fold_col_name="fold", expr_col_name="means", na.rm=FALSE, color_col_name=color_col, is_vulc=TRUE)

            if (!is.null(self$ev_obj)) {
                p_plt <- self$ev_obj$pvalhist(rdf$pvals)
                density_plt_group <- self$ev_obj$sample_dist(sdf, color_col=as.factor(ddf[[groupCol]]))
                density_plt_batch <- self$ev_obj$sample_dist(sdf, color_col=as.factor(ddf[[batchCol]]))
            }
            else {
                warning("ev_obj not available, no P-value histogram or density plots shown")
            }

            pca_plt <- private$pca(sdf, ddf$sample, as.factor(ddf[[groupCol]]), pcs=c(1,2))

            test_plt <- grid.arrange(
                ma_plt, vulc_plt, p_plt, pca_plt,
                ncol=2)

            test_plt <- grid.arrange(
                density_plt_group, density_plt_batch,
                ncol=2)
            
            if (!is.null(spike_pattern)) {
                
                roc_plt <- private$roc(list(pvals=rdf$pvals), grepl(spike_pattern, annot_col))
                test_plt <- grid.arrange(roc_plt, ncol=1)
            }
        }
    ),
    private = list(
        
        ma = function(table, x_lims=NULL, y_lims=NULL, sig_thres=0.1,
                      title="Expression pattern", sig_col_name="is_sig",
                      fold_col_name="fold", expr_col_name="expr", color_col_name="spikein", p_col_name=NULL, na.rm=FALSE, is_vulc=FALSE) {
            
            plot_df <- as.data.frame(table)
            plot_df <- plot_df[with(plot_df, order(pvals)),]
            plot_df$log_p <- -log10(plot_df$pvals)
            
            if (!is_vulc) {
                p <- ggplot(plot_df, aes_string(expr_col_name, fold_col_name, colour=color_col_name))
            }
            else {
                p <- ggplot(plot_df, aes_string(fold_col_name, "log_p", colour=color_col_name))
            }
            
            p <- p + 
                ggtitle(title) + 
                geom_point(size=2, na.rm=na.rm, aes_string(shape=sig_col_name)) +
                scale_color_brewer(palette="Dark2") +
                theme_classic()
            
            
            if (!is.null(x_lims)) {
                p <- p + xlim(x_lims)
            }
            
            if (!is.null(y_lims)) {
                p <- p + ylim(y_lims)
            }
            
            p            
        },
        
        pca = function(data_df, samples, group, contrast_levels, pcs=c(1,2)) {
            
            df_pca <- prcomp(t(data_df[complete.cases(data_df),]), scale=T, center=T)
            df_out <- as.data.frame(df_pca$x)
            percentage <- round(df_pca$sdev ^ 2 / sum(df_pca$sdev ^ 2) * 100, 2)
            percentage <- paste0(colnames(df_out), " (", paste0(as.character(percentage), "%)"))
            df_out$sample <- samples
            
            pc1 <- paste0("PC", pcs[1])
            pc2 <- paste0("PC", pcs[2])
            
            plt <- ggplot(df_out, aes_string(x=pc1, y=pc2, color="group", label="sample")) + 
                geom_text() + 
                scale_color_brewer(palette="Dark2") +
                theme_classic() +
                ggtitle("PCA") +
                xlab(percentage[pcs[1]]) +
                ylab(percentage[pcs[2]])
            
            plt
        },
        
        roc = function(pval_list, truth_vector, title="ROC") {
            
            if (is.null(self$roc_obj)) {
                stop("Need to assign 'roc_obj' in constructor to use this method")
            }
            
            roc_plt <- self$roc_obj$roc(
                truth_vector=truth_vector, 
                pvals_list=pval_list,
                reversed_size_importance=TRUE,
                show_auc=TRUE,
                legend_name="Processing type",
                title=title,
                round=4
            ) + 
                ylim(0.8, 1) +
                theme_classic()
            
            roc_plt
        },
        
        accuracy = function() {}
    )
)


simbatch <- SimulateBatchUtils$new()
sb <- simbatch
print("Loading module to 'simbatch' and 'sb'")
