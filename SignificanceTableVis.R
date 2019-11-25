library(ggplot2)
library(ggdendro)
library(ggfortify)
library(RColorBrewer)
library(R6)

SignificanceTableVis <- R6Class(
    public = list(

        ma = function(table, x_lims=NULL, y_lims=NULL, sig_thres=0.1, title="Expression pattern", sig_col_name="adj.P.Val", fold_col_name="logFC", expr_col_name="AveExpr", na.rm=FALSE) {
            
            sig_col <- table[, sig_col_name] < sig_thres
            plot_df <- data.frame(expr=table[, expr_col_name], fold=table[, fold_col_name], sig_val=table[, sig_col_name], sig=sig_col)
            plot_df <- plot_df[with(plot_df, order(sig)),]

            p <- ggplot(plot_df, aes(expr, fold, colour=sig)) +
                geom_point(size=1, alpha=0.5, na.rm=na.rm) +
                ggtitle(title)

            if (!is.null(x_lims)) {
                p <- p + xlim(x_lims)
            }

            if (!is.null(y_lims)) {
                p <- p + ylim(y_lims)
            }

            p
        },

        vulc = function(table, x_lims=NULL, y_lims=NULL, sig_thres=0.1, title="Expression pattern", sig_col_name="adj.P.Val", fold_col_name="logFC", expr_col_name="AveExpr", na.rm=FALSE) {
            
            sig_col <- table[, sig_col_name] < sig_thres
            plot_df <- data.frame(expr=table[, expr_col_name], fold=table[, fold_col_name], sig_val=table[, sig_col_name], sig=sig_col)
            plot_df <- plot_df[with(plot_df, order(sig)),]

            p <- ggplot(plot_df, aes(fold, -log10(sig_val), colour=sig)) +
                geom_point(size=1, alpha=0.5, na.rm=na.rm) +
                ggtitle(title)

            if (!is.null(x_lims)) {
                p <- p + xlim(x_lims)
            }

            if (!is.null(y_lims)) {
                p <- p + ylim(y_lims)
            }

            p
        },
        
        adf_vulc = function(adf, name, title, show_ids=0, ylim=NULL, xlim=NULL, sig_thres=0.05, fold_cutoff=NULL, p_cutoff=FALSE, id_col=NULL, id_trim=NULL, colors=c("#AAAAAA", "#0C81A2"), ...) {
            
            logfc_col <- sprintf("%s.logFC", name)
            pval_col <- sprintf("%s.P.Value", name)
            fdr_col <- sprintf("%s.adj.P.Val", name)
            
            adf <- adf %>% filter(!is.na(UQ(as.name(pval_col)))) %>% arrange(UQ(as.name(sprintf("%s.P.Value", name))))
            adf$min_log10_pval <- -log10(adf[[pval_col]])
            adf$min_log10_qval <- -log10(adf[[fdr_col]])
            
            if (p_cutoff) {
                stat_col <- pval_col
                legend_label <- "P-val"
            }
            else {
                stat_col <- fdr_col
                legend_label <- "FDR"
            }
            
            if (is.null(fold_cutoff)) {
                adf$is_sig <- adf[[stat_col]] < sig_thres
                fold_label <- ""
            }
            else {
                adf$in_fold_range <- abs(adf[[logfc_col]]) > fold_cutoff
                adf$is_sig <- adf$in_fold_range & adf[[stat_col]] < sig_thres
                fold_label <- sprintf(", log2 fold > %s", fold_cutoff)
            }
            
            plt <- ggplot(adf, aes_string(x=logfc_col, y="min_log10_pval", color="is_sig")) + 
                geom_point(alpha=0.4) + 
                ggtitle(title) + 
                xlab("log2 fold") + 
                ylab("-log10(P-value)") + 
                labs(color=sprintf("%s < %s%s", legend_label, sig_thres, fold_label)) +
                scale_color_manual(values=colors)
            
            if (!is.null(xlim)) {
                plt <- plt + xlim(xlim)
            }
            
            if (!is.null(ylim)) {
                plt <- plt + ylim(ylim)
            }
            
            if (show_ids > 0) {
                adf$labels <- adf[[id_col]]
                if (!is.null(id_trim)) {
                    adf$labels <- adf$labels %>% gsub(sprintf("%s.*", id_trim), "", .)
                }
                adf$labels[seq(show_ids+1, nrow(adf))] <- ""
                plt + geom_text_repel(data=adf, aes(label=labels), min.segment.length = 0, label_size=0.15)
            }
            else {
                plt
            }
        },
        
        threeway_venn = function(sig_table, contrasts, thres_col_base, thres, fold_col_base="logFC",
                                 check_greater_than=FALSE) {
            
            thres_cols <- paste(contrasts, thres_col_base, sep=".")
            fold_cols <- paste(contrasts, fold_col_base, sep=".")
            
            sig_vals <- apply(sig_table[, thres_cols], 2, function(col, thres, inverse) {
                if (!inverse) {
                    col < thres
                }
                else {
                    col > thres
                }
            }, thres=thres, inverse=check_greater_than)

            # 1 if same sign, -1 if contra sign            
            fold_signs <- apply(sig_table[, fold_cols], 2, function(col) {
                2 * as.numeric(col > 0) - 1
            })
            
            colnames(sig_vals) <- paste(c("A", "B", "C"), "sig", sep="_")
            colnames(fold_signs) <- paste(c("A", "B", "C"), "same", sep="_")
            
            venn_stat_df <- private$calculate_stat_table(sig_vals, fold_signs)
            tot_counts <- private$get_tot_counts(venn_stat_df)
            contra_counts <- private$get_contra_counts(venn_stat_df)
            same_counts <- tot_counts - contra_counts
            
            circle_spacing <- 0.6
            df.circles <- data.frame(
                x = c(-circle_spacing, circle_spacing, 0),
                y = c(circle_spacing, circle_spacing, -circle_spacing),
                labels = contrasts
            )
            
            x_number_positions <- c(-1.5, 1.5, 0, 0, 1, -1, 0)
            y_number_positions <- c(0.8, 0.8, -1.6, 1.3, -0.3, -0.3, 0.15) 
            
            df.data <- data.frame(
                x=x_number_positions,
                y=y_number_positions,
                tot_counts=tot_counts,
                contra_counts=contra_counts,
                same_counts=same_counts
            )
            df.data$label <- c(
                tot_counts[1:3],
                paste(df.data$same_counts[4:7], df.data$contra_counts[4:7], sep="\n")
            )
            
            ggplot(data=df.circles) +
                ggforce::geom_circle(
                    aes_string(x0="x", y0="y", r=1.5, fill="labels"),
                    alpha=0.2, size=0.5, color="darkgray"
                ) +
                coord_fixed() +
                theme_void() +
                theme(legend.position="right", legend.direction="vertical") +
                scale_fill_manual(values = c("#009900", "#000099", "#990000")) +
                scale_color_manual(values = c("#009900", "#000099", "#990000"), guide=FALSE) +
                labs(fill=NULL) +
                annotate("text", x=df.data$x, y=df.data$y, label=df.data$label, size=5) +
                ggtitle("") + theme(plot.title = element_text(hjust=0.5))
            
        },
        
        plot_comp_venns = function(sig_table, contrasts, sig_thres=0.1, base_sig_col="adj.P.Val", 
                                   log2_fold_thres=0, base_fold_name="logFC", check_greater_than=FALSE) {
            
            contrast_sig <- list()
            folds <- list()
            for (contrast in contrasts) {
                
                sig_col <- paste(contrast, base_sig_col, sep=".")
                fold_name <- paste(contrast, base_fold_name, sep=".")
                if (!check_greater_than) {
                    sigs <- sig_table[[sig_col]] < sig_thres
                }
                else {
                    sigs <- sig_table[[sig_col]] > sig_thres
                }
                contrast_sig[[contrast]] <- sigs
                folds[[contrast]] <- sig_table[[fold_name]]
            }
            
            comp_count <- length(contrasts)
            plts <- list()
            index <- 0
            
            for (ind_out in seq_len(comp_count-1)) {
                for (ind_in in (ind_out+1):comp_count) {
                    
                    # Credit goes to:
                    # https://scriptsandstatistics.wordpress.com/2018/04/26/how-to-plot-venn-diagrams-using-r-ggplot2-and-ggforce/
                    
                    sig_inds <- as.data.frame(contrast_sig[c(ind_in, ind_out)])
                    
                    both_sig_inds <- which(apply(sig_inds, 1, all))
                    contra_count <- length(which(sign(folds[[ind_in]][both_sig_inds]) != sign(folds[[ind_out]][both_sig_inds])))
                    
                    vdc <- limma::vennCounts(sig_inds)
                    class(vdc) <- 'matrix'
                    df.vdc <- data.frame(vdc[-1, ])
                    df.vdc <- rbind(df.vdc, c(-1, -1, contra_count))
                    df.vdc$x <- c(-1.5, 1.5, 0, 0)
                    df.vdc$y <- c(0, 0, 0, -0.3)
                    df.vdc$label <- df.vdc$Counts
                    df.vdc$label[3] <- paste0(df.vdc$label[3] - df.vdc$label[4], " same")
                    df.vdc$label[4] <- paste0(df.vdc$label[4], " contra")
                    
                    df.venn <- data.frame(x = c(-0.5, 0.5),
                                          y = c(0, 0),
                                          labels = c(contrasts[ind_out], contrasts[ind_in]))
                    
                    if (!check_greater_than) {
                        title <- paste0(base_sig_col, " < ", sig_thres)
                    }
                    else {
                        title <- paste0(base_sig_col, " > ", sig_thres)
                    }
                    
                    if (log2_fold_thres != 0) {
                        title <- paste0(title, ", |log2 fold| >= ", log2_fold_thres)
                    }
                    
                    plt <- ggplot2::ggplot(data=df.venn) +
                        ggforce::geom_circle(
                            ggplot2::aes_string(x0 = "x", y0 = "y", r = 1.5, fill = "labels"), 
                            alpha = .3, size = 0.5, colour = 'darkgray') +
                        ggplot2::coord_fixed() +
                        ggplot2::theme_void() +
                        ggplot2::theme(legend.position = 'bottom', legend.direction='vertical') +
                        ggplot2::scale_fill_manual(values = c('cornflowerblue', 'firebrick',  'gold')) +
                        ggplot2::scale_colour_manual(values = c('cornflowerblue', 'firebrick', 'gold'), guide = FALSE) +
                        ggplot2::labs(fill = NULL) +
                        ggplot2::annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$label, size = 5) +
                        ggplot2::ggtitle(title)
                    
                    index <- index + 1
                    plts[[index]] <- plt
                }
            }
            plts
        }
    ),
    private = list(
        get_tot_counts = function(stat_df) {
            tot_counts <- c(
                stat_df %>% filter(A_sig & !ABC_sig & !AB_sig & !AC_sig) %>% nrow(),
                stat_df %>% filter(B_sig & !ABC_sig & !AB_sig & !BC_sig) %>% nrow(),
                stat_df %>% filter(C_sig & !ABC_sig & !AC_sig & !BC_sig) %>% nrow(),
                stat_df %>% filter(AB_sig & !ABC_sig) %>% nrow(),
                stat_df %>% filter(BC_sig & !ABC_sig) %>% nrow(),
                stat_df %>% filter(AC_sig & !ABC_sig) %>% nrow(),
                stat_df %>% filter(ABC_sig) %>% nrow()
            )
            tot_counts
        },
        
        get_contra_counts = function(stat_df) {
            contra_counts <- c(
                0,
                0,
                0,
                stat_df %>% filter(AB_sig & !ABC_sig) %>% filter(!AB_same) %>% nrow(),
                stat_df %>% filter(BC_sig & !ABC_sig) %>% filter(!BC_same) %>% nrow(),
                stat_df %>% filter(AC_sig & !ABC_sig) %>% filter(!AC_same) %>% nrow(),
                stat_df %>% filter(ABC_sig) %>% filter(!ABC_same) %>% nrow()
            )
            contra_counts
        },
        
        calculate_stat_table = function(sig_vals, fold_signs) {
            venn_stat_df <- data.frame(cbind(sig_vals, fold_signs))
            
            venn_stat_df$AB_sig <- venn_stat_df$A_sig & venn_stat_df$B_sig
            venn_stat_df$AC_sig <- venn_stat_df$A_sig & venn_stat_df$C_sig
            venn_stat_df$BC_sig <- venn_stat_df$B_sig & venn_stat_df$C_sig
            venn_stat_df$ABC_sig <- venn_stat_df$A_sig & venn_stat_df$B_sig & venn_stat_df$C_sig
            
            venn_stat_df$AB_same <- venn_stat_df$A_same == venn_stat_df$B_same
            venn_stat_df$AC_same <- venn_stat_df$A_same == venn_stat_df$C_same
            venn_stat_df$BC_same <- venn_stat_df$B_same == venn_stat_df$C_same
            venn_stat_df$ABC_same <- venn_stat_df$A_same == venn_stat_df$B_same & venn_stat_df$B_same == venn_stat_df$C_same
            venn_stat_df
        }
    )
)

sigtabvis <- SignificanceTableVis$new()
stv <- sigtabvis
print("Loading module to 'sigtabvis' and 'stv'")
