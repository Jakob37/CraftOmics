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
    )
)

sigtabvis <- SignificanceTableVis$new()
stv <- sigtabvis
print("Loading module to 'sigtabvis' and 'stv'")
