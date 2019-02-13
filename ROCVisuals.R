library(R6)
library(plotROC)
library(ggplot2)
library(tidyverse)

ROCVisuals <- R6Class(
    public = list(

        roc = function(truth_vector, pvals_list, annot=NULL, reversed_size_importance=T, 
                       cuts=0, round=3, ymin=0, ymax=1, legend_name="Stat type", title="ROC", show_auc=F) {
            
            if (typeof(pvals_list) != "list") {
                warning("pvals_list argument should be provided as a list of vectors.",
                        "Now found: ", typeof(pvals_list), " will probably misbehave")
            }
            
            raw_df <- cbind(truth=truth_vector, 
                        as.data.frame(pvals_list))
            
            if (!is.null(annot)) {
                raw_df <- cbind(raw_df, annot=annot)
            }
            
            shaped_df <- raw_df %>% gather(stat_type, val, names(pvals_list))
            
            if (reversed_size_importance) {
                shaped_df$val_importance <- 1 - shaped_df$val
            }
            else {
                shaped_df$val_importance <- shaped_df$val
            }
            
            plt <- ggplot(shaped_df, aes(d=truth, m=val_importance, color=stat_type)) + 
                geom_roc(n.cuts=cuts, labelround=round) + 
                ylim(ymin, ymax) + 
                guides(color=guide_legend(title=legend_name)) +
                ggtitle(title) +
                xlab("FPR") +
                ylab("TPR")
            
            if (show_auc) {
                auc <- self$auc(plt)$AUC                

                plt <- plt + 
                    scale_color_brewer(palette="Dark2") +
                    scale_color_hue(labels=paste0(names(pvals_list), " (", round(auc, round), ")"))
            }
            else {
                plt <- plt
            }
                
            plt
        },
        
        auc = function(plt_obj) {
            plotROC::calc_auc(plt_obj)
        },
        
        # Input format:
        # -> List containing data frames
        # -> Each data frame contains Normalyzer dataframe
        # Out of the data frame, true positives (spike-ins) and background are selected based on patterns in one of the annotation columns
        # Significance levels are then checked from selected p- or q-value column
        generate_roc_plot = function(entry_dfs, annot_col, spikein_pattern, sig_col,
                                      use_fdr=T, only_normal=F, debug=F, sig_coord_thres=0.1, y_cutoff=0, title="ROC curve") {

            warning("This code is deprecated. Use the 'roc' instead")
            
            comb_df <- data.frame(pvals=c(), type=c(), method=c())
            sig_coord_df <- data.frame(x=c(), y=c(), method=c())
            plot_ylabs <- "y-label"
            aucs <- c()

            for (method_i in 1:length(entry_dfs)) {

                name <- names(entry_dfs)[method_i]
                target_df <- entry_dfs[[method_i]]

                df_spikein <- target_df[which(grepl(spikein_pattern, target_df[, annot_col])),]
                df_background <- target_df[which(!grepl(spikein_pattern, target_df[, annot_col])),]

                ci <- CurveInfo$new(target_sig=df_spikein[, sig_col], background_sig=df_background[, sig_col])
                sig_coord <- NULL

                sig_numbers <- cbind(pvals=ci$target_sig, type=1, sample=method_i)
                background_numbers <- cbind(pvals=ci$background_sig, type=0, sample=method_i)

                unsorted_df <- rbind(sig_numbers, background_numbers)
                rownames(unsorted_df) <- seq(1, nrow(unsorted_df))
                sorted_df <- unsorted_df[order(unsorted_df[, "pvals"]),]

                for (i in 1:nrow(sorted_df)) {

                    row <- sorted_df[i,]
                    if (row["type"] == 0) {
                        ci$new_neg_found()
                    }
                    else {
                        ci$new_pos_found()
                    }

                    if (is.null(sig_coord) && row["pvals"] >= sig_coord_thres) {
                        sig_coord <- data.frame(x=ci$last_x_coord, y=ci$last_y_coord, method=name)
                        sig_coord_df <- rbind(sig_coord_df, sig_coord)
                    }
                }

                coord_df <- data.frame(FPR=ci$x_vals, TPR=ci$y_vals, sample=name)
                comb_df <- rbind(comb_df, coord_df)

                aucs <- c(aucs, self$calculate_auc_(ci$x_vals, ci$y_vals))
            }

            orig_lvs <- levels(comb_df$sample)
            new_lvs <- paste(round(aucs, 3), orig_lvs, sep=", ")
            levels(comb_df$sample) <- new_lvs
            sig_coord_df$method <- new_lvs

            plt <- self$plot_curve_(comb_df, sig_coord_df, title=title, y_cutoff=y_cutoff)

            return(plt)
        },

        calculate_auc_ = function(x_coord, y_coord) {
            
            warning("This code is deprecated")

            tot_sum <- 0
            prev_x <- x_coord[1]
            for (i in 2:length(x_coord)) {

                x <- x_coord[i]
                y <- y_coord[i]
                x_diff <- x - prev_x

                if (x_diff != 0) {
                    diff_area <- y * x_diff
                    tot_sum <- tot_sum + diff_area
                }
                prev_x <- x
            }

            tot_sum
        },

        plot_curve_ = function(comb_df, sig_coord_df, title, y_cutoff, debug=F) {
            
            warning("This code is deprecated")

            colorCount <- length(unique(comb_df$sample))
            getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

            plt <- ggplot() +
                geom_line(data=comb_df, aes_string("FPR", "TPR", color="sample")) +
                ggtitle(title) +
                geom_point(data=sig_coord_df, aes_string(x="x", y="y", color="method")) +
                ylim(y_cutoff, NA) +
                scale_color_manual(values=getPalette(colorCount))

            return(plt)
        }
    ),
    private = list()
)

CurveInfo = R6Class("CurveInfo",

                    public = list(

                        target_sig = NULL,
                        background_sig = NULL,

                        found_pos=0,
                        found_neg=0,
                        last_x_coord=0,
                        last_y_coord=0,
                        x_vals=c(0),
                        y_vals=c(0),

                        initialize = function(target_sig, background_sig) {

                            warning("This code is deprecated")
                            self$target_sig <- target_sig
                            self$background_sig <- background_sig
                        },
                        new_pos_found = function() {
                            self$found_pos <- self$found_pos + 1
                            self$last_y_coord <- self$found_pos / self$get_tot_true()
                            private$update_coords()
                        },
                        new_neg_found = function() {
                            self$found_neg <- self$found_neg + 1
                            self$last_x_coord <- self$found_neg / self$get_tot_false()
                            private$update_coords()
                        },
                        get_tot_true = function() {
                            length(self$target_sig)
                        },
                        get_tot_false = function() {
                            length(self$background_sig)
                        }
                    ),
                    private = list(
                        update_coords = function() {
                            self$x_vals <- c(self$x_vals, self$last_x_coord)
                            self$y_vals <- c(self$y_vals, self$last_y_coord)
                        }
                    )
)

roc <- ROCVisuals$new()
print("Loading module to 'roc'")
