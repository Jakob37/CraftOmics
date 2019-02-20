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
            
            raw_df <- cbind(truth=truth_vector, as.data.frame(pvals_list))
            
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

                warning("\nKnown to have been trouble in the AUC labelling before leading to shuffling of labels, double check\n")
                
                plt <- plt + 
                    scale_color_brewer(palette="Dark2") +
                    scale_color_hue(labels=paste0(sort(names(pvals_list)), " (", round(auc, round), ")"))
            }
            else {
                plt <- plt
            }
                
            plt
        },
        
        auc = function(plt_obj) {
            plotROC::calc_auc(plt_obj)
        }
    ),
    private = list()
)

roc <- ROCVisuals$new()
print("Loading module to 'roc'")
