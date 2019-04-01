library(R6)
library(plotROC)
library(ggplot2)
library(tidyverse)

ROCVisuals <- R6Class(
    public = list(

        roc = function(truth_vector, pvals_list, annot=NULL, reversed_size_importance=TRUE, 
                       cuts=0, round=3, ymin=0, ymax=1, legend_name="Stat type", title="ROC", 
                       show_auc=FALSE, subset=NULL) {
            
            if (!all(unname(unlist(lapply(pvals_list, length))) == length(truth_vector))) {
                stop(
                    "truth_vector and all pvals_lists must be same length\n",
                    "Found truth_vector: ", length(truth_vector), " and pvals_list: ", 
                    paste(unname(unlist(lapply(scores, length))), collapse=", ")
                )
            }
            
            if (!is.null(subset)) {
                warning("Subsetting, only for debugging purposes")
                truth_vector <- head(truth_vector, n=subset)
                pvals_list <- lapply(pvals_list, head, n=subset)
            }
            
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
        },
        
        perf_measures = function(confusion_matrix, true_pattern, false_pattern) {
            warning("DEPRECATED use the one in the MachineLearningTools module instead")
            
            tn <- confusion_matrix$table[false_pattern, false_pattern]
            tp <- confusion_matrix$table[true_pattern, true_pattern]
            fn <- confusion_matrix$table[false_pattern, true_pattern]
            fp <- confusion_matrix$table[true_pattern, false_pattern]
            
            sensitivity <- tp / (tp + fn)
            specificity <- tn / (tn + fp)
            pos_pred <- tp / (tp + fp)
            neg_pred <- tn / (tn + fn)
            accuracy <- (tp + tn) / (tp + fp + fn + tn)
            mcc <- (tp * tn - fp * fn) / (sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
            
            data.frame(
                sensitivity,
                specificity,
                pos_pred,
                neg_pred,
                accuracy,
                mcc,
                tp,
                tn,
                fp,
                fn,
                tot=tp+tn+fp+fn
            )
        },
        
        show_measures = function(measure_df, ignores=NULL) {
            
            warning("DEPRECATED use the one in MachineLearningTools module instead")
            
            measure_df$analysis <- rownames(measure_df)
            long_df <- tidyr::gather(measure_df, "measure", "value", -analysis, -one_of(ignores))
            ggplot(long_df) +
                geom_bar(aes(x=analysis, y=value, fill=measure), stat="identity", position=position_dodge()) +
                ggtitle("Stat. measures")
        }
    ),
    private = list()
)

roc <- ROCVisuals$new()
print("Loading module to 'roc'")
