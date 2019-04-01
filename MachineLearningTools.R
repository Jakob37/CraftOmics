library(R6)
library(caret)

MachineLearningTools <- R6Class(
    public = list(
        caret_roc_curves = function(trained_models, truth_pattern, roc_module) {
            truth_and_scores <- self$get_roc_lists(trained_models, truth_pattern)
            roc_plt <- roc_module$roc(
                truth_vector=truth_and_scores$truth[[1]],
                pvals_list=truth_and_scores$scores,
                reversed_size_importance=TRUE,
                show_auc=TRUE
            ) + theme_classic()
            roc_plt
        },
        
        caret_performance_df = function(trained_models, high_pattern, low_pattern, roc_aucs=NULL) {
            measures_list <- lapply(trained_models, function(model) {
                conf_mat <- caret::confusionMatrix(model)
                self$perf_measures(conf_mat, high_pattern, low_pattern)
            })
            measures_mat <- do.call("rbind", measures_list)
            sorted_measures_mat <- measures_mat[sort(rownames(measures_mat)), ]
            if (!is.null(roc_aucs)) {
                sorted_measures_mat <- cbind(sorted_measures_mat, auc=roc_aucs$AUC)
            }
            sorted_measures_mat
        },
        
        show_measures = function(measure_df, ignores=c("fp", "fn", "tn", "tp", "tot"), y_to_one=TRUE) {
            
            measure_df$analysis <- rownames(measure_df)
            long_df <- tidyr::gather(measure_df, "measure", "value", -analysis, -one_of(ignores))
            plt <- ggplot(long_df) +
                geom_bar(aes(x=analysis, y=value, fill=measure), stat="identity", position=position_dodge()) +
                ggtitle("Stat. measures") + theme_classic()
            if (y_to_one) {
                plt <- plt + ylim(c(NA, 1))
            }
            plt
        },
        
        perf_measures = function(caret_confusion_matrix, true_pattern, false_pattern) {
            
            tn <- caret_confusion_matrix$table[false_pattern, false_pattern]
            tp <- caret_confusion_matrix$table[true_pattern, true_pattern]
            fn <- caret_confusion_matrix$table[false_pattern, true_pattern]
            fp <- caret_confusion_matrix$table[true_pattern, false_pattern]
            
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
        
        get_variable_importances = function(trained_models) {
            varimp_plts <- list()
            for (model_name in names(trained_models)) {
                varimp_plts[[model_name]] <- ggplot(caret::varImp(trained_models[[model_name]])) +
                    theme_classic() +
                    ggtitle(model_name)
            }
            varimp_plts
        },
        # trained_models is assumed to be a named list of trained Caret models
        # Requires the Caret runs to have been saved by specifying option 'savePredictions'
        get_roc_lists = function(trained_models, truth_pattern) {
            
            all_truths <- list()
            all_scores <- list()
            
            for (model_name in names(trained_models)) {
                model <- trained_models[[model_name]]
                
                tune_names <- names(model$bestTune)
                best_tune_values <- model$bestTune
                all_predictions <- model$pred
                
                # Slicing out values corresponding to best predictive run
                best_predictions_contrast <- apply(
                    all_predictions[tune_names], 
                    1, 
                    function(row, best_tune) { all(row == best_tune_values ) },
                    best_tune=best_tune_values
                )
                
                all_scores[[model_name]] <- model$pred[best_predictions_contrast, ][[truth_pattern]]
                all_truths[[model_name]] <- model$pred[best_predictions_contrast, ]$obs
            }
            list(truth=all_truths, scores=all_scores)
        }
        
    ),
    private = list(
    )
)

mltools <- MachineLearningTools$new()
mlt <- mltools
print("Loading module to 'mltools' and 'mlt")
