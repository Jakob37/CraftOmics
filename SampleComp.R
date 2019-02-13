library(R6)

SampleComp <- R6Class(

    public = list(

        get_fold = function(table, fold_col="logFC") {
            fold_text <- function(fold) {
                if (is.na(fold)) {
                    NA
                }
                else if (fold > 0) {
                    paste("UP", round(fold, 3))
                }
                else {
                    paste("DOWN", round(fold, 3))
                }
            }

            sorted_table <- table[ order(as.numeric(as.character(rownames(table)))), ]
            folds <- sapply(sorted_table[, fold_col], fold_text)
            names(folds) <- rownames(sorted_table)
            folds
        },

        get_shared_features = function(raw_df, annot_df, annot_key, limma_table1, limma_table2, limma_table3=NULL,
                                       sig_thres=0.1, label1="label1", label2="label2", label3="label3", sig_col="adj.P.Val") {

            var1_rows <- private$get_sig_rows(limma_table1, sig_thres=sig_thres, sig_col=sig_col)
            var2_rows <- private$get_sig_rows(limma_table2, sig_thres=sig_thres, sig_col=sig_col)
            var12_rows <- intersect(var1_rows, var2_rows)

            if (is.null(limma_table3)) {
                folds_df <- cbind(self$get_fold(limma_table1),
                                  self$get_fold(limma_table2))
                colnames(folds_df) <- paste("regulation_", c(label1, label2), sep="")
            }
            else {
                folds_df <- cbind(self$get_fold(limma_results),
                                  self$get_fold(limma_results),
                                  self$get_fold(limma_results))
                colnames(folds_df) <- paste("regulation_", c(label1, label2, label3), sep="")
            }

            if (!is.null(limma_table3)) {
                var3_rows <- private$get_sig_rows(limma_table3, sig_thres=sig_thres, sig_col=sig_col)
                var3_folds <- self$get_fold(limma_table3)
                var123_rows <- intersect(var12_rows, var3_rows)
                target_rows <- var123_rows
            }
            else {
                target_rows <- var12_rows
            }

            if (length(target_rows) == 0) {
                return(NULL)
            }

            target_folds <- folds_df[as.character(target_rows), ]
            target_raw_rows <- raw_df[target_rows, ]
            target_annot_col <- target_raw_rows[, annot_key]

            # target_annot_col_trimmed <- sapply(as.character(target_annot_col), function(annot) { unlist(strsplit(annot, ";"))[1] })
            # target_full_annotation <- annot_df[target_annot_col_trimmed, ]
            target_full_annotation <- annot_df[target_rows,]

            cbind(target_folds, target_full_annotation, sig1=limma_table1[target_rows, sig_col], sig2=limma_table2[target_rows, sig_col])
        },

        generate_peptide_table = function(raw_df, full_de_table, sig_thres, annot_names, included_sample_cols=NULL,
                                          annot_key=NULL, annot_df=NULL, omit_data=F, sig_col="adj.P.Val") {
            # generate_peptide_table = function(raw_df, design_m, full_de_table, sig_thres, annot_names,
            #                                   annot_key=NULL, annot_df=NULL, omit_data=F, sig_col="adj.P.Val") {

            if (is.null(included_sample_cols)) {
                all_cols <- c(annot_names)
            }
            else {
                included_sample_cols <- as.character(included_sample_cols)
                all_cols <- c(annot_names, included_sample_cols)
            }

            print(head(full_de_table))

            sig_table <- full_de_table[which(full_de_table[, sig_col] < sig_thres), ]
            sig_rows <- rownames(sig_table)
            sig_full_m <- raw_df[sig_rows, ]

            final_df <- cbind(data_id=sig_rows, sig_full_m[, all_cols], sig_table)
            # final_df <- cbind(data_id=sig_rows, sig_full_m[, all_cols], sig_table, sig_full_annotation)

            if (!is.null(annot_key) && !is.null(annot_df)) {
                sig_annot_col <- sig_full_m[, annot_key]
                if (length(sig_annot_col) > 0) {
                    sig_annot_col_trimmed <- sapply(as.character(sig_annot_col), function(annot) { unlist(strsplit(annot, ";"))[1] })
                }
                else {
                    sig_annot_col_trimmed <- c()
                }
                sig_full_annotation <- annot_df[sig_annot_col_trimmed, ]
                final_df <- cbind(final_df, sig_full_annotation)
            }

            final_df
        }
    ),
    private = list(
        get_sig_rows = function(full_de_table, sig_thres, sig_col="adj.P.Val") {

            sig_table <- full_de_table[which(full_de_table[, sig_col] < sig_thres), ]
            sig_rows <- rownames(sig_table)
            sig_rows
        }
    )
)

comp <- SampleComp$new()
print("Loading module to 'comp'")
