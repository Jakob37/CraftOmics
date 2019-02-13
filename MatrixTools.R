library(R6)
library(grid)

MatrixTools <- R6Class(
    public = list(
        calc_cond_rows = function(data_m, cond_col, filter_val, filter_all_na=F) {

            # Calculate number of partial or empty cases for condition

            target_cols <- cond_col == filter_val
            data_m_target <- data_m[, target_cols]

            if (filter_all_na) {
                nrow(na.omit(data_m_target))
            }
            else {
                has_value <- apply(data_m_target, 1, function(row) { length(na.omit(row)) > 0 })
                data_m_w_value <- data_m_target[has_value, ]
                nrow(data_m_w_value)
            }
        },
        get_conds_rows_df = function(data_m, cond_col, filter_vals) {

            # Generate matrix with number of complete, partial and no values for each provided condition

            vals_list <- list()
            name_list <- c()
            for (filter_val in filter_vals) {

                all_rows <- nrow(data_m)
                least_one_rows <- self$calc_cond_rows(data_m, cond_col, filter_val, filter_all_na=F)
                no_na_rows <- self$calc_cond_rows(data_m, cond_col, filter_val, filter_all_na=T)
                row_vals <- list(c(all=all_rows, w_value=least_one_rows, no_na=no_na_rows))
                vals_list <- c(vals_list, row_vals)
                name_list <- c(name_list, filter_val)
            }
            df <- as.data.frame(vals_list)
            colnames(df) <- name_list
            df
        },
        remove_decoys = function(full_df, annot_col, decoy_pattern) {
            is_decoy <- grepl(decoy_pattern, full_df[, annot_col])
            print(paste("Removing", length(which(is_decoy)), "decoys using pattern", decoy_pattern))
            full_df[!is_decoy,]
        },
        write_table = function(df, out_fp, with_rows=F) {
            if (!with_rows) {
                write.table(gen_table, file=out_fp, sep="\t", row.names=F)
            }
            else {
                write.table(gen_table, file=out_fp, sep="\t", col.names=NA)
            }
            print(paste(nrow(gen_table), "entries written to", out_fp))
            print(paste("Dimensions:", paste(dim(gen_table), collapse=" ")))

            print("First row:")
            head(gen_table, 1)
        }
    ),
    private = list()
)


mattool <- MatrixTools$new()
mt <- mattool
print("Loading module to 'mattool' and 'mt'")
