library(limma)
library(tidyverse)
library(R6)
library(wCorr)

StatTools <- R6Class(
    public = list(

        calc_pearson_table = function(df, ref_levels, min_vals=5, weights=NULL, retain_status=FALSE) {

            corr_val_table <- apply(df, 1, function(row) {
                non_na <- which(!is.na(row))
                # Too few valid values
                if (length(non_na) < min_vals) {
                    corr_val <- NA
                    corr_pval <- NA
                    wcorr_val <- NA
                    status <- "too_few"
                }
                # Outcome variable only has one level
                else if (length(unique(ref_levels[non_na])) == 1) {
                    corr_val <- NA
                    corr_pval <- NA
                    wcorr_val <- NA
                    status <- "one_ref"
                }
                else {
                    corr_val <- cor(x=row[non_na], y=ref_levels[non_na], method="pearson", use="complete.obs")
                    corr_pval <- cor.test(x=row[non_na], y=ref_levels[non_na], method="pearson", use="complete.obs")$p.value
                    if (!is.null(weights)) {
                        wcorr_val <- wCorr::weightedCorr(x=row[non_na], y=ref_levels[non_na], method="pearson", weights=weights[non_na])
                    }
                    status <- "OK"
                }
                
                if (is.null(weights)) {
                    c(corr_val, corr_pval, status)
                }
                else {
                    c(corr_val, wcorr_val, corr_pval, status)
                }
            })
                        
            corr_val_df <- data.frame(t(corr_val_table))
            if (is.null(weights)) {
                colnames(corr_val_df) <- c("pearson_corr", "pearson_p", "status")
            }
            else {
                colnames(corr_val_df) <- c("pearson_corr", "pearson_w_corr", "pearson_p", "status")
            }

            corr_val_df$pearson_q <- p.adjust(corr_val_df$pearson_p, method = "BH")
            
            message("Calculation status:")
            print(table(corr_val_df$status))

            if (!retain_status) {
                corr_val_df[, !(colnames(corr_val_df) %in% "status")]
            }
            else {
                corr_val_df
            }
        },

        reduce_technical_replicates_for_matrices = function(dataMat, designMat, techRepGroups) {
            
            uniqueGroups <- unique(techRepGroups)
            indices <- lapply(uniqueGroups, function(i) { which(techRepGroups %in% i) })

            collData <- lapply(
                indices, 
                function(inds) { rowMeans(dataMat[, inds, drop=FALSE], na.rm = TRUE) })
            collDataMat <- as.matrix(data.frame(collData))
            colnames(collDataMat) <- uniqueGroups
            first_indices <- sapply(indices, function(ind_group){ind_group[1]})
            collDesignMat <- designMat[first_indices, ]
            collDesignMat$sample <- uniqueGroups

            list(data=collDataMat, design=collDesignMat)
        },

        # Filters out rows where a given number of values isn't present
        # in all replicate groups
        filter_low_rep = function(df, groups, least_rep=2) {

            # row_meet_thres_contrast <- apply(df, 1, self$all_replicates_have_values, groups=groups, min_count=least_rep)
            row_meet_thres_contrast <- self$filter_low_rep_contrast(df, groups, least_rep)
            filtered_df <- df[row_meet_thres_contrast, ]
            filtered_df
        },

        filter_low_rep_contrast = function(df, groups, least_rep=2) {
            
            row_meet_thres_contrast <- apply(df, 1, self$all_replicates_have_values, groups=groups, min_count=least_rep)
            row_meet_thres_contrast
        },
        
        all_replicates_have_values = function(row, groups, min_count) {
            names(row) <- groups
            rep_counts <- table(names(na.omit(row)))
            length(rep_counts) == length(unique(groups)) && min(rep_counts) >= min_count
        },
        
        contrasts = function(data_m, conditions, contrasts, verbose=FALSE, batch_factor=NULL) {
            
            stat_tables <- lapply(contrasts, function(contr){
                if (verbose) print(paste("Calculating: ", paste(contr, collapse=", ")))
                tabl <- self$contrast(data_m, conditions, contr, verbose=verbose, batch_factor=batch_factor)
                tabl
            })
        },

        contrast = function(data_m, level_factor, target_levels, type="limma", batch_factor=NULL, 
                            filter_factor=NULL, filter_level=NULL, verbose=FALSE) {

            if (length(c(filter_factor, filter_level)) == 1) {
                stop("Need to assign either neither or both variables 'filter_factor' and 'filter_level'")
            }

            if (!is.null(filter_factor)) {
                target_samples <- filter_factor == filter_level
                data_m <- data_m[, target_samples]
                level_factor <- level_factor[target_samples]
            }
            
            Variable <- factor(level_factor)
            contrast_levels <- target_levels

            if (is.null(batch_factor)) {
                model <- ~0+Variable
            }
            else {
                Batch <- factor(batch_factor)
                model <- ~0+Variable+Batch
            }

            my_contrast <- paste0("Variable", target_levels[1], "-", "Variable", target_levels[2])
            design <- model.matrix(model)

            if (type == "limma") {
                stats_table <- self$calculate_limma_table(data_m, design, my_contrast)
            }
            else if (type %in% c("anova", "kruskal-wallis")) {
                warning("ANOVA/Kruskal-Wallis mode is not thoroughly debugged yet!")
                if (!is.null(batch_factor)) {
                    stats_table <- self$calculate_anova_table(data_m, Variable, contrast_levels, Batch=Batch, type=type)
                }
                else {
                    stats_table <- self$calculate_anova_table(data_m, Variable, contrast_levels, type=type)
                }
            }
            else {
                stop(paste("Unknown type:", type))
            }

            if (verbose) {
                st$summarize(stats_table)
            }
            
            stats_table
        },

        calculate_anova_table = function(data_m, Variable, contrast_levels, Batch=NULL, type=type) {

            p_values <- apply(data_m, 1, private$single_anova, contrast_levels=contrast_levels, Variable=Variable, Batch=Batch, type=type)
            fdr_values <- p.adjust(p_values, "BH")

            ave_expr <- apply(data_m, 1, mean, na.rm=T)
            row_var <- apply(data_m, 1, var, na.rm=T)
            na_count <- apply(data_m, 1, function(row) { sum(is.na(unlist(row))) })

            sign_table <- cbind(AveExpr=ave_expr, Variance=row_var, FullRowNAs=na_count, P.Value=p_values, adj.P.Val=fdr_values)
            rownames(sign_table) <- rownames(data_m)
            as.data.frame(sign_table)
        },

        calculate_limma_table = function(data_m, design, my_contrast, p_sorted=F) {

            fit <- lmFit(data_m, design)
            contrast.matrix <- makeContrasts(contrasts=c(my_contrast), levels=design)
            fit_contrasts <- contrasts.fit(fit, contrast.matrix)
            fit_bayes <- eBayes(fit_contrasts)

            limma_table <- topTable(fit_bayes, coef=1, number=Inf, sort="none")
            limma_table$row_nbr <- seq_len(nrow(limma_table))
            
            if (p_sorted) {
                limma_table <- limma_table %>% arrange(P.Value)
            }

            limma_table
        },

        summarize = function(ltab, label=NULL, p_col_name="P.Value", q_col_name="adj.P.Val") {

            p005 <- ltab[which(ltab[, p_col_name] < 0.05),]
            q01 <- ltab[which(ltab[, q_col_name] < 0.1),]

            cat(paste("Features p < 0.05:", nrow(p005),
                      "features q < 0.1:", nrow(q01),
                      "total features:", nrow(ltab),
                      label, "\n", sep="\t"))
        },

        variance_reduce = function(data_m, fraction) {

            # Get starting order of row names
            orig_order <- rownames(data_m)

            # Calculate variance for each row
            vars <- apply(data_m, 1, function(row) {var(unlist(row), na.rm=T)})
            names(vars) <- orig_order

            # Get variance sorted in falling size
            sort_vars <- sort(vars, decreasing=T)

            # Set cutting position based on fraction, rounding up
            cut_pos <- as.integer(fraction * length(vars)) + 1

            # Get names for entries with largest variance
            big_var_names <- names(sort_vars[1:cut_pos])

            # Extract IDs in original order (rather than sorting, which would reorder unnumbered rows)
            big_var_names_orig_order <- which(orig_order %in% big_var_names)

            # Extract the subset with original ordering
            var_subset_m <- data_m[big_var_names_orig_order, ]

            var_subset_m
        }
        
    ),
    private = list(

        single_anova = function(values, contrast_levels, Variable, Batch=NULL, type="anova") {

            if (!self$all_replicates_have_values(unlist(values), Variable, 1)) {
                return(NaN)
            }

            if (!is.null(Batch)) {
                row_df <- data.frame(Data=unlist(values), Variable=Variable, Batch=Batch)

                test_levels <- levels(Variable)
                MyComp <- c(rep(0, length(test_levels)))
                MyComp[which(test_levels == contrast_levels[1])] <- -1
                MyComp[which(test_levels == contrast_levels[2])] <- 1

                mat <- cbind(MyComp)
                contrasts(row_df$Variable) <- mat

                if (type == "anova") {
                    my_model <- lm(Data ~ Variable + Batch, data=row_df)
                    summary(my_model)$coefficients["VariableMyComp", "Pr(>|t|)"]
                }
                else if (type == "kruskal-wallis") {
                    my_model <- kruskal.test(Data ~ Variable + Batch, data=row_df)
                    my_model$p.value
                }
                else {

                }
            }
            else {
                row_df <- data.frame(Data=unlist(values), Variable=Variable)

                test_levels <- levels(Variable)
                MyComp <- c(rep(0, length(test_levels)))
                MyComp[which(test_levels == contrast_levels[1])] <- -1
                MyComp[which(test_levels == contrast_levels[2])] <- 1

                mat <- cbind(MyComp)
                contrasts(row_df$Variable) <- mat

                if (type == "anova") {
                    my_model <- lm(Data ~ Variable, data=row_df)
                    summary(my_model)$coefficients["VariableMyComp", "Pr(>|t|)"]
                }
                else if (type == "kruskal-wallis") {
                    my_model <- kruskal.test(Data ~ Variable, data=row_df)
                    my_model$p.value
                }
                else {
                    stop(paste("Unknown stats type:", type))
                }
            }
        }
    )
)

stattools <- StatTools$new()
st <- stattools
print("Loading module to 'stattools' and 'st'")
