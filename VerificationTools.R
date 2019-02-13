library(R6)

VerificationTools <- R6Class(
    public = list(
        do_verification_plots = function(id_list, data_m, contrast_df, samples, cond_col, color_col=NULL, cols=4, decs=4) {

            scatters <- list()
            for (i in 1:length(id_list)) {

                id <- id_list[i]
                cond <- cond_col
                samples <- as.character(samples)
                intensities <- unlist(data_m[id, samples])
                q_val <- contrast_df[id, "adj.P.Val"]
                if (is.null(color_col)) {
                    plt <- qplot(cond, intensities,
                                 main=paste("FDR:", round(q_val, decs), "\nID: ", id),
                                 data=data.frame(cond=cond, intensities=intensities),
                                 ylab="intensity")
                }
                else {
                    plt <- qplot(cond, intensities,
                                 main=paste("FDR:", round(q_val, decs), "\nID: ", id),
                                 data=data.frame(cond=cond, intensities=intensities),
                                 ylab="intensity", color=color_col)
                }
                scatters[[i]] <- plt
            }
            scatters
        },

        get_random_names = function(df, count, sig_thres=NULL, pick_lower=T, sig_col="adj.P.Val") {

            if (!is.null(sig_thres)) {
                if (pick_lower) {
                    df <- df[which(df[, sig_col] < sig_thres), ]
                }
                else {
                    df <- df[which(df[, sig_col] > sig_thres), ]
                }
            }

            df_names <- rownames(df)

            index_picks <- sample(1:length(df_names), count, replace=F)
            name_picks <- df_names[index_picks]
            name_picks
        }
    ),
    private = list()
)

vertools <- VerificationTools$new()
vt <- vertools
print("Loading module to  'vertools' and 'vt'")
