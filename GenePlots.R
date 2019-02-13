library(R6)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggthemes)

GenePlots <- R6Class(
    public = list(
        prepare_long_df = function(smat, sample_annots_list, annot_col, source_name=NULL) {
            
            data_wide <- cbind(annot=annot_col, smat)
            data_long <- data_wide %>% data.frame() %>% gather("Sample", "Expression", -annot)
            data_long$Expression <- as.numeric(data_long$Expression)

            for (sample_annot_name in names(sample_annots_list)) {
                sample_annot <- sample_annots_list[[sample_annot_name]]
                sample_cond_map <- sample_annot
                names(sample_cond_map) <- colnames(smat)
                data_long[[sample_annot_name]] <- vapply(
                    data_long$Sample %>% unlist(),
                    function(sample, map) { map[[sample]] },
                    "",
                    map=sample_cond_map
                )
            }
            
            if (!is.null(source_name)) {
                data_long$Source <- source_name
            }

            data_long
        },

        make_stripplot = function(long_df, sample_col, value_col, color_col, split_col=NULL,
                                  custom_colors=NULL, nudges=c(-0.1, 0.1), title="Default title") {
        
            plt <- ggplot()

            if (!is.null(split_col)) {
                iteration <- 1

                for (split_level in unique(long_df[, split_col])) {
                    sub_df <- long_df[long_df[, split_col] == split_level, ]
                    sub_df$conds <- paste0(sub_df[, split_col], "_", sub_df[, color_col])

                    nudge_dist <- nudges[iteration]
                    plt <- plt + geom_jitter(
                        data=sub_df,
                        aes_string(x=sample_col, y=value_col, color="conds"),
                        position=position_nudge(x = nudge_dist)
                    )

                    iteration <- iteration + 1
                }
            }

            else {
                plt <- plt + geom_jitter(
                    long_df, 
                    aes_string(x=sample_col, y=value_col, color=color_col)
                )
            }

            if (!is.null(custom_colors)) {
                plt <- plt + scale_color_manual(values=custom_colors)
            }

            plt <- plt + theme(axis.text.x = element_text(angle = 90, hjust = 1))

            plt
        }
    ),
    private = list()
)


geneplots <- GenePlots$new()
gplt <- geneplots
print("Loading module to 'geneplots' and 'gplt'")
