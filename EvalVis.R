library(R6)
library(tidyverse)
library(reshape2)

EvalVis <- R6Class(
    public = list(
        sample_dist = function(data_m, title="Sample densities", alpha=0.001, color_col=NULL, max_count=NULL) {
            
            samples <- colnames(data_m)
            if (is.null(color_col)) {
                color_col <- samples
            }
            
            if (!is.null(max_count)) {
                data_m <- data_m %>% data.frame() %>% sample_n(max_count, replace=FALSE)
                title <- paste(title, "Subset:", max_count)
            }
            
            out_m <- data_m %>%
                t() %>%
                as.data.frame() %>%
                add_column(sample=samples, color=color_col, .before=1) %>%
                gather(row, value, -sample, -color)
            
            plt <- ggplot(out_m, aes(x=value, group=sample, color=color)) +
                geom_density(alpha=alpha, na.rm=TRUE) +
                ggtitle(title)
            private$style_plt(plt)
        },
        

        abundance_bars = function(data_m, title="Sample abundances", color_col=NULL, 
                                  show_missing=FALSE, show_average=FALSE, cond_order=FALSE) {

            if (!show_missing) {
                if (!show_average) {
                    values <- colSums(data_m, na.rm=TRUE)
                }
                else {
                    values <- colMeans(data_m, na.rm=TRUE)
                }
            }
            else {
                if (!show_average) {
                    values <- colSums(is.na(data_m))
                }
                else {
                    values <- colMeans(is.na(data_m))
                }
            }

            if (!is.null(color_col)) {

                plt_df <- data.frame(Sample=colnames(data_m), vals=values, level=color_col)

                if (cond_order) {
                    plt_df <- plt_df %>%
                        arrange(level) %>%
                        mutate(Sample=as.character(Sample))

                    plt_df$Sample <- reorder(plt_df$Sample, as.numeric(plt_df$level))
                }
                plt <- ggplot(plt_df, aes(x=Sample, y=vals, fill=level)) +
                    geom_bar(stat="identity")
            }
            else {
                plt_df <- data.frame(sample=colnames(data_m), vals=values)
                plt <- ggplot(plt_df, aes(x=sample, y=vals)) +
                    geom_bar(stat="identity")
            }

            if (!show_missing) {
                plt <- plt + ylab("Total intensity")
            }
            else {
                plt <- plt + ylab("Missing values")
            }

            plt <- plt + 
                theme_classic() + 
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + 
                ggtitle(title)
            
            plt
        },

        pvalhist = function(pvals, title="P value histogram", binwidth=0.01, vline=NULL, na.rm=FALSE) {
            
            plt <- ggplot(
                data.frame(pvals=pvals), 
                aes(pvals)) + 
                geom_histogram(binwidth=binwidth, na.rm=na.rm, fill="#268BD2") + 
                ggtitle(title) +
                xlim(-0.1, 1.1)

            if (!is.null(vline)) {
                plt <- plt + geom_vline(xintercept=vline)
            }
            
            plt <- plt + theme_classic() 

            plt
        },

        qq = function(df, title="QQ", max_count=NULL, cond_col=NULL) {
            
            if (!is.null(max_count)) {
                df <- df %>% data.frame() %>% sample_n(max_count, replace=FALSE)
                title <- paste(title, "Subset:", max_count)
            }
            
            df_data <- df %>% 
                data.frame() %>% 
                gather("id", "expression") %>% 
                mutate(id=as.character(id))
            
            if (!is.null(cond_col)) {
                cond_map <- cbind(id=colnames(df), cond=cond_col) %>% data.frame()
                print(head(df_data))
                print(head(cond_map))
                df_data <- df_data %>% inner_join(y=cond_map, by="id")
            }
            
            vec <- df_data$expression
            y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
            x <- qnorm(c(0.25, 0.75))
            slope <- diff(y) / diff(x)
            int <- y[1L] - slope * x[1L]

            if (!is.null(cond_col)) {
                print("color")
                base_plt <- ggplot(df_data, aes(sample=expression, color=cond))
            }
            else {
                print("non color")
                base_plt <- ggplot(df_data, aes(sample=expression))
            }
            
            qq_plt <- base_plt + geom_point(stat="qq", size=0.5, na.rm=TRUE) + 
                ggtitle(title) + 
                geom_abline(slope=slope, intercept=int)
                
            private$style_plt(qq_plt)
        }
    ),

    private = list(
        style_plt = function(plt) {
            palette <- "Dark2"
            plt <- plt + 
                theme_classic()
            # scale_color_brewer(palette=palette) + 
            # scale_fill_brewer(palette=palette)
            plt
        }
    )
)

evalvis <- EvalVis$new()
ev <- evalvis
print("Loading module to 'evalvis' and 'ev'")
