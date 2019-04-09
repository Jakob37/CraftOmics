library(R6)
library(tidyverse)

library(ggpubr)

MasterWidgetPlotFuncs <- R6Class(
    
    public = list(
        
        annotate = function(plts, input) {
            
            if (input$plot_title != "") {
                title <- input$plot_title
            }
            else {
                title <- "A title can be assigned under 'Display'"
            }
            
            plot_rows <- ceiling(length(plts) / input$plot_cols)
            
            if (input$legend_color != "" || input$legend_fill != "") {
                plts <- lapply(
                    plts, 
                    function(plt, color, fill) { plt + labs(color=color, fill=fill) },
                    color=input$legend_color,
                    fill=input$legend_fill
                )
            }
            
            if (input$tabs != "Venns") {
                plts <- lapply(
                    plts,
                    function(plt, title_size, axis_size, ticks_size) { 
                        plt + theme(
                            plot.title=element_text(size=title_size),
                            axis.text=element_text(size=ticks_size),
                            axis.title=element_text(size=axis_size)
                        ) 
                    },
                    title_size=input$subtitle_size,
                    axis_size=input$axis_size,
                    ticks_size=input$ticks_size
                )
            }
            
            arr_obj <- ggpubr::ggarrange(
                plotlist=plts, 
                ncol=input$plot_cols, 
                nrow=plot_rows, 
                common.legend=input$legend_common,
                legend = input$legend_pos
            )
            
            if (input$legend_common) {
                arr_obj <- arr_obj 
            }
            
            subtext <- strwrap(input$plot_subtext, input$subtext_wrap, simplify=FALSE)
            subtext_parsed <- sapply(subtext, paste, collapse="\n")
            
            annotate_figure(
                arr_obj, 
                top = text_grob(title, color = "black", face = "bold", size = input$title_size),
                bottom = text_grob(subtext_parsed, color = "darkgray", size = input$title_size),
                fig.lab.size = input$subtitle_size
            )
        },
        
        do_bar = function(datasets, input, outlier_sets) {
            dobs <- self$get_preproc_list(datasets, c(input$data1, input$data2), input$checkgroup, outlier_sets)
            plts <- list()
            for (i in seq_len(2)) {
                plts[[i]] <- ev$abundance_bars(
                    dobs[[i]]$sdf, 
                    color_col=as.factor(dobs[[i]]$ddf[[input[[paste0("cond", i)]]]]), 
                    title=dobs[[i]]$title, 
                    show_missing=input$show_na, 
                    show_average=input$show_mean
                )
            }
            
            plts
        },
        
        do_density = function(datasets, input, outlier_sets) {
            dobs <- self$get_preproc_list(datasets, c(input$data1, input$data2), input$checkgroup, outlier_sets)
            if (input$fulldata) usecount <- NULL
            else usecount <- input$subset
            
            plts <- list()
            for (i in seq_len(2)) {
                plts[[i]] <- ev$sample_dist(
                    dobs[[i]]$sdf, 
                    color_col=as.factor(dobs[[i]]$ddf[[input[[paste0("cond", i)]]]]), 
                    title=dobs[[1]]$title, 
                    max_count=usecount
                )
            }
            
            list(plts[[1]], plts[[2]])
        },
        
        do_qq = function(datasets, input, outlier_sets) {
            dobs <- self$get_preproc_list(datasets, c(input$data1, input$data2), input$checkgroup, outlier_sets)
            if (input$fulldata) usecount <- NULL
            else usecount <- input$subset
            
            plts <- list()
            for (i in seq_len(2)) {
                plts[[i]] <- ev$qq(
                    dobs[[i]]$sdf, 
                    title=dobs[[i]]$title, 
                    cond_col=as.factor(dobs[[i]]$ddf[[input[[paste0("cond", i)]]]]), 
                    max_count=usecount
                )
            }
            
            list(plts[[1]], plts[[2]])
        },
        
        do_pca = function(datasets, input, outlier_sets) {
            
            dobs <- self$get_preproc_list(datasets, c(input$data1, input$data2), input$checkgroup, outlier_sets)
            if (!input$as_label) label <- NULL
            else label <- dobs[[1]]$ddf[, input$text_labels]
            
            plts <- list()
            for (i in seq_len(2)) {
                plts[[i]] <- mv$pca(
                    dobs[[i]]$sdf,
                    as.factor(dobs[[i]]$ddf[[input[[paste0("cond", i)]]]]),
                    pcs=c(as.numeric(input[[paste0("pc1_plt", i)]]), as.numeric(input[[paste0("pc2_plt", i)]])),
                    label=label
                ) + ggtitle(dobs[[i]]$title)
            }
            
            if (!input$pca_hide_loadings) {
                scree_plt <- mv$plot_component_fraction(dobs[[1]]$sdf, max_comps=input$pc_comps)
                plts[[length(plts)+1]] <- scree_plt
            }
            
            plts
        },
        
        do_cluster = function(datasets, input, outlier_sets) {
            
            dobs <- self$get_preproc_list(datasets, c(input$data1, input$data2), input$checkgroup, outlier_sets)
            if (!input$as_label) label <- NULL
            else label <- dobs[[1]]$ddf[, input$text_labels]
            
            plts <- list()
            for (i in seq_len(2)) {
                plts[[i]] <- mv$dendogram(
                    dobs[[i]]$sdf, 
                    as.factor(dobs[[i]]$ddf[[input[[paste0("cond", i)]]]])) + ggtitle(dobs[[i]]$title)
            }
            
            plts
        },
        
        do_hists = function(datasets, input, contrast_suffix, outlier_sets) {
            
            dobs <- self$get_preproc_list(datasets, input$stat_data, input$checkgroup, outlier_sets)
            contrasts <- private$get_contrasts_from_suffix(colnames(dobs[[1]]$adf), contrast_suffix)
            
            plts <- list()
            for (contrast in contrasts) {
                
                target_col <- paste(contrast, input$hist_target, sep=".")
                
                plt <- ev$pvalhist(dobs[[1]]$adf[[target_col]], na.rm=TRUE, bincount=input$hist_bins) +
                    theme_classic() +
                    scale_fill_brewer(palette="Dark2") +
                    ggtitle(contrast) +
                    xlab(input$hist_target) +
                    ylab("Count")
                
                plts[[contrast]] <- plt
            }
            
            if (input$scale_x_axis) {
                plts <- private$scale_axis(plts, "x")
            }
            
            if (input$scale_y_axis) {
                plts <- private$scale_axis(plts, "y")
            }
            
            plts
        },
        
        do_venns = function(datasets, input, contrast_suffix, outlier_sets) {
            
            dobs <- self$get_preproc_list(datasets, input$stat_data, input$checkgroup, outlier_sets)
            contrasts <- private$get_contrasts_from_suffix(colnames(dobs[[1]]$adf), contrast_suffix)
            
            if (length(contrasts) < 2) {
                stop("Must have at least two contrast-levels to show Venn")
            }
            
            if (!input$threeway_venn) {
                out <- stv$plot_comp_venns(
                    dobs[[1]]$adf, 
                    contrasts, 
                    base_sig_col = input$venn_type,
                    base_fold_name = input$venn_fold,
                    sig_thres = input$venn_thres,
                    check_greater_than = input$venn_inverse
                    # log2_fold_thres = input$fold
                )
                out
                # grid.arrange(grobs=out, ncol=3, top=paste0("Dataset: ", input$stat_data))
            }
            else {
                target_contrasts <- contrasts[1:3]
                if (length(target_contrasts) != length(contrasts)) {
                    warning("Performing threeway-venn only for first three contrasts")
                }
                plt <- stv$threeway_venn(
                    dobs[[1]]$adf,
                    contrasts,
                    thres_col_base=input$venn_type,
                    thres=input$venn_thres,
                    fold_col_base=input$venn_fold,
                    check_greater_than=input$venn_inverse
                )
                list(plt)
            }
        },
        
        do_general_scatter = function(datasets, input, contrast_suffix, outlier_sets) {
            
            dobs <- self$get_preproc_list(datasets, input$stat_data, input$checkgroup, outlier_sets)
            contrasts <- private$get_contrasts_from_suffix(colnames(dobs[[1]]$adf), contrast_suffix)
            
            plts <- list()
            for (contrast in contrasts) {
                
                target_cols <- paste(contrast, c(input$scatter_x, input$scatter_y, input$scatter_color), sep=".")
                
                x_col <- target_cols[1]
                y_col <- target_cols[2]
                sig_col <- target_cols[3]
                
                tbl <- dobs[[1]]$adf[, target_cols]
                tbl <- tbl[complete.cases(tbl), ]
                
                if (input$scatter_minuslog_y) {
                    tbl[[y_col]] <- -log10(tbl[[y_col]])
                }
                
                tbl$is_sig <- tbl[, sig_col] < input$scatter_color_cutoff
                tbl <- tbl %>% arrange(UQ(as.name(sig_col)))
                
                plt <- ggplot(tbl, aes_string(x_col, y_col, color="is_sig")) +
                    geom_point(size=1, alpha=0.5, na.rm=TRUE) +
                    ggtitle(paste("Dataset:", input$stat_data)) +
                    theme_classic()
                
                plts[[contrast]] <- plt
            }
  
            
            if (input$scalexaxis) {
                plts <- private$scale_axis(plts, "x")
            }
            
            if (input$scaleyaxis) {
                plts <- private$scale_axis(plts, "y")
            }
            
            plts
        },
        
        do_spotcheck = function(datasets, input) {
            ggplot() + ggtitle("Spotcheck currently not implemented") + theme_classic()
        },
        
        # do_table = function(datasets, input) {
        #     ggplot() + ggtitle("Table currently not implemented") + theme_classic()
        #     
        #     thedata <- reactive({
        #         
        #         retained <- rowData(stat_data[[input$data]]) %>% data.frame()
        #         if (length(input$filters) > 0) {
        #             
        #             if (input$exclusive) {
        #                 unique_retained <- retained
        #                 
        #                 # Only include features passing all filters
        #                 for (filter in input$filters) {
        #                     unique_retained <- unique_retained %>% filter(UQ(as.name(filter)) < input$filterthres)
        #                 }
        #                 retained <- unique_retained
        #             }
        #             else {
        #                 all_retained <- NULL
        #                 for (filter in input$filters) {
        #                     # Include features passing at least one filter
        #                     filter_retained <- retained %>% filter(UQ(as.name(filter)) < input$filterthres)
        #                     all_retained <- rbind(all_retained, filter_retained)
        #                 }
        #                 retained <- all_retained %>% distinct()
        #             }
        #         }
        #         
        #         filtered_selected <- retained %>% 
        #             dplyr::select(input$fields) %>%
        #             data.frame()
        #         
        #         filtered_selected
        #     })
        #     
            # observe({
            #     new_choices <- colnames(rowData(stat_data[[input$data]]))
            #     updateSelectInput(
            #         session,
            #         "fields",
            #         choices=new_choices,
            #         selected=input$fields
            #     )
            # })
            # 
        #     output$table = DT::renderDataTable({
        #         
        #         thedata() %>% 
        #             datatable(options=list(
        #                 pageLength=10, 
        #                 scrollX=TRUE,
        #                 # autoWidth=TRUE,
        #                 columnDefs=list(list(width="10px", targets="_all"))
        #             )) %>%
        #             DT::formatRound(columns=input$fields, digits=input$decimals)
        #     })
        # },
        
        do_table = function(datasets, input, outlier_sets) {

            
            
            dobs <- self$get_preproc_list(datasets, input$stat_data, input$checkgroup, outlier_sets)
            filter_adf <- dobs[[1]]$adf
            
            all_filters <- unique(c(
                input$table_filters, 
                unique(private$get_contrast_fields(filter_adf, input$contrast_filters))
            ))
            
            if (length(all_filters) > 0) {
                filter_adf <- private$multi_filter_table(filter_adf, all_filters, input$table_filterthres, input$table_filter_less_than, exclusive=input$table_exclusive_filter, verbose=TRUE)
            }

            select_fields <- unique(c(
                input$table_fields, 
                unique(private$get_contrast_fields(filter_adf, input$table_contrast_fields))
            ))

            target_adf <- filter_adf %>% dplyr::select(select_fields) %>% data.frame() 
            target_adf
        },
        
        do_profile = function(datasets, input) {
            ggplot() + ggtitle("Profile currently not implemented") + theme_classic()
        },
        
        get_preproc_list = function(datasets, dataset_names, checkgroup, outlier_sets) {
            
            outliers <- unname(unlist(outlier_sets[checkgroup]))
            
            dataset_objs <- list()
            for (i in seq_len(length(dataset_names))) {
                
                dataset_name <- dataset_names[i]
                dataset_obj <- self$parse_dataset(datasets[[dataset_name]], outliers)
                dataset_obj$title <- paste("Dataset:", dataset_name)
                dataset_obj$outliers <- outliers
                dataset_objs[[i]] <- dataset_obj
            }
            
            dataset_objs
        },
        
        parse_dataset = function(dataset, outliers=NULL, target_assay=1) {
            
            if (typeof(target_assay) == "character" && !(target_assay %in% names(assays(dataset)))) {
                stop("Unknown character target_assay: ", target_assay)
            }
            
            if (!is.null(outliers)) {
                non_outliers <- colnames(dataset)[!colnames(dataset) %in% outliers]
            }
            else {
                non_outliers <- colnames(dataset)
            }
            
            sdf <- assays(dataset)[[target_assay]][, non_outliers]
            ddf <- data.frame(colData(dataset)) %>% filter(sample %in% non_outliers)
            adf <- rowData(dataset) %>% data.frame()
            list("sdf"=sdf, "ddf"=ddf, "adf"=adf)
        }
    ),
    private = list(
        get_contrasts_from_suffix = function(row_data_cols, contrast_suffix, include_suffix=FALSE) {
            
            contrasts <- make.names(gsub(
                paste0(".", contrast_suffix), "", 
                row_data_cols[grepl(paste0(contrast_suffix, "$"), row_data_cols)]))
            if (include_suffix) {
                contrasts <- paste(contrasts, contrast_suffix, sep=".")
            }
            contrasts
        },
        
        scale_axis = function(plts, axis, padding_fraction=0.01) {
            
            if (axis == "x") {
                lim_func <- ggplot2::xlim
            }
            else if (axis == "y") {
                lim_func <- ggplot2::ylim
            }
            else {
                stop("Unknown axis: '", axis, "', only allowed are: 'x' and 'y'")
            }
            
            min_vals <- lapply(plts, function(plt) {
                layer_scales(plt)[[axis]]$range$range[1]
            })
            max_vals <- lapply(plts, function(plt) {
                layer_scales(plt)[[axis]]$range$range[2]
            })
            plts <- lapply(plts, function(plt) {
                min_val <- min(unlist(min_vals))
                max_val <- max(unlist(max_vals))
                plt <- plt + lim_func(min_val - 0.01 * abs(min_val), max_val + 0.01 * abs(max_val))
            })
            plts
        },
        
        multi_filter_table = function(raw_df, filters, filter_thres, filter_less_than, exclusive=FALSE, verbose=TRUE) {
            
            raw_df$temp_id <- seq_len(nrow(raw_df))
            filter_adf <- raw_df
            
            if (verbose) {
                message("Before filtering: ", nrow(filter_adf), " rows")            
            }
            
            for (filter in filters) {
                
                non_exclusive_ids <- list()
                
                if (exclusive) {
                    filter_adf <- private$do_filter(filter_adf, filter, filter_thres, greater_than=!filter_less_than)
                    if (verbose) {
                        message("Filtering: ", filter, " retained rows: ", nrow(filter_adf))
                    }
                }
                else {
                    
                    temp_filter_adf <- private$do_filter(raw_df, filter, filter_thres, greater_than=filter_less_than)
                    non_exclusive_ids[[filter]] <- temp_filter_adf$temp_id
                }
                
            }
            
            if (!exclusive) {
                all_non_exclusive_ids <- sort(unique(unlist(non_exclusive_ids)))
                raw_df %>% dplyr::filter(temp_id %in% all_non_exclusive_ids) %>% dplyr::select(-temp_id)
            }
            else {
                filter_adf %>% dplyr::select(-temp_id)
            }
        },
        
        do_filter = function(df, field, value, greater_than=FALSE) {
            if (greater_than) {
                df %>% dplyr::filter(UQ(as.name(field)) > value)
            }
            else {
                df %>% dplyr::filter(UQ(as.name(field)) < value)
            }
        },
        
        get_contrast_fields = function(adf, contrast_filters) {
            fields <- c()
            for (contrast_suffix in contrast_filters) {
                contrast_fields <- private$get_contrasts_from_suffix(colnames(adf), contrast_suffix, include_suffix=TRUE)
                fields <- unique(c(fields, contrast_fields))
            }
            fields
        }
    )
)