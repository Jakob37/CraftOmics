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
                title <- paste("Dataset:", input$data1)
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
            
            annotate_figure(
                arr_obj, 
                top = text_grob(title, color = "black", face = "bold", size = input$title_size),
                fig.lab.size = input$subtitle_size
            )
        },
        
        do_bar = function(datasets, input) {
            dobs <- self$get_preproc_list(datasets, c(input$data1, input$data2), input$checkgroup)
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
        
        do_density = function(datasets, input) {
            dobs <- self$get_preproc_list(datasets, c(input$data1, input$data2), input$checkgroup)
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
        
        do_qq = function(datasets, input) {
            dobs <- self$get_preproc_list(datasets, c(input$data1, input$data2), input$checkgroup)
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
        
        do_pca = function(datasets, input) {
            
            dobs <- self$get_preproc_list(datasets, c(input$data1, input$data2), input$checkgroup)
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
            scree <- mv$plot_component_fraction(dobs[[1]]$sdf, max_comps=input$pc_comps)
            
            list(plts[[1]], plts[[2]], scree)
            
        },
        
        do_cluster = function(datasets, input) {
            
            dobs <- self$get_preproc_list(datasets, c(input$data1, input$data2), input$checkgroup)
            if (!input$as_label) label <- NULL
            else label <- dobs[[1]]$ddf[, input$text_labels]
            
            plts <- list()
            for (i in seq_len(2)) {
                plts[[i]] <- mv$dendogram(
                    dobs[[i]]$sdf, 
                    as.factor(dobs[[i]]$ddf[[input[[paste0("cond", i)]]]])) + ggtitle(dobs[[i]]$title)
            }
            
            plts
            # grid.arrange(plts[[1]], plts[[2]], ncol=2)
        },
        
        do_hists = function(datasets, input, contrast_suffix) {
            
            dobs <- self$get_preproc_list(datasets, input$stat_data, input$checkgroup)
            contrasts <- private$get_contrasts_from_suffix(colnames(dobs[[1]]$adf), contrast_suffix)
            
            plts <- list()
            for (contrast in contrasts) {
                
                target_col <- paste(contrast, input$hist_target, sep=".")
                
                plt <- ev$pvalhist(dobs[[1]]$adf[[target_col]], na.rm=TRUE, bincount=input$hist_bins) +
                    theme_classic() +
                    scale_fill_brewer(palette="Dark2") +
                    ggtitle(contrast) 
                
                if (input$limitxaxis) {
                    xmin <- min(as.numeric(dobs[[1]]$adf[[target_col]]), na.rm = TRUE) - 0.01
                    xmax <- max(as.numeric(dobs[[1]]$adf[[target_col]]), na.rm = TRUE) + 0.01
                    plt <- plt + xlim(xmin, xmax)
                }
                
                plts[[contrast]] <- plt
            }
            
            if (input$scaleaxis) {
                max_y_vals <- lapply(plts, function(plt) {
                    layer_scales(plt)$y$range$range[2]
                })
                
                plts <- lapply(plts, function(plt) {
                    plt <- plt + 
                        ylim(0, 1.01 * max(unlist(max_y_vals)))
                })
            }
            
            plts
            # grid.arrange(grobs=plts, ncol=1, top=dobs[[1]]$title)
        },
        
        do_venns = function(datasets, input, contrast_suffix) {
            dobs <- self$get_preproc_list(datasets, input$stat_data, input$checkgroup)
            contrasts <- private$get_contrasts_from_suffix(colnames(dobs[[1]]$adf), contrast_suffix)
            
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
                stv$threeway_venn(
                    dobs[[1]]$adf,
                    contrasts,
                    thres_col_base=input$venn_type,
                    thres=input$venn_thres,
                    fold_col_base=input$venn_fold,
                    check_greater_than=input$venn_inverse
                )
                
            }
            
        },
        
        do_ma = function(datasets, input, contrast_suffix) {
            ggplot() + ggtitle("MA currently not implemented") + theme_classic()
        },
        
        do_vulc = function(datasets, input, contrast_suffix) {
            ggplot() + ggtitle("Vulcano currently not implemented") + theme_classic()
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
        
        do_table = function(datasets, input) {
            table <- datasets[[input$data1]] %>% rowData() %>% data.frame()
            DT::renderDataTable({
                table
            })
        },
        
        do_profile = function(datasets, input) {
            ggplot() + ggtitle("Profile currently not implemented") + theme_classic()
        },
        
        get_preproc_list = function(datasets, dataset_names, checkgroup) {
            
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
        get_contrasts_from_suffix = function(row_data_cols, contrast_suffix) {
            
            make.names(gsub(
                paste0(".", contrast_suffix), "", 
                row_data_cols[grepl(paste0(contrast_suffix, "$"), row_data_cols)]))
        }
    )
)