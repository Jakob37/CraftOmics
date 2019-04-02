MasterWidgetPlotFuncs <- R6Class(
    
    public = list(
        
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
            
            grid.arrange(plts[[1]], plts[[2]])
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
            
            grid.arrange(plts[[1]], plts[[2]])
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
            
            grid.arrange(plts[[1]], plts[[2]])
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
            
            grid.arrange(plts[[1]], plts[[2]], scree, ncol=2)
            
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
            
            grid.arrange(plts[[1]], plts[[2]], ncol=2)
        },
        
        do_hists = function(datasets, input, contrast_suffix) {
            
            dobs <- self$get_preproc_list(datasets, c(input$data1, input$data2), input$checkgroup)
            contrasts <- private$get_contrasts_from_suffix(colnames(dobs[[1]]$adf), contrast_suffix)
            
            plts <- list()
            for (contrast in contrasts) {
                
                target_col <- paste(contrast, input$target_col, sep=".")
                
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
            
            grid.arrange(grobs=plts, ncol=1, top=dobs[[1]]$title)
        },
        
        do_venns = function(datasets, input, contrast_suffix) {
            dobs <- self$get_preproc_list(datasets, c(input$data1, input$data2), input$checkgroup)
            contrasts <- private$get_contrasts_from_suffix(colnames(dobs[[1]]$adf), contrast_suffix)
            
            out <- stv$plot_comp_venns(
                dobs[[1]]$adf, 
                contrasts, 
                base_sig_col = input$type,
                base_fold_name = input$fold_base,
                sig_thres = input$thres, 
                log2_fold_thres = input$fold
            )
            grid.arrange(grobs=out, ncol=3, top=paste0("Dataset: ", input$data))
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