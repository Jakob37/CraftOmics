
library(shiny)

# Datasets expected format: Named list linked to SummarizedExperiment instances

MyWidgets <- R6Class(
    public = list(
        total_intensity_widget = function(name, datasets, outlier_sets, height=800, default_cond=NULL) {
            
            dataset_names <- names(datasets)
            default_name <- dataset_names[1]
            dataset <- datasets[[1]]
            
            if (is.null(default_cond)) {
                default_cond <- colnames(colData(dataset))[1]
            }
            
            shinyApp(
                ui = fluidPage(
                    tags$head(tags$style(HTML(".shiny-split-layout > div { overflow: visible; }"))),
                    splitLayout(
                        selectInput("data", "Dataset:", selected=default_name, choices=dataset_names),
                        selectInput("cond", "Condition:", choices = colnames(colData(dataset)), selected=default_cond)
                    ),
                    splitLayout(
                        checkboxGroupInput("checkgroup", "Remove group", choices=names(outlier_sets), selected=NULL),
                        checkboxInput("show_na", "Show NA counts:", value=FALSE)
                    ),
                    plotOutput("plot")
                ),
                server = function(input, output) {
                    output$plot = renderPlot({
                        
                        dataset <- datasets[[input$data]]
                        outliers <- unname(unlist(outlier_sets[input$checkgroup]))
                        parsed <- private$parse_dataset(dataset, outliers)

                        plt <- ev$abundance_bars(
                            parsed$sdf, 
                            color_col=as.factor(parsed$ddf[[input$cond]]), 
                            title=name, 
                            show_missing=input$show_na)
                        plt
                    })
                },
                options=list(height=height)
            )
        },
        
        qq_widget = function(name, datasets, outlier_sets, height=800, default_cond=NULL) {
            
            dataset_names <- names(datasets)
            default_name <- dataset_names[1]
            dataset <- datasets[[1]]
            
            if (is.null(default_cond)) {
                default_cond <- colnames(colData(dataset))[1]
            }
            
            shinyApp(
                ui = fluidPage(
                    tags$head(tags$style(HTML(".shiny-split-layout > div { overflow: visible; }"))),
                    splitLayout(
                        selectInput("data", "Dataset:", selected=default_name, choices=dataset_names),
                        selectInput("cond", "Condition:", choices = colnames(colData(dataset)), selected=default_cond)
                    ),
                    splitLayout(
                        checkboxGroupInput("checkgroup", "Remove group", choices=names(outlier_sets), selected=NULL),
                        checkboxInput("show_na", "Show NA counts:", value=FALSE)
                    ),
                    numericInput("subset", "Partial data:", value=1000, min=100, max=nrow(assay(dataset))),
                    plotOutput("plot")
                ),
                server = function(input, output) {
                    output$plot = renderPlot({
                        
                        dataset <- datasets[[input$data]]
                        outliers <- unname(unlist(outlier_sets[input$checkgroup]))
                        parsed <- private$parse_dataset(dataset, outliers)
                        
                        plt <- ev$qq(
                            parsed$sdf, 
                            title=name, 
                            max_count=input$subset, 
                            cond_col=parsed$ddf[[input$cond]]
                        )
                        
                        plt
                    })
                },
                options=list(height=height)
            )
        },
        
        pca_widget = function(name, datasets, outlier_sets, height=1300, default_cond1=NULL, default_cond2=NULL) {
            
            dataset_names <- names(datasets)
            default_name <- dataset_names[1]
            dataset <- datasets[[1]]
            
            if (is.null(default_cond1)) {
                default_cond1 <- colnames(colData(dataset))[1]
            }
            if (is.null(default_cond2)) {
                default_cond2 <- colnames(colData(dataset))[2]
            }
            
            shinyApp(
                ui = fluidPage(
                    tags$head(tags$style(HTML(".shiny-split-layout > div { overflow: visible; }"))),
                    splitLayout(
                        selectInput("data", "Dataset:", selected=default_name, choices=dataset_names),
                        checkboxGroupInput("checkgroup", "Remove group", choices=names(outlier_sets), selected=NULL)
                    ),
                    checkboxInput("as_label", "Show as text"),
                    splitLayout(
                        selectInput("pc1_plt1", "PC1 (plot1):", choices=1:10, selected=1),
                        selectInput("pc1_plt2", "PC1 (plot2):", choices=1:10, selected=1)
                    ),
                    splitLayout(
                        selectInput("pc2_plt1", "PC2 (plot1):", choices=1:10, selected=2),
                        selectInput("pc2_plt2", "PC2 (plot2):", choices=1:10, selected=2)
                    ),
                    splitLayout(
                        selectInput("cond_plt1", "Condition (plot1):", selected=default_cond1, choices = colnames(colData(dataset))),
                        selectInput("cond_plt2", "Condition (plot2):", selected=default_cond2, choices = colnames(colData(dataset)))
                    ),
                    plotOutput("pca")
                ),
                server = function(input, output) {
                    output$pca = renderPlot({
                        
                        dataset <- datasets[[input$data]]
                        outliers <- unname(unlist(outlier_sets[input$checkgroup]))
                        parsed <- private$parse_dataset(dataset, outliers)
                        
                        title1 <- paste0("Cond: ", input$cond_plt1, " PCs: ", input$pc1_plt1, ", ", input$pc2_plt1)
                        plt1 <- mv$pca(
                            parsed$sdf, 
                            as.factor(parsed$ddf[[input$cond_plt1]]), 
                            pcs=c(as.numeric(input$pc1_plt1), as.numeric(input$pc2_plt1)), 
                            label=input$as_label) + 
                                ggtitle(title1)
                        
                        title2 <- paste0("Cond: ", input$cond_plt2, " PCs: ", input$pc1_plt2, ", ", input$pc2_plt2)
                        plt2 <- mv$pca(
                            parsed$sdf, 
                            as.factor(parsed$ddf[[input$cond_plt2]]), 
                            pcs=c(as.numeric(input$pc1_plt2), as.numeric(input$pc2_plt2)), 
                            label=input$as_label) + 
                                ggtitle(title2)
                        
                        scree <- mv$plot_component_fraction(parsed$sdf)
                        
                        grid.arrange(plt1, plt2, scree, ncol=2)
                    }, height = 800)
                },
                options = list(height=height)
            )
        },
        
        clustering_widget = function(name, datasets, outlier_sets, height=1500, default_cond=NULL) {
            
            plot_height <- 1200
            
            dataset_names <- names(datasets)
            default_name <- dataset_names[1]
            dataset <- datasets[[1]]
            
            if (is.null(default_cond)) {
                default_cond <- colnames(colData(dataset))[1]
            }
            
            shinyApp(
                ui = fluidPage(
                    tags$head(tags$style(HTML(".shiny-split-layout > div { overflow: visible; }"))),
                    splitLayout(
                        selectInput("data", "Dataset:", selected=default_name, choices=dataset_names),
                        selectInput("cond", "Condition:", choices = colnames(colData(dataset)), selected=default_cond)
                    ),
                    checkboxGroupInput("checkgroup", "Remove group", choices=names(outlier_sets), selected=NULL),
                    plotOutput("plot", width="100%")
                ),
                server = function(input, output) {
                    output$plot = renderPlot({
                        
                        dataset <- datasets[[input$data]]
                        outliers <- unname(unlist(outlier_sets[input$checkgroup]))
                        parsed <- private$parse_dataset(dataset, outliers)
                        mv$dendogram(parsed$sdf, as.factor(parsed$ddf[[input$cond]]))
                        
                    }, height = plot_height)
                },
                options = list(height=height)
            )
        },
        
        phist_widget = function(name, stat_data, contrasts, p_col="P.Value", q_col="adj.P.Val", height=700) {
            
            p_cols <- paste(contrasts, p_col, sep=".")
            dataset_names <- names(stat_data)
            default_name <- dataset_names[1]
            dataset <- stat_data[[1]]
            
            shinyApp(
                ui = fluidPage(
                    selectInput("data", "Dataset:", selected=default_name, choices=dataset_names),
                    numericInput("bin", "P-hist bin width", value=0.02, min=0.0001, max=0.1),
                    checkboxInput("fdr", "Q-value histogram"),
                    splitLayout(
                        checkboxInput("vline", "Draw vertical line"),
                        sliderInput("vline_val", "Vertical line position", min=0, max=1, step=0.01, value=0.1)
                    ),
                    checkboxInput("scaleaxis", "Scale Y axes"),
                    plotOutput("plots", height=200)
                ),
                server = function(input, output) {
                    
                    output$plots = renderPlot({
                        
                        rdf <- data.frame(rowData(stat_data[[input$data]]))
                        
                        plts <- list()
                        for (contrast in contrasts) {
                            if (!input$fdr) {
                                target_col <- paste(contrast, p_col, sep=".")
                            }
                            else {
                                target_col <- paste(contrast, q_col, sep=".")
                            }

                            if (input$vline) {
                                vline <- input$vline_val
                            }
                            else {
                                vline <- NULL
                            }
                            
                            plt <- ev$pvalhist(rdf[[target_col]], na.rm=TRUE, binwidth=input$bin, vline=vline) +
                                theme_classic() +
                                scale_fill_brewer(palette="Dark2") +
                                ggtitle(contrast)
                            plts[[contrast]] <- plt
                        }
                        
                        if (input$scaleaxis) {
                            max_vals <- lapply(plts, function(plt) {
                                layer_scales(plt)$y$range$range[2]
                            })
                            plts <- lapply(plts, function(plt) {
                                plt <- plt + ylim(0, 1.01 * max(unlist(max_vals)))
                            })
                        }
                        
                        grid.arrange(grobs=plts, ncol=1)
                    }, height = 500)
                },
                options=list(height=height)
            )
        },
        
        scatter_widgets = function(name, stat_data, contrasts, p_col="P.Value", q_col="adj.P.Val", 
                                   fold_col="logFC", expr_col="AveExpr", height=1100) {
            
            p_cols <- paste(contrasts, p_col, sep=".")
            dataset_names <- names(stat_data)
            default_name <- dataset_names[1]
            dataset <- stat_data[[1]]
            
            shinyApp(
                ui = fluidPage(
                    selectInput("data", "Dataset:", selected=default_name, choices=dataset_names),
                    selectInput("type", "Type:", selected="adj.P.Val", choices = c("adj.P.Val", "P.Value")),
                    checkboxInput("scaleaxis", "Scale Y axes"),
                    sliderInput("thres", "Threshold:", value=0.1, min=0, max=1, step=0.01),
                    plotOutput("vulcs"),
                    plotOutput("mas")
                ),
                server = function(input, output) {
                    
                    output$vulcs = renderPlot({
                        
                        rdf <- data.frame(rowData(stat_data[[input$data]]))
                        plts <- list()
                        for (contrast in contrasts) {
                            
                            tbl <- data.frame(
                                P.Value=rdf[[paste(contrast, p_col, sep=".")]],
                                adj.P.Val=rdf[[paste(contrast, q_col, sep=".")]],
                                logFC=rdf[[paste(contrast, fold_col, sep=".")]],
                                AveExpr=rdf[[paste(contrast, expr_col, sep=".")]]
                            )
                            plt <- stv$vulc(tbl, sig_thres = input$thres, sig_col_name = input$type, na.rm=TRUE) +
                                theme_classic() +
                                ggtitle(contrast)
                            plts[[contrast]] <- plt
                        }
                        
                        if (input$scaleaxis) {
                            max_vals <- lapply(plts, function(plt) {
                                layer_scales(plt)$y$range$range[2]
                            })
                            plts <- lapply(plts, function(plt) {
                                plt <- plt + ylim(0, 1.01 * max(unlist(max_vals)))
                            })
                        }
                        
                        grid.arrange(grobs=plts, ncol=3)
                    }, height=300)
                    
                    output$mas = renderPlot({
                        
                        rdf <- data.frame(rowData(stat_data[[input$data]]))
                        plts <- list()
                        for (contrast in contrasts) {
                            
                            tbl <- data.frame(
                                P.Value=rdf[[paste(contrast, p_col, sep=".")]],
                                adj.P.Val=rdf[[paste(contrast, q_col, sep=".")]],
                                logFC=rdf[[paste(contrast, fold_col, sep=".")]],
                                AveExpr=rdf[[paste(contrast, expr_col, sep=".")]]
                            )
                            plt <- stv$ma(tbl, sig_thres = input$thres, sig_col_name = input$type, na.rm=TRUE) +
                                theme_classic() + 
                                ggtitle(contrast)
                            plts[[contrast]] <- plt
                        }
                        
                        if (input$scaleaxis) {
                            max_vals <- lapply(plts, function(plt) { layer_scales(plt)$y$range$range[2] })
                            min_vals <- lapply(plts, function(plt) { layer_scales(plt)$y$range$range[1] })
                            plts <- lapply(plts, function(plt) { 
                                plt <- plt + ylim(1.01 * min(unlist(min_vals)), 1.01 * max(unlist(max_vals))) 
                            })
                        }
                        
                        grid.arrange(grobs=plts, ncol=3)
                    }, height=300)
                },
                options=list(height=height)
            )
        },
        venn_widgets = function(name, stat_data, contrasts, p_col="P.Value", q_col="adj.P.Val", 
                                fold_col="logFC", height=1100) {
            
            p_cols <- paste(contrasts, p_col, sep=".")
            dataset_names <- names(stat_data)
            default_name <- dataset_names[1]
            dataset <- stat_data[[1]]
            
            shinyApp(
                ui = fluidPage(
                    selectInput("data", "Dataset:", selected=default_name, choices=dataset_names),
                    selectInput("type", "Type:", selected="adj.P.Val", choices = c("adj.P.Val", "P.Value")),
                    sliderInput("thres", "Threshold:", value=0.1, min=0, max=1, step=0.01),
                    sliderInput("fold", "Fold threshold:", value=0, min=0, max=5, step=0.25),
                    plotOutput("plots")
                ),
                server = function(input, output) {
                    output$plots = renderPlot({

                        rdf <- data.frame(rowData(stat_data[[input$data]]))
                        out <- stv$plot_comp_venns(
                            rdf, 
                            contrasts, 
                            base_sig_col = input$type,
                            base_fold_name = fold_col,
                            sig_thres = input$thres, 
                            log2_fold_thres = input$fold
                        )
                        grid.arrange(grobs=out, ncol=3)
                    })
                },
                options=list(height=800)
            )
        },
        
        table_widget = function(name, stat_data, height=900, default_selected=NULL) {
            
            if (is.null(default_selected)) {
                default_selected <- colnames(rowData(stat_data[[1]]))
            }
            dataset_names <- names(stat_data)
            default_name <- dataset_names[1]
            full_annotation <- rowData(stat_data[1])

            shinyApp(
                ui = fluidPage(
                    selectInput("data", "Dataset:", selected=default_name, choices=dataset_names),
                    selectInput("fields", "Shown fields", choices=colnames(full_annotation), multiple=TRUE, selected=default_selected),
                    splitLayout(
                        sliderInput("fdrthres", "FDR thres.", 0.1, min=0, max=1, step=0.01),
                        sliderInput("decimals", "Decimals", 2, min=0, max=10, step=1)
                    ),
                    checkboxInput("exclusive", "Only show significant in all groups"),
                    downloadButton('download', "Download Table"),
                    DT::dataTableOutput("table")
                ),
                server = function(input, output) {
                    
                    thedata <- reactive({
                        
                            unfiltered <- rowData(stat_data[[input$data]]) %>% 
                                data.frame() %>%
                                select(input$fields[order(match(input$fields, default_selected))])
                            
                            if (!input$exclusive) {
                                unfiltered %>% 
                                    filter(batch1.adj.P.Val < input$fdrthres | batch2.adj.P.Val < input$fdrthres | batch2.adj.P.Val < input$fdrthres)
                            }
                            else {
                                unfiltered %>%
                                    filter(batch1.adj.P.Val < input$fdrthres & batch2.adj.P.Val < input$fdrthres & batch3.adj.P.Val < input$fdrthres)
                            }
                        })
                    
                    output$table = DT::renderDataTable({
                        
                        thedata() %>% 
                            datatable(options=list(
                                pageLength=10, 
                                scrollX=TRUE, 
                                autoWidth=TRUE,
                                columnDefs=list(list(width="10px", targets="_all"))
                            )) %>%
                            DT::formatRound(columns=input$fields[input$fields %in% numeric_cols], digits=input$decimals)
                    })
                    
                    output$download <- downloadHandler(
                        filename = function() {"result_table.tsv"},
                        content = function(fname) {
                            write_tsv(thedata(), fname)
                        } 
                    )
                },
                options=list(height=900)
            )
            
        },
        spotcheck_widget = function(name, stat_data, id_col, split_col, split_vals, contrast_cond,
                                    height=1100, default_data=NULL, default_gene=NULL, color_cond=NULL) {
            
            if (is.null(default_data)) {
                default_data <- colnames(rowData(stat_data[[1]]))
            }
            
            if (is.null(default_gene)) {
                default_gene <- colnames(rowData(stat_data[[1]]))
            }
            
            if (is.null(color_cond)) {
                color_cond <- colnames(colData(dataset))[1]
            }
            
            dataset_names <- names(stat_data)
            default_name <- dataset_names[1]
            dataset <- stat_data[[1]]
            
            # selectInput("cond", "Condition:", choices = colnames(colData(dataset)), selected=default_cond)
            
            shinyApp(
                ui = fluidPage(
                    selectInput("data", "Dataset:", selected=default_name, choices=dataset_names),
                    selectInput("rowid", "Row ID", selected=default_gene, choices = row_ids),
                    selectInput("color", "Coloring category", selected=color_cond, choices = colnames(colData(dataset))),
                    checkboxInput("colorscatter", "Color boxplot scatter on cond"),
                    plotOutput("scatters"),
                    plotOutput("contrast")
                ),
                server = function(input, output) {
                    output$scatters = renderPlot({
                        
                        target <- comb_se[rowData(stat_data[[input$data]])[[id_col]] == input$rowid, ]
                        make_scatter <- function(row_ses, title) {
                            fert_vals <- colData(row_ses)$Fertility
                            expr_vals <- assay(row_ses)[1, ]
                            color <- colData(row_ses)[[input$color]]
                            ggplot(data.frame(expr=assay(row_ses)[1,], fert=fert_vals, color=color), aes(expr, fert, color=color)) + 
                                geom_point() +
                                xlab("Expression") +
                                ylab("Fertility") +
                                theme_classic() +
                                ggtitle(title)
                        }
                        
                        b1_plt <- make_scatter(target[, colData(target)[[split_col]] == split_vals[1]], "Batch 1")
                        b2_plt <- make_scatter(target[, colData(target)[[split_col]] == split_vals[2]], "Batch 2")
                        b3_plt <- make_scatter(target[, colData(target)[[split_col]] == split_vals[3]], "Batch 3")
                        
                        ggarrange(b1_plt, b2_plt, b3_plt, ncol=3, common.legend=TRUE, legend="bottom")
                    })
                    
                    output$contrast = renderPlot({
                        
                        target <- comb_se[rowData(stat_data[[input$data]])[[id_col]] == input$rowid, ]
                        
                        make_box <- function(row_ses, title) {
                            expr_vals <- assay(row_ses)[1, ]
                            high_fert <- colData(row_ses[1, ])[[contrast_cond]]
                            stat_df <- data.frame(expr=expr_vals, high_fert=high_fert, color=colData(row_ses[1, ])[[input$color]])
                            plt <- ggplot(stat_df, aes(x=high_fert, y=expr_vals, fill=high_fert)) + geom_boxplot(alpha=0.5) + theme_classic()
                            
                            if (input$colorscatter) plt + geom_point(aes(color=color))
                            else plt + geom_point()
                        }
                        
                        b1_plt <- make_box(target[, colData(target)[[split_col]] == split_vals[1]], split_vals[1])
                        b2_plt <- make_box(target[, colData(target)[[split_col]] == split_vals[2]], split_vals[2])
                        b3_plt <- make_box(target[, colData(target)[[split_col]] == split_vals[3]], split_vals[3])
                        
                        ggarrange(b1_plt, b2_plt, b3_plt, ncol=3, common.legend=TRUE, legend="bottom")
                    })
                },
                options=list(height=height)
            )
            
            
        }
    ),
    private = list(
        parse_dataset = function(dataset, outliers) {
            non_outliers <- colnames(dataset)[!colnames(dataset) %in% outliers]
            sdf <- assay(dataset)[, non_outliers]
            ddf <- data.frame(colData(dataset)) %>% filter(sample %in% non_outliers)
            rdf <- rowData(dataset) %>% data.frame()
            list("sdf"=sdf, "ddf"=ddf, "rdf"=rdf)
        },
        make_scale = 1.01
    )
) 

mywidgets <- MyWidgets$new()
mw <- mywidgets
print("Loading module to 'mywidgets' and 'mw'")
