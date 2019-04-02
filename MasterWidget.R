library(R6)
library(shiny)
library(tidyverse)
library(gridExtra)
library(DT)
library(ggpubr)

# Datasets expected format: Named list linked to SummarizedExperiment instances
# Does it make sense to weave everything into a single widget?

pf <- MasterWidgetPlotFuncs$new()

MasterWidget <- R6Class(
    public = list(
        general_qual_widget = function(datasets, outlier_sets=NULL, height=800, default_plot="bar", default_cond=NULL, 
                                       interactive=TRUE, contrast_suffix=NULL) {
            
            dataset_names <- names(datasets)
            default_name <- dataset_names[1]
            dataset <- datasets[[1]]
            
            if (is.null(default_cond)) {
                default_cond <- colnames(colData(dataset))[1]
            }
            
            if (!interactive) {
                plt <- self$make_sample_dist(datasets, default_name, default_cond, outliers=outlier_sets, usecount=NULL)
                return(plt)
            }
            
            show_cols=c("P.Value", "adj.P.Val", "logFC", "AveExpr")
            
            all_settings <- c(
                "tabs",
                "data1", "data2", "checkgroup",
                "cond1", "cond2",
                "show_na", "show_mean",
                "subset", "fulldata",
                "output_path", "output_width", "output_height", "output_dpi"
            )
            
            explore_plots <- c("Barplot", "QQ", "Density", "PCA", "Cluster")
            
            height_step_size <- 50
            
            # ("Barplot", numericInput(inputId="Barplot_height", "Plot height", value=500, step=height_step_size)),
            # tabPanel("QQ", numericInput(inputId="QQ_height", "Plot height", value=500, step=height_step_size)),
            # tabPanel("Density", numericInput(inputId="Density_height", "Plot height", value=500, step=height_step_size)),
            # tabPanel("PCA", numericInput(inputId="PCA_height", "Plot height", value=500, step=height_step_size)),
            # tabPanel("Cluster", numericInput(inputId="Cluster_height", "Plot height", value=500, step=height_step_size)),
            # tabPanel("Hists", numericInput(inputId="Hists_height", "Plot height", value=500, step=height_step_size)),
            # tabPanel("Venns",
            
            shinyApp(
                ui = fluidPage(
                    
                    tags$head(tags$style(HTML(".shiny-split-layout > div { overflow: visible; }"))),
                    
                    titlePanel("Explomics"),
                    
                    # id = "tabs",
                    # type = "tabs",
                    # tabPanel("Barplot", 
                    
                    sidebarLayout(
                        sidebarPanel(
                            
                            tabsetPanel(
                                id = "options_tabs",
                                type = "tabs",
                                tabPanel(
                                    "Settings",
                                    # General
                                    conditionalPanel(
                                        condition = "input.tabs == 'Barplot' || input.tabs == 'QQ' || input.tabs == 'Density' || input.tabs == 'PCA' || input.tabs == 'Cluster'",
                                        selectInput("data1", "Dataset 1:", selected=default_name, choices=dataset_names),
                                        selectInput("data2", "Dataset 2:", selected=default_name, choices=dataset_names),
                                        checkboxGroupInput("checkgroup", "Remove group", choices=names(outlier_sets), selected=NULL),
                                        selectInput("cond1", "Condition:", choices = colnames(colData(dataset)), selected=default_cond),
                                        selectInput("cond2", "Condition:", choices = colnames(colData(dataset)), selected=default_cond)
                                    ),
                                    conditionalPanel(
                                        condition = "input.tabs == 'Venn' || input.tabs == 'Hists'",
                                        selectInput("data1", "Dataset 1:", selected=default_name, choices=dataset_names)
                                    ),
                                    
                                    conditionalPanel(
                                        condition = "input.tabs == 'Barplot'",
                                        checkboxInput("show_na", "Show NA counts", value=FALSE),
                                        checkboxInput("show_mean", "Show mean", value=FALSE)
                                    ),
                                    
                                    conditionalPanel(
                                        condition = "input.tabs == 'QQ' || input.tabs == 'Density'",
                                        numericInput("subset", "Partial data:", value=500, min=100, max=nrow(assay(dataset))),
                                        checkboxInput("fulldata", "Use full data", value = FALSE)
                                    ),
                                    
                                    # Principal component
                                    conditionalPanel(
                                        condition = "input.tabs == 'PCA'",
                                        checkboxInput("as_label", "Show as text"),
                                        selectInput("text_labels", "Text labels:", selected=default_cond, choices=colnames(colData(dataset))),
                                        numericInput("pc_comps", "Max PCs in Scree", value=6, min=1),
                                        selectInput("pc1_plt1", "PC1 (plot1):", choices=1:8, selected=1),
                                        selectInput("pc1_plt2", "PC1 (plot2):", choices=1:8, selected=1),
                                        selectInput("pc2_plt1", "PC2 (plot1):", choices=1:8, selected=2),
                                        selectInput("pc2_plt2", "PC2 (plot2):", choices=1:8, selected=2)
                                    ),
                                    
                                    # Histograms
                                    conditionalPanel(
                                        condition = "input.tabs == 'Hists'",
                                        numericInput("hist_bins", "Bin count", value=50, min=1, max=100),
                                        selectInput("target_col", "Target columns", selected=show_cols[1], choices=show_cols),
                                        checkboxInput("scaleaxis", "Scale Y axes"),
                                        checkboxInput("limitxaxis", "Limit X axes")
                                    ),
                                    
                                    # Venns
                                    conditionalPanel(
                                        condition = "input.tabs == 'Venns'",
                                        selectInput("type", "Type:", selected="adj.P.Val", choices = c("adj.P.Val", "P.Value")),
                                        sliderInput("thres", "Threshold:", value=0.1, min=0, max=1, step=0.01),
                                        sliderInput("fold", "Fold threshold:", value=0, min=0, max=5, step=0.25),
                                        textInput("fold_base", "Fold base:", value="logFC")
                                    )
                                 ),
                                tabPanel(
                                    "Download",
                                    fluidRow(
                                        textInput(inputId="output_path", label="Path"),
                                        numericInput(inputId="output_width", label="Width", value=2000),
                                        numericInput(inputId="output_height", label="Height", value=1000),
                                        numericInput(inputId="output_dpi", label="DPI", value=150),
                                        downloadButton(outputId="download", label="Download plot"),
                                        downloadButton(outputId="download_params", label="Download parameters")
                                    )
                                )
                            )
                        ),
                        mainPanel(
                            # Technical
                            tabsetPanel(
                                id = "tabs",
                                type = "tabs",
                                tabPanel("Barplot", numericInput(inputId="Barplot_height", "Plot height", value=500, step=height_step_size)),
                                tabPanel("QQ", numericInput(inputId="QQ_height", "Plot height", value=500, step=height_step_size)),
                                tabPanel("Density", numericInput(inputId="Density_height", "Plot height", value=500, step=height_step_size)),
                                tabPanel("PCA", numericInput(inputId="PCA_height", "Plot height", value=500, step=height_step_size)),
                                tabPanel("Cluster", numericInput(inputId="Cluster_height", "Plot height", value=500, step=height_step_size)),
                                tabPanel("Hists", numericInput(inputId="Hists_height", "Plot height", value=500, step=height_step_size)),
                                tabPanel("Venns", numericInput(inputId="Venns_height", "Plot height", value=500, step=height_step_size))
                            ),
                            uiOutput("BarplotUI")
                        )
                    )
                ),
                server = function(session, input, output) {
                    
                    if (is.null(contrast_suffix)) {
                        hideTab(inputId="tabs", target="Hists")
                    }
                    
                    plotHeight <- reactive({
                        target_name <- paste0(input$tabs, "_height")
                        input[[target_name]]
                    })
                    
                    output$BarplotUI <- renderUI({
                        plotOutput(input$tabs, height=plotHeight())
                    })

                    output$Barplot = renderPlot({
                        pf$do_bar(datasets, input)
                    })
                    
                    output$QQ = renderPlot({
                        pf$do_qq(datasets, input)
                    })
                    
                    output$Density = renderPlot({
                        pf$do_density(datasets, input)
                    })
                    
                    output$PCA = renderPlot({
                        pf$do_pca(datasets, input)
                    })
                    
                    output$Cluster = renderPlot({
                        pf$do_cluster(datasets, input)
                    })
                    
                    output$Hists = renderPlot({
                        pf$do_hists(datasets, input, contrast_suffix)
                    })
                    
                    output$Venns = renderPlot({
                        pf$do_venns(datasets, input, contrast_suffix)
                    })
                    
                    output$download <- downloadHandler(
                        
                        filename = function() {
                            input$output_path
                        },
                        content = function(file) {
                            
                            if (input$tabs == "Barplot") {
                                curr_plot <- pf$do_bar
                            }
                            else if (input$tabs == "QQ") {
                                curr_plot <- pf$do_qq
                            }
                            else if (input$tabs == "Density") {
                                curr_plot <- pf$do_density
                            }
                            else if (input$tabs == "PCA") {
                                curr_plot <- pf$do_pca
                            }
                            else if (input$tabs == "Cluster") {
                                curr_plot <- pf$do_cluster
                            }
                            else if (input$tabs == "Hists") {
                                curr_plot <- pf$do_hists
                            }
                            else if (input$tabs == "Venns") {
                                curr_plot <- pf$do_venns
                            }
                            else {
                                stop("Unknown option: ", input$tabs)
                            }
                            
                            plt <- curr_plot(datasets, input)
                            
                            inch_width <- input$output_width / input$output_dpi
                            inch_height <- input$output_height / input$output_dpi
                            ggsave(file, plot=plt, device="png", width=inch_width, height=inch_height, units="in", dpi=input$output_dpi)
                        }
                    )
                    
                    output$download_params <- downloadHandler(

                        filename = function() {
                            paste(input$output_path, "params", sep=".")
                        },
                        content = function(file) {
                            
                            settings <- list()
                            for (setting in all_settings) {
                                settings[[setting]] <- c(setting, input[[setting]])
                            }
                            
                            setting_table <- data.frame(do.call("rbind", settings))
                            colnames(setting_table) <- c("setting", "value")
                            write_tsv(setting_table, path=file)
                        }
                    )

                    observe({
                        self$update_input_choices(session, datasets, "cond1", input$data1, input$cond1)
                        self$update_input_choices(session, datasets, "cond2", input$data2, input$cond2)
                        self$update_input_choices(session, datasets, "text_labels", input$data1, input$text_labels)
                    })
                },
                options=list(height=height)
            )
        },
        
        update_input_choices = function(session, datasets, target, data_name, cond_name) {
            
            choices <- colnames(colData(datasets[[data_name]]))
            if (cond_name %in% choices) selected <- cond_name
            else selected <- default_cond
            
            updateSelectInput(
                session,
                target,
                choices=choices,
                selected=selected
            )
        },
        
        scatter_widgets = function(stat_data, contrast_suffix, p_col="P.Value", q_col="adj.P.Val", 
                                   fold_col="logFC", expr_col="AveExpr", height=1100) {
            
            dataset_names <- names(stat_data)
            default_name <- dataset_names[1]

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
                        contrasts <- private$get_contrasts_from_suffix(stat_data[[input$data]], contrast_suffix)
                        
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
                        
                        grid.arrange(grobs=plts, ncol=3, top=paste0("Dataset: ", input$data))
                    }, height=300)
                    
                    output$mas = renderPlot({
                        
                        rdf <- data.frame(rowData(stat_data[[input$data]]))
                        plts <- list()
                        contrasts <- private$get_contrasts_from_suffix(stat_data[[input$data]], contrast_suffix)
                        
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
                        
                        grid.arrange(grobs=plts, ncol=3, top=paste0("Dataset: ", input$data))
                    }, height=300)
                },
                options=list(height=height)
            )
        },
        # venn_widgets = function(stat_data, contrast_suffix, p_col="P.Value", q_col="adj.P.Val", 
        #                         fold_col="logFC", height=1100) {
        #     
        #     dataset_names <- names(stat_data)
        #     default_name <- dataset_names[1]
        # 
        #     shinyApp(
        #         ui = fluidPage(
        #             selectInput("data", "Dataset:", selected=default_name, choices=dataset_names),
        #             selectInput("type", "Type:", selected="adj.P.Val", choices = c("adj.P.Val", "P.Value")),
        #             sliderInput("thres", "Threshold:", value=0.1, min=0, max=1, step=0.01),
        #             sliderInput("fold", "Fold threshold:", value=0, min=0, max=5, step=0.25),
        #             plotOutput("plots")
        #         ),
        #         server = function(input, output) {
        #             output$plots = renderPlot({
        # 
        #                 contrasts <- private$get_contrasts_from_suffix(stat_data[[input$data]], contrast_suffix)
        #                 
        #                 rdf <- data.frame(rowData(stat_data[[input$data]]))
        #                 out <- stv$plot_comp_venns(
        #                     rdf, 
        #                     contrasts, 
        #                     base_sig_col = input$type,
        #                     base_fold_name = fold_col,
        #                     sig_thres = input$thres, 
        #                     log2_fold_thres = input$fold
        #                 )
        #                 grid.arrange(grobs=out, ncol=3, top=paste0("Dataset: ", input$data))
        #             })
        #         },
        #         options=list(height=800)
        #     )
        # },
        
        table_widget = function(stat_data, height=1200, default_selected=NULL, default_filter=NULL) {
            
            if (is.null(default_selected)) {
                default_selected <- colnames(rowData(stat_data[[1]]))
            }
            dataset_names <- names(stat_data)
            default_name <- dataset_names[1]
            full_annotation <- rowData(stat_data[[1]])

            shinyApp(
                ui = fluidPage(
                    selectInput("data", "Dataset:", selected=default_name, choices=dataset_names),
                    selectInput("fields", "Shown fields", choices=colnames(full_annotation), multiple=TRUE, selected=default_selected),
                    selectInput("filters", "Filter fields", choices=colnames(full_annotation), multiple=TRUE, selected=default_filter),
                    splitLayout(
                        sliderInput("filterthres", "Filter thres.", 0.1, min=0, max=1, step=0.01),
                        sliderInput("decimals", "Decimals", 2, min=0, max=10, step=1)
                    ),
                    checkboxInput("exclusive", "Only show significant in all groups"),
                    downloadButton('download', "Download Table"),
                    DT::dataTableOutput("table")
                ),
                server = function(input, output, session) {
                    
                    thedata <- reactive({
                        
                        retained <- rowData(stat_data[[input$data]]) %>% data.frame()
                        if (length(input$filters) > 0) {
                            
                            if (input$exclusive) {
                                unique_retained <- retained
                                
                                # Only include features passing all filters
                                for (filter in input$filters) {
                                    unique_retained <- unique_retained %>% filter(UQ(as.name(filter)) < input$filterthres)
                                }
                                retained <- unique_retained
                            }
                            else {
                                all_retained <- NULL
                                for (filter in input$filters) {
                                    # Include features passing at least one filter
                                    filter_retained <- retained %>% filter(UQ(as.name(filter)) < input$filterthres)
                                    all_retained <- rbind(all_retained, filter_retained)
                                }
                                retained <- all_retained %>% distinct()
                            }
                        }
                        
                        filtered_selected <- retained %>% 
                            dplyr::select(input$fields) %>%
                            data.frame()

                        filtered_selected
                    })
                    
                    observe({
                        new_choices <- colnames(rowData(stat_data[[input$data]]))
                        updateSelectInput(
                            session,
                            "fields",
                            choices=new_choices,
                            selected=input$fields
                        )
                    })
                    
                    output$table = DT::renderDataTable({
                        
                        thedata() %>% 
                            datatable(options=list(
                                pageLength=10, 
                                scrollX=TRUE,
                                # autoWidth=TRUE,
                                columnDefs=list(list(width="10px", targets="_all"))
                            )) %>%
                            DT::formatRound(columns=input$fields, digits=input$decimals)
                    })
                    
                    output$download <- downloadHandler(
                        filename = function() {"result_table.tsv"},
                        content = function(fname) {
                            write_tsv(thedata(), fname)
                        } 
                    )
                },
                options=list(height=height)
            )
            
        },
        spotcheck_widget = function(stat_data, id_col, contrast_suffix, contrast_cond, split_col=NULL,
                                    base_height=600, default_data=NULL, default_gene=NULL, color_cond=NULL, corr_col=NULL,
                                    outlier_sets=NULL, default_label=NULL) {

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
            dataset <- stat_data[[1]]
            row_ids <- rowData(dataset)[[id_col]]

            if (!is.null(corr_col)) {
                height <- base_height + 500
            }
            else {
                height <- base_height + 500
            }
            
            shinyApp(
                ui = fluidPage(
                    tags$head(tags$style(HTML(".shiny-split-layout > div { overflow: visible; }"))),
                    splitLayout(
                        selectInput("data", "Dataset:", selected=default_data, choices=dataset_names),
                        selectInput("color", "Coloring category", selected=color_cond, choices = colnames(colData(dataset)))
                    ),
                    splitLayout(
                        selectInput("rowid", "Row ID", selected=default_gene, choices = row_ids),
                        selectInput("label_type", "Label category", selected=default_label, choices=colnames(colData(dataset)))
                    ),
                    splitLayout(
                        fluidPage(
                            checkboxInput("colorscatter", "Color boxplot scatter on cond"),
                            checkboxInput("showlabels", "Show labels instead of scatter dots")
                        ),
                        checkboxGroupInput("outliers", "Remove group", choices=names(outlier_sets), selected=NULL)
                    ),
                    selectInput("contrast", "Contrast column:", selected=contrast_cond, choices=colnames(colData(dataset))),
                    plotOutput("contrast"),
                    plotOutput("scatters")
                ),
                server = function(input, output) {
                    
                    output$contrast = renderPlot({
                        
                        contrast_cond <- input$contrast
                        
                        se <- stat_data[[input$data]]
                        target <- se[rowData(stat_data[[input$data]])[[id_col]] == input$rowid, ]
                        
                        outliers <- unname(unlist(outlier_sets[input$outliers]))
                        non_outliers <- colnames(target)[!colnames(target) %in% outliers]
                        target <- target[, non_outliers]
                        
                        if (!is.null(split_col)) {
                            split_col <- colData(target)[[split_col]]
                            split_vals <- sort(unique(split_col))
                            plts <- lapply(split_vals, function(split_val, data, split_col, contrast_cond, label_type, color, show_labels, color_scatter) { 
                                data_part <- data[, split_col == split_val] 
                                name <- paste("Split:", split_val)
                                plt <- private$make_box(data_part, name, contrast_cond, label_type, color, show_labels, color_scatter)
                            }, 
                            data=target, split_col=split_col, contrast_cond=contrast_cond, label_type=input$label_type, color=input$color, show_labels=input$showlabels, color_scatter=input$colorscatter)
                        }
                        else {
                            plts <- list()
                            plts[[1]] <- private$make_box(target, "Contrast", contrast_cond, input$label_type, input$color, input$showlabels, input$colorscatter)
                        }
                        
                        fig <- ggpubr::ggarrange(plotlist=plts, ncol=3, common.legend=TRUE, legend="bottom")
                        ggpubr::annotate_figure(fig, top=text_grob(paste0("Dataset: ", input$data)))
                    })
                    
                    output$scatters = renderPlot({
                        
                        if (is.null(corr_col)) {
                            return()
                        }
                        
                        contrast_cond <- input$contrast
                        
                        se <- stat_data[[input$data]]
                        target <- se[rowData(stat_data[[input$data]])[[id_col]] == input$rowid, ]
                        
                        outliers <- unname(unlist(outlier_sets[input$outliers]))
                        non_outliers <- colnames(target)[!colnames(target) %in% outliers]
                        target <- target[, non_outliers]
                        
                        if (!is.null(split_col)) {
                            split_col <- colData(target)[[split_col]]
                            split_vals <- sort(unique(split_col))
                            plts <- lapply(split_vals, function(split_val, data, split_col, contrast_cond, label_type, color, show_labels, color_scatter, corr_col) { 
                                data_part <- data[, split_col == split_val] 
                                name <- paste("Split:", split_val)
                                plt <- private$make_scatter(data_part, name, show_labels, color, label_type, corr_col)
                            }, 
                            data=target, split_col=split_col, contrast_cond=contrast_cond, label_type=input$label_type, color=input$color, show_labels=input$showlabels, corr_col=corr_col)
                        }
                        else {
                            plts <- list()
                            plts[[1]] <- private$make_scatter(target, "Scatter", input$showlabels, input$color, input$label_type, corr_col)
                        }
                        
                        fig <- ggpubr::ggarrange(plotlist=plts, ncol=3, common.legend=TRUE, legend="bottom")
                        ggpubr::annotate_figure(fig, top=text_grob(paste0("Dataset: ", input$data)))
                    })
                },
                options=list(height=height)
            )
        }
    ),
    private = list(
        make_scale = 1.01,
        make_scatter = function(row_ses, title, show_labels, color, label_type, corr_col) {
            
            fert_vals <- colData(row_ses)[[corr_col]]
            color <- colData(row_ses)[[color]]
            label_names <- colData(row_ses)[, label_type]
            
            plt <- ggplot(data.frame(expr=assay(row_ses)[1,], fert=fert_vals, color=color, label=label_names), 
                          aes(expr, fert, color=color, label=label)) + 
                xlab("Expression") +
                ylab("Fertility") +
                theme_classic() +
                ggtitle(title)
            
            if (show_labels) target_geom <- geom_text
            else target_geom <- geom_point
            
            plt + target_geom(na.rm=TRUE)
        },
        make_box = function(row_ses, title, contrast_cond, label_type, color, showlabels, colorscatter) {
            expr_vals <- assay(row_ses)[1, ]
            high_fert <- colData(row_ses[1, ])[[contrast_cond]]
            label_names <- colData(row_ses)[, label_type]

            stat_df <- data.frame(
                expr=expr_vals, 
                high_fert=high_fert, 
                color=colData(row_ses[1, ])[[color]],
                labels=label_names
            )
            
            plt <- ggplot(
                stat_df, 
                aes(x=high_fert, y=expr_vals, fill=high_fert, label=labels)
            ) + 
                geom_boxplot(alpha=0.5, na.rm=TRUE) + 
                theme_classic()
            
            if (showlabels) target_geom <- geom_text
            else target_geom <- geom_point
            
            if (colorscatter) plt + target_geom(aes(color=color))
            else plt + target_geom(na.rm=TRUE)
        }
    )
) 

masterwidget <- MasterWidget$new()
print("Loading module to 'masterwidget'")
