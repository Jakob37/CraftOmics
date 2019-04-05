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
        general_qual_widget = function(datasets, outlier_sets=NULL, height=800, default_cond=NULL, 
                                       interactive=TRUE, contrast_suffix=NULL) {
            
            dataset_names <- names(datasets)
            default_name <- dataset_names[1]
            dataset <- datasets[[1]]
            annot_col_names <- colnames(rowData(dataset))
            if (!is.null(contrast_suffix)) {
                contrasts <- private$get_contrasts_from_suffix(annot_col_names, contrast_suffix)
                contrast_suffixes <- private$get_suffixes_from_contrast(annot_col_names, contrasts[1])
            }

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
            
            shinyApp(
                ui = fluidPage(
                    
                    tags$head(tags$style(HTML(".shiny-split-layout > div { overflow: visible; }"))),
                    
                    titlePanel("Ponderomics"),
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
                                        selectInput("cond1", "Condition 1:", choices = colnames(colData(dataset)), selected=default_cond),
                                        selectInput("cond2", "Condition 2:", choices = colnames(colData(dataset)), selected=default_cond)
                                    ),
                                    
                                    conditionalPanel(
                                        condition = "input.tabs == 'Venns' || input.tabs == 'Hists' || input.tabs == 'Scatter' || input.tabs == 'Table'",
                                        selectInput("stat_data", "Dataset:", selected=default_name, choices=dataset_names, size=20, selectize=FALSE)
                                    ),
                                    
                                    # Histograms
                                    conditionalPanel(
                                        condition = "input.tabs == 'Hists'",
                                        selectInput("hist_target", "Target col:", selected=contrast_suffix, choices = contrast_suffixes),
                                        numericInput("hist_bins", "Bin count", value=100, min=1),
                                        checkboxInput("scale_x_axis", "Share X axis"),
                                        checkboxInput("scale_y_axis", "Share Y axis")
                                    ),
                                    
                                    # Venns
                                    conditionalPanel(
                                        condition = "input.tabs == 'Venns'",
                                        selectInput("venn_type", "Type:", selected=contrast_suffix, choices = contrast_suffixes),
                                        selectInput("venn_fold", "Fold col:", selected="logFC", choices = contrast_suffixes),
                                        sliderInput("venn_thres", "Threshold:", value=0.1, min=0, max=1, step=0.01),
                                        checkboxInput("venn_inverse", "Check if greater than", value=FALSE),
                                        checkboxInput("threeway_venn", "Three-way Venn", value=FALSE)
                                    ),
                                    
                                    # General scatter
                                    conditionalPanel(
                                        condition = "input.tabs == 'Scatter'",
                                        selectInput("scatter_x", "X-axis:", selected=contrast_suffixes[2], choices=contrast_suffixes),
                                        selectInput("scatter_y", "Y-axis:", selected=contrast_suffixes[1], choices=contrast_suffixes),
                                        selectInput("scatter_color", "Color cond:", selected="adj.P.Val", choices=contrast_suffixes),
                                        numericInput("scatter_color_cutoff", "Color cutoff:", value=0.05, min=0, max=1),
                                        checkboxInput("scatter_minuslog_y", "Minus-log y:", value=FALSE)
                                    ),
                                    
                                    # Table
                                    conditionalPanel(
                                        condition = "input.tabs == 'Table'",
                                        #             selectInput("data", "Dataset:", selected=default_name, choices=dataset_names),
                                        #             selectInput("fields", "Shown fields", choices=colnames(full_annotation), multiple=TRUE, selected=default_selected),
                                        #             selectInput("filters", "Filter fields", choices=colnames(full_annotation), multiple=TRUE, selected=default_filter),
                                        #             splitLayout(
                                        #                 sliderInput("filterthres", "Filter thres.", 0.1, min=0, max=1, step=0.01),
                                        #                 sliderInput("decimals", "Decimals", 2, min=0, max=10, step=1)
                                        #             ),
                                        #             checkboxInput("exclusive", "Only show significant in all groups"),
                                        #             downloadButton('download', "Download Table"),
                                        #             DT::dataTableOutput("table")
                                        selectInput("fields", "Shown fields", choices=annot_col_names, multiple=TRUE, selected=annot_col_names[1:5])
                                    )
                                 ),
                                tabPanel(
                                    "Customize",
                                    # Principal component
                                    conditionalPanel(
                                        condition = "input.tabs == 'PCA'",
                                        checkboxInput("as_label", "Show as text"),
                                        selectInput("text_labels", "Text labels:", selected=default_cond, choices=colnames(colData(dataset))),
                                        numericInput("pc_comps", "Max PCs in Scree", value=6, min=1),
                                        selectInput("pc1_plt1", "PC1 (plot1):", choices=1:8, selected=1),
                                        selectInput("pc2_plt1", "PC2 (plot1):", choices=1:8, selected=2),
                                        selectInput("pc1_plt2", "PC1 (plot2):", choices=1:8, selected=1),
                                        selectInput("pc2_plt2", "PC2 (plot2):", choices=1:8, selected=2),
                                        checkboxInput("pca_hide_loadings", "Hide loadings", value=FALSE)
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
                                    )
                                ),
                                
                                tabPanel(
                                    "Display",
                                    fluidRow(
                                        textInput(inputId="plot_title", label="Title", value=""),
                                        textInput(inputId="plot_subtext", label="Subtext", value=""),
                                        numericInput(inputId="plot_cols", label="Columns", value=2, step=1, min=1),
                                        textInput(inputId="legend_color", label="Leg. lab. Color", value=""),
                                        textInput(inputId="legend_fill", label="Leg. lab. Fill", value=""),
                                        selectInput(inputId="legend_pos", label="Legend: Pos", selected="bottom", 
                                                    choices = c("top", "bottom", "left", "right")),
                                        checkboxInput(inputId="legend_common", label="Legend: Common", value=TRUE),
                                        numericInput(inputId="title_size", label="Title size", value=14, step=1, min=4),
                                        numericInput(inputId="subtitle_size", label="Subtitle size", value=14, step=1, min=4),
                                        numericInput(inputId="axis_size", label="Label size", value=12, step=1, min=4),
                                        numericInput(inputId="ticks_size", label="Tick label size", value=10, step=1, min=4),
                                        numericInput(inputId="subtext_size", label="Subtext size", value=10, step=1, min=4),
                                        numericInput(inputId="subtext_wrap", label="Subtext wrap length", value=50, step=5),
                                        numericInput(inputId="sidebar_width", label="Sidebar width", value=4)
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
                            ), width=4
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
                                tabPanel("Venns", numericInput(inputId="Venns_height", "Plot height", value=500, step=height_step_size)),
                                tabPanel("Scatter", numericInput(inputId="Scatter_height", "Plot height", value=500, step=height_step_size)),
                                tabPanel("Spotcheck", numericInput(inputId="Spotcheck_height", "Plot height", value=500, step=height_step_size)),
                                tabPanel("Table", numericInput(inputId="Table_height", "Plot height", value=500, step=height_step_size)),
                                tabPanel("Profile", numericInput(inputId="Profile_height", "Plot height", value=500, step=height_step_size))
                            ),
                            uiOutput("PlotUI")
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
                    
                    output$PlotUI <- renderUI({
                        if (input$tabs == "Table") {
                            dataTableOutput("Table")
                        }
                        else {
                            plotOutput(input$tabs, height=plotHeight())
                        }
                        
                    })

                    output$Barplot = renderPlot({
                        plts <- pf$do_bar(datasets, input, outlier_sets)
                        pf$annotate(plts, input)
                    })
                    
                    output$QQ = renderPlot({
                        plt <- pf$do_qq(datasets, input, outlier_sets)
                        pf$annotate(plt, input)
                    })
                    
                    output$Density = renderPlot({
                        plt <- pf$do_density(datasets, input, outlier_sets)
                        pf$annotate(plt, input)
                    })
                    
                    output$PCA = renderPlot({
                        plt <- pf$do_pca(datasets, input, outlier_sets)
                        pf$annotate(plt, input)
                    })
                    
                    output$Cluster = renderPlot({
                        plt <- pf$do_cluster(datasets, input, outlier_sets)
                        pf$annotate(plt, input)
                    })
                    
                    output$Hists = renderPlot({
                        plt <- pf$do_hists(datasets, input, contrast_suffix, outlier_sets)
                        pf$annotate(plt, input)
                    })
                    
                    output$Venns = renderPlot({
                        plt <- pf$do_venns(datasets, input, contrast_suffix, outlier_sets)
                        pf$annotate(plt, input)
                    })
                    
                    output$Scatter = renderPlot({
                        plts <- pf$do_general_scatter(datasets, input, contrast_suffix, outlier_sets)
                        pf$annotate(plts, input)
                    })
                    
                    output$Spotcheck = renderPlot({
                        pf$do_spotcheck(datasets, input)
                    })
                    
                    output$Table = renderDT({
                        # iris
                        pf$do_table(datasets, input, outlier_sets)
                    }, options = list(lengthChange = TRUE))
                    
                    output$Profile = renderPlot({
                        pf$do_profile(datasets, input)
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
                        private$update_input_choices(session, datasets, "cond1", input$data1, input$cond1)
                        private$update_input_choices(session, datasets, "cond2", input$data2, input$cond2)
                        private$update_input_choices(session, datasets, "text_labels", input$data1, input$text_labels)
                        private$update_input_range(session, datasets, "venn_thres", input$data1, paste(contrasts, input$venn_type, sep="."))
                        
                        # For dynamic UI, maybe check:
                        # https://shiny.rstudio.com/articles/dynamic-ui.html
                    })
                },
                options=list(height=height)
            )
        },
        


        # table_widget = function(stat_data, height=1200, default_selected=NULL, default_filter=NULL) {
        #     
        #     if (is.null(default_selected)) {
        #         default_selected <- colnames(rowData(stat_data[[1]]))
        #     }
        #     dataset_names <- names(stat_data)
        #     default_name <- dataset_names[1]
        #     full_annotation <- rowData(stat_data[[1]])
        # 
        #     shinyApp(
        #         ui = fluidPage(
        #             selectInput("data", "Dataset:", selected=default_name, choices=dataset_names),
        #             selectInput("fields", "Shown fields", choices=colnames(full_annotation), multiple=TRUE, selected=default_selected),
        #             selectInput("filters", "Filter fields", choices=colnames(full_annotation), multiple=TRUE, selected=default_filter),
        #             splitLayout(
        #                 sliderInput("filterthres", "Filter thres.", 0.1, min=0, max=1, step=0.01),
        #                 sliderInput("decimals", "Decimals", 2, min=0, max=10, step=1)
        #             ),
        #             checkboxInput("exclusive", "Only show significant in all groups"),
        #             downloadButton('download', "Download Table"),
        #             DT::dataTableOutput("table")
        #         ),
        #         server = function(input, output, session) {
        #             
        #             thedata <- reactive({
        #                 
        #                 retained <- rowData(stat_data[[input$data]]) %>% data.frame()
        #                 if (length(input$filters) > 0) {
        #                     
        #                     if (input$exclusive) {
        #                         unique_retained <- retained
        #                         
        #                         # Only include features passing all filters
        #                         for (filter in input$filters) {
        #                             unique_retained <- unique_retained %>% filter(UQ(as.name(filter)) < input$filterthres)
        #                         }
        #                         retained <- unique_retained
        #                     }
        #                     else {
        #                         all_retained <- NULL
        #                         for (filter in input$filters) {
        #                             # Include features passing at least one filter
        #                             filter_retained <- retained %>% filter(UQ(as.name(filter)) < input$filterthres)
        #                             all_retained <- rbind(all_retained, filter_retained)
        #                         }
        #                         retained <- all_retained %>% distinct()
        #                     }
        #                 }
        #                 
        #                 filtered_selected <- retained %>% 
        #                     dplyr::select(input$fields) %>%
        #                     data.frame()
        # 
        #                 filtered_selected
        #             })
        #             
        #             observe({
        #                 new_choices <- colnames(rowData(stat_data[[input$data]]))
        #                 updateSelectInput(
        #                     session,
        #                     "fields",
        #                     choices=new_choices,
        #                     selected=input$fields
        #                 )
        #             })
        #             
        #             output$table = DT::renderDataTable({
        #                 
        #                 thedata() %>% 
        #                     datatable(options=list(
        #                         pageLength=10, 
        #                         scrollX=TRUE,
        #                         # autoWidth=TRUE,
        #                         columnDefs=list(list(width="10px", targets="_all"))
        #                     )) %>%
        #                     DT::formatRound(columns=input$fields, digits=input$decimals)
        #             })
        #             
        #             output$download <- downloadHandler(
        #                 filename = function() {"result_table.tsv"},
        #                 content = function(fname) {
        #                     write_tsv(thedata(), fname)
        #                 } 
        #             )
        #         },
        #         options=list(height=height)
        #     )
        #     
        # },
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
        
        update_input_range = function(session, datasets, target, dataset_name, column_names) {
            # private$update_input_range(session, datasets, "venn_thres", input$venn_type)
            # browser()
            col_data <- rowData(datasets[[dataset_name]])[column_names] %>% unlist()
            min_val <- round(min(col_data, na.rm=TRUE), digits=2)
            max_val <- round(max(col_data, na.rm=TRUE), digits=2)
            updateNumericInput(
                session,
                target,
                min=min_val,
                max=max_val
            )
        },
        
        get_contrasts_from_suffix = function(row_data_cols, contrast_suffix) {
            make.names(
                gsub(
                    paste0(".", contrast_suffix), "", 
                    row_data_cols[grepl(paste0(contrast_suffix, "$"), row_data_cols)]
                )
            )
        },
        get_suffixes_from_contrast = function(row_data_cols, contrast) {
            
            matching_fields <- row_data_cols[grepl(contrast, row_data_cols)]
            make.names(gsub(paste0(contrast, "."), "", matching_fields))
        }
    )
) 

masterwidget <- MasterWidget$new()
print("Loading module to 'masterwidget'")
