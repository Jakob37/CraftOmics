library(R6)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggthemes)

VisualizationUtils <- R6Class(
    public = list(
        pltsize = function(width, height) {
            options(repr.plot.width=width, repr.plot.height=height)
        },

        multiplot = function(..., plotlist=NULL, cols=1, layout=NULL) {

            # Make a list from the ... arguments and plotlist
            plots <- c(list(...), plotlist)

            numPlots = length(plots)

            # If layout is NULL, then use 'cols' to determine layout
            if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
            }

            if (numPlots==1) {
                print(plots[[1]])

            }
            else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                    # Get the i,j matrix positions of the regions that contain this subplot
                    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

                    print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                    layout.pos.col = matchidx$col))
                }
            }
        },
        
        grobplot = function(plotlist, ncol=1, file=NULL, headAsList=FALSE, size=NULL, ylim=NULL, xlim=NULL, title=NULL, legend_title=NULL, return_val=FALSE, header=NULL) {
            
            if (!is.null(header)) {
                warning("Option 'header' deprecated, use 'title' instead")
                title <- header
            }
            
            if (!is.null(size)) {
                self$pltsize(size[1], size[2])
            }
            
            annotate_plots <- function(plotname, plotlist, headAsList=FALSE, xlim=NULL, ylim=NULL, legend_title=NULL) {
                
                plt <- plotlist[[plotname]]
                if(headAsList) {
                    plot_head <- plotname
                }
                else {
                    plot_head <- NULL
                }
                plt <- self$annotate_plot(plt, plotname=plot_head, xlim=xlim, ylim=ylim, legend=legend_title)
                plt
            }
            
            if (is.null(names(plotlist))) {
                names(plotlist) <- seq_len(length(plotlist))
            }
            
            plotlist <- lapply(names(plotlist), annotate_plots, plotlist=plotlist, headAsList=headAsList, ylim=ylim, xlim=xlim, legend_title=legend_title)
            
            grobs <- arrangeGrob(grobs=plotlist, ncol=ncol, top=title)
            grid.draw(grobs)

            if (!is.null(file)) {
                message("Writing plot to path: ", file)
                ggsave(grobs, file=file, width=size[1], height=size[2])
            }
            
            if (return_val) {
                grobs
            }
        },

        save_plotlist = function(plotlist, filePath, cols=1, width=4, height=4, type="pdf", res=150) {
            if (type == "pdf") {
                pdf(filePath, width=width, height=height)
                self$multiplot(plotlist=plotlist, cols=cols)
                dev.off()
            }
            else if (type == "png") {
                png(filePath, width=width, height=height, units="in", res=res)
                self$multiplot(plotlist=plotlist, cols=cols)
                dev.off()
            }
            else {
                stop(paste("Unknown file type supplied:", type))
            }
            print(paste("Figure written to:", filePath))
        },

        annotate_plot = function(plt, plotname=NULL, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, legend=NULL, legend_cols=1, number_colors=NULL, number_fill=NULL) {

            if (!is.null(plotname)) {
                plt <- plt + ggtitle(plotname)
            }
            
            if (!is.null(xlab)) {
                plt <- plt + xlab(xlab)
            }

            if (!is.null(ylab)) {
                plt <- plt + ylab(ylab)
            }
            
            if (!is.null(xlim)) {
                plt <- plt + xlim(xlim)
            }
            
            if (!is.null(ylim)) {
                plt <- plt + ylim(ylim)
            }

            if (!is.null(legend)) {
                plt <- plt + guides(fill=guide_legend(title=legend, ncol=legend_cols))
                plt <- plt + guides(color=guide_legend(title=legend, ncol=legend_cols))
            }

            if (!is.null(number_colors)) {
                getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
                plt <- plt + scale_color_manual(values=getPalette(number_colors))
            }

            if (!is.null(number_fill)) {
                getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
                plt <- plt + scale_fill_manual(values=getPalette(number_fill))
            }
            plt
        },
        
        active_levels_factor = function(full_factor, active_levels_list) {
            
            active_inds_list <- list()
            
            for (level_name in names(active_levels_list)) {
                active_level <- active_levels_list[[level_name]]
                active_inds <- which(full_factor %in% active_level)
                
                for (ind in active_inds) {
                    active_inds_list[[as.character(ind)]] <- level_name
                }
            }
            
            get_level <- function(index, active_inds_list) {
                if (index %in% names(active_inds_list)) {
                    active_inds_list[[as.character(index)]]
                }
                else {
                    "other"
                }
            }
            
            sub_level <- vapply(seq_len(length(full_factor)), get_level, "", active_inds_list=active_inds_list)
            sub_level
        }
    ),
    private = list()
)


visutil <- VisualizationUtils$new()
vu <- visutil
print("Loading module to 'visutil' and 'vu'")
