library(R6)
library(grid)
library(gridExtra)

VisualizationUtils <- R6Class(
    public = list(
        
        generate_pdf_from_list = function(plot_list, file) {
            
            # Figures directly printed
            # Room for annotation
            
            pdf(file=file, paper="a4")  
            
            for (plt_i in seq_len(length(plot_list))) {
                plt <- plot_list[[plt_i]]
                grid.draw(plt)
                if (plt_i < length(plot_list)) {
                    grid.newpage()
                }
            }
            
            dev.off() 
            
            message("Report written to: ", file)
        }
        
    ),
    private = list()
)


reportgen <- VisualizationUtils$new()
rg <- reportgen
print("Loading module to 'reportgen' and 'rg'")
