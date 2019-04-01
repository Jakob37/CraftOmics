library(ggplot2)
library(ggdendro)
library(ggfortify)
library(RColorBrewer)
library(R6)
library(MASS)

MultivarVis <- R6Class(
    public = list(

        style_plt = function(plt) {
            palette <- "Dark2"
            plt <- plt + 
                theme_classic()
                # scale_color_brewer(palette=palette) + 
                # scale_fill_brewer(palette=palette)
            plt
        },
        
        pca = function(expr_m, color_factor, title="PCA", pcs=c(1,2), label=NULL) {

            # Inspired by:
            # https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
            
            expr_m_nona <- expr_m[complete.cases(expr_m), ]
            pca_obj <- prcomp(t(expr_m_nona), scale=TRUE, center=TRUE)
            style_df <- data.frame(color=color_factor)
            
            if (!is.null(label)) {
                rownames(style_df) <- make.names(label, unique=TRUE)
            }

            if (is.null(label)) shape <- NULL else shape <- FALSE

            plt <- autoplot(
                pca_obj, 
                data=style_df, 
                colour="color", 
                label=!is.null(label), 
                shape=shape,
                loadings=FALSE, 
                loadings.label=FALSE, 
                x=pcs[1], 
                y=pcs[2]
            ) + ggtitle(title)
            
            plt <- self$style_plt(plt)
            return(plt)
        },
        
        lda = function(data, color_factor, title="LDA", lds=c(1,2), label=NULL) {
            
            
            
            lda_obj <- MASS::lda(
                
            )
            
        },

        plotMDS = function(expr_m, levels, comp1=1, comp2=2, title="no title") {

            labels <- colnames(expr_m)
            d <- stats::dist(scale(t(stats::na.omit(expr_m)), center=TRUE, scale=TRUE))
            fit <- stats::cmdscale(d, eig=TRUE, k=2)
            x <- fit$points[, comp1]
            y <- fit$points[, comp2]
            graphics::plot(x, y, type="n", main=title, xlab="", ylab="")
            graphics::text(fit$points[, 1], fit$points[, 2], col=levels, labels=labels)
        },

        get_component_fraction = function(expr_m) {

            # Retrieves a vector with percentage contributions for each PC
            
            expr_m_nona <- expr_m[complete.cases(expr_m),]
            pca_object <- prcomp(t(expr_m_nona), scale=TRUE, center=TRUE)
            # pca_object <- self$get_pca_object(expr_m)

            percentVar <- pca_object$sdev^2 / sum(pca_object$sdev^2 )
            names(percentVar) <- colnames(pca_object$x)
            return(percentVar)
        },

        plot_component_fraction = function(expr_m, max_comps=NULL) {

            # Directly outputs the PCA numbers together with PC fractions

            comp_perc <- self$get_component_fraction(expr_m) * 100
            if (!is.null(max_comps)) {
                comp_perc <- head(comp_perc, max_comps)
            }
            
            plot_df <- data.frame(x=paste0("", seq_len(length(comp_perc))), y=comp_perc)
            ggplot(plot_df, aes(x, y)) + geom_bar(stat="identity") + 
                theme_classic() +
                ggtitle("PCA loadings") + ylab("Variance (%)") + xlab("Principal component")
            # plot(100 * comp_perc, main="PC loadings", xlab="PC", ylab="Perc. var")
        },

        dendogram = function(data_m, color_levels, labels=NULL, pick_top_variance=null, title="Dendogram") {

            samples <- colnames(data_m)
            
            if (is.null(labels)) {
                labels <- samples
            }
            
            # Setup data
            expr_m_nona <- data_m[complete.cases(data_m),]

            # Calculate tree
            scaledTransposedMatrix <- scale(t(expr_m_nona), center=TRUE, scale=TRUE)
            hc <- stats::hclust(stats::dist(scaledTransposedMatrix), "ave")
            dhc <- as.dendrogram(hc)
            # Note - Label order is shuffled within this object! Be careful with coloring.
            ddata <- dendro_data(dhc, type="rectangle")

            # Prepare for plotting
            cluster_label_order <- match(ddata$labels$label, samples)
            ddata$labels$color <- color_levels[cluster_label_order]
            ddata$labels$label <- labels[cluster_label_order]

            # Visualize
            plt <- ggplot(segment(ddata)) +
                geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
                theme_dendro() +
                geom_text(data=label(ddata),
                          aes(x=x, y=y, label=label, color=color),
                          vjust=0.5, hjust=0, size=3) +
                coord_flip() +
                scale_y_reverse(expand=c(0.2, 0)) +
                scale_x_continuous(expand=c(0,1)) +
                ggtitle(title)
            plt
        },
        
        table_heatmap = function(data_m, row_annot) {
            
            df <- cbind(annot=row_annot, data_m)
            df.melt <- df %>% gather("Sample", "Level", -annot)
            plt <- ggplot(df.melt, aes(Sample, annot)) + 
                geom_tile(aes(fill=Level, color="white")) +
                scale_fill_gradient(low="white", high="steelblue")
            plt
        }
    )
)

multvis <- MultivarVis$new()
mv <- multvis
print("Loading module to 'multvis' and 'mv'")
