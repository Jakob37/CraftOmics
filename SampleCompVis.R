library(VennDiagram)
library(gridExtra)
library(R6)

SampleCompVis <- R6Class(
    public = list(

        simple_venn = function(col1, col2, dereplicate=TRUE, title="Venn", colors=NULL, verbose=FALSE, labels=c("left", "right"), legend.position="right") {
            
            if (dereplicate) {
                col1 <- unique(col1)
                col2 <- unique(col2)
            }
            
            df.venn <- data.frame(
                x = c(-0.5, 0.5), y = c(0, 0),
                labels = labels
            )
            
            joint_count <- length(intersect(col1, col2))
            left_count <- length(col1) - joint_count
            right_count <- length(col2) - joint_count
            # ses$prot_des_prot %>% get_uniques_from_rdf("Protein")
            if (verbose) {
                message("Overlapping features: ", intersect(col1, col2))
                message("Left only features: ", setdiff(col1, col2))
                message("Right only features: ", setdiff(col2, col1))
            }
            
            df.vdc <- data.frame(
                x=c(-1.5, 1.5, 0),
                y=c(0, 0, 0),
                label=c(left_count, right_count, joint_count)
            )

            self$do_ggvenn(df.venn, df.vdc, legend.position=legend.position, title=title)
        },
        
        count_venn = function(left_count, right_count, joint_count, title="Venn", colors=NULL, verbose=FALSE, labels=c("left", "right"), legend.position="right") {
            
            df.venn <- data.frame(
                x = c(-0.5, 0.5), y = c(0, 0),
                labels = labels
            )
            
            df.vdc <- data.frame(
                x=c(-1.5, 1.5, 0),
                y=c(0, 0, 0),
                label=c(left_count-joint_count, right_count-joint_count, joint_count)
            )

            self$do_ggvenn(df.venn, df.vdc, legend.position=legend.position, title=title)

        },
        
        do_ggvenn = function(df.venn, df.vdc, legend.position='none', colors=NULL, title="No title") {
            
            
            if (is.null(colors)) {
                colors <- c("#A0D398", "#FDF9BD")
            }
            
            plt <- ggplot2::ggplot(data=df.venn) +
                ggforce::geom_circle(
                    ggplot2::aes_string(x0 = "x", y0 = "y", r = 1.5, fill = "labels"),
                    alpha = 0.3,
                    size = 0.5,
                    colour = 'darkgray'
                ) +
                ggplot2::coord_fixed() +
                ggplot2::theme_void() +
                ggplot2::scale_fill_manual(values = colors) +
                ggplot2::scale_colour_manual(values = colors, guide = FALSE) +
                ggplot2::labs(fill = NULL) +
                ggplot2::annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$label, size = 5) +
                ggplot2::ggtitle(title) + 
                theme(plot.title = element_text(hjust = 0.5))
            
            if (!is.null(labels)) {
                plt <- plt + ggplot2::theme(legend.position = legend.position, legend.direction='vertical')
            }
            
            plt
        },
        
        do_paired_expression_venn = function(adf, contrast1_base, contrast2_base, sig_thres=0.05, fold_thres=NULL, title="", colors=c("#AAAAAA", "#0C81A2")) {

            warning("Known to mix up coloring sides of Venn, double check numbers")
            
            if (is.null(fold_thres)) {
                contrast1_sig_contrast <- adf[[paste0(contrast1_base, ".adj.P.Val")]] < sig_thres
                contrast2_sig_contrast <- adf[[paste0(contrast2_base, ".adj.P.Val")]] < sig_thres
            }
            else {
                contrast1_sig_contrast <- adf[[paste0(contrast1_base, ".adj.P.Val")]] < sig_thres & abs(adf[[paste0(contrast1_base, ".logFC")]]) > fold_thres
                contrast2_sig_contrast <- adf[[paste0(contrast2_base, ".adj.P.Val")]] < sig_thres & abs(adf[[paste0(contrast2_base, ".logFC")]]) > fold_thres
            }
            
            contrast1_fold <- adf[[paste0(contrast1_base, ".logFC")]]
            contrast2_fold <- adf[[paste0(contrast2_base, ".logFC")]]
            
            both_sig_inds <- which(contrast1_sig_contrast & contrast2_sig_contrast)
            contra_count <- length(which(sign(contrast1_fold[both_sig_inds]) != sign(contrast2_fold[both_sig_inds])))
            
            print(table(contrast1_sig_contrast))
            print(table(contrast2_sig_contrast))
            
            vdc <- limma::vennCounts(cbind(contrast1_sig_contrast, contrast2_sig_contrast))
            class(vdc) <- 'matrix'
            df.vdc <- data.frame(vdc[-1, ])
            df.vdc <- rbind(df.vdc, c(-1, -1, contra_count))
            df.vdc$x <- c(1.5, -1.5, 0, 0)
            df.vdc$y <- c(0, 0, 0, -0.3)
            df.vdc$label <- df.vdc$Counts
            df.vdc$label[3] <- paste0(df.vdc$label[3] - df.vdc$label[4], " same")
            df.vdc$label[4] <- paste0(df.vdc$label[4], " contra")
            
            df.venn <- data.frame(x = c(-0.5, 0.5),
                                  y = c(0, 0),
                                  labels = c(contrast1_base, contrast2_base))
            
            print(df.venn)
            
            plt <- ggplot2::ggplot(data=df.venn) +
                ggforce::geom_circle(
                    ggplot2::aes_string(x0 = "x", y0 = "y", r = 1.5, fill = "labels"), 
                    alpha = .3, 
                    size = 0.5, 
                    colour = 'darkgray'
                ) +
                ggplot2::coord_fixed() +
                ggplot2::theme_void() +
                ggplot2::theme(legend.position = 'bottom', legend.direction='vertical') +
                ggplot2::scale_fill_manual(values = colors) +
                ggplot2::scale_colour_manual(values = colors, guide = FALSE) +
                ggplot2::labs(fill = NULL) +
                ggplot2::annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$label, size = 5) +
                ggplot2::ggtitle(title)
            
            plt
        },
        
        make_threeway_venn = function(
            id1_lab, id2_lab, id3_lab, id1, id2, id3, id12, id13, id23, id123, 
            show_perc=FALSE, dummy_path="/dev/null", title="Venn diagram", title_ypos=0.93, colors=c("blue", "red", "green"), scaled=FALSE) {

            pdf(file=dummy_path)
            venn.plot <- draw.triple.venn(
                area1 = id1,
                area2 = id2,
                area3 = id3,
                n12 = id12,
                n23 = id23,
                n13 = id13,
                n123 = id123,
                category = c(id1_lab, id2_lab, id3_lab),
                print.mode = perc,
                fill = colors,
                scaled=scaled)
            dev.off()
            
            venn_gtree <- gTree(children=gList(children=venn.plot, textGrob(title, y=unit(title_ypos, "npc"))))
            
            return(venn_gtree)
        },

        make_pairwise_venn = function(
            id1_lab, id2_lab, id1, id2, id12, 
            show_perc=FALSE, dummy_path="/dev/null", title="Venn diagram", title_ypos=0.93, colors=c("red", "blue"), scaled=FALSE) {

            labs <- c(id1_lab, id2_lab)

            if (!show_perc) {
                perc <- "raw"
            }
            else {
                perc <- c("raw", "percent")
            }
            
            if (id2 > id1) {
                rotation <- 180
            }
            else {
                rotation <- 0
            }

            pdf(file=dummy_path)
            venn.plot <- draw.pairwise.venn(
                area1 = id1,
                area2 = id2,
                cross.area = id12,
                cex = 1,
                lwd = 0.5,
                category = labs,
                cat.pos = c(180, 180),
                cat.cex = 0.9,
                fill = colors,
                print.mode = perc,
                rotation.degree = rotation,
                scaled = scaled)
            dev.off()
            
            venn_gtree <- gTree(children=gList(children=venn.plot, textGrob(title, y=unit(title_ypos, "npc"))))
            
            return(venn_gtree)
        },

        do_col_venn = function(col1, col2, name1, name2, sig_thres=0.1, show_perc=F, dummy_path='/dev/null', ...) {

            intersect <- intersect(col1, col2)

            out1 <- length(col1)
            out2 <- length(col2)
            split <- length(intersect)

            self$make_pairwise_venn(name1, name2, out1, out2, split, show_perc=show_perc, dummy_path=dummy_path, ...)
        },

        do_limma_venn = function(table1, table2, name1, name2, sig_thres=0.1, sig_col="adj.P.Val", dummy_path='/dev/null', show_perc=F, assign_rows=TRUE, ...) {

            if (assign_rows) {
                
                if (nrow(table1) != nrow(table2)) {
                    stop("Different number of rows found for table 1 and 2: ", nrow(table1), ", ", nrow(table2))
                }
                
                warning("Assumes that input tables are of same order")
                rownames(table1) <- as.character(seq_len(nrow(table1)))
                rownames(table2) <- as.character(seq_len(nrow(table2)))
            }
            
            table1_names <- rownames(table1[which(table1[, sig_col] < sig_thres),])
            table2_names <- rownames(table2[which(table2[, sig_col] < sig_thres),])
            self$do_col_venn(table1_names, table2_names, name1, name2, sig_thres=sig_thres, dummy_path=dummy_path, show_perc=show_perc, ...)
        }
    ),
    private = list()
)

samplecompvis <- SampleCompVis$new()
scv <- samplecompvis
print("Loading module to 'samplecompvis' and 'scv'")
