library(VennDiagram)
library(gridExtra)
library(R6)

SampleCompVis <- R6Class(
    public = list(

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
