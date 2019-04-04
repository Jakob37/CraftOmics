library(R6)
library(outliers)
library(tidyverse)

ProteinRollup <- R6Class(
    public = list(
        qpline_inferno = function(qpline_path, peptide_fp, design_fp, output_dir, sample_header, group_header, peptide_header, protein_header, 
                                  quiet=FALSE, return_df=FALSE, overwrite_existing=FALSE) {
            
            # peptide_header and protein_header refers to columns in the raw data containing peptide information and protein information ('External.IDs' in Proteios)
            
            base_name <- gsub(".txt$", "", basename(peptide_fp))
            
            remove_path <- function(path) {
                message("Removing: ", path)
                unlink(path, recursive = TRUE)
            }
            
            if (overwrite_existing) {
                remove_path(paste0(output_dir))
            }
            
            options <- c(
                "--input", peptide_fp,
                "--groups", design_fp,
                "--outputdir", output_dir,
                "-sample_header", sample_header,
                "-group_header", group_header,
                "-peptide_header", peptide_header,
                "-protein_header", protein_header,
                "-run", "inferno"
            )
            
            # Note: This running is the Jupyter local Python, at least when executing notebook on Bruce
            command <- c("python3", qpline_path, options)
            if (!quiet) {
                message(paste("Running command: ", paste(command, collapse=" ")))
            }
            
            run_status <- system2(command, stdout=!quiet, stderr=!quiet, wait = TRUE)
            out <- lapply(run_status, function(line) {message(line)})
            
            if (return_df) {
                self$load_qpline_inferno_from_output_dir(peptide_fp, output_dir)
            }
        },
        qpline_help = function(qpline_path) {
            command <- c("python3", qpline_path, "--help")
            run_status <- system2(command, stdout=TRUE, stderr=TRUE, wait = TRUE)
            out <- lapply(run_status, function(line) { message(line) })
        },
        
        load_qpline_inferno_from_output_dir = function(peptide_fp, output_dir) {
            base_name <- gsub(".txt$", "", basename(peptide_fp))
            out_path <- paste0(output_dir, "/Inferno_Output/", base_name, "_RRollup_proteinreport.txt")
            if (!quiet) {
                message("Reading: ", out_path)
            }
            df <- read.csv(out_path, "\t")
            df
        },
        
        protein_rollup_on_matrix = function(design_fp, data_fp, protein_col, peptide_col, sample_col, out_fp, rollup_func="rrollup") {
            
            ddf <- read_tsv(design_fp, col_types=cols())
            raw_rdf <- read_tsv(data_fp, col_types=cols())
            rdf <- raw_rdf %>% filter(!is.na(UQ(as.name(protein_col))))
            
            trimmed_count <- nrow(raw_rdf) - nrow(rdf)
            if (trimmed_count != 0) {
                message("Trimmed ", trimmed_count, " proteins with missing IDs")
            }
            
            sdf <- rdf %>% dplyr::select(one_of(ddf[[sample_col]]))
            
            protein_data <- rdf %>% dplyr::select(protein_col) %>% unlist()
            peptide_data <- rdf %>% dplyr::select(peptide_col) %>% unlist()

            message("Performing rollup for ", length(peptide_data), " peptides to ", length(unique(protein_data)), " proteins")
            
            if (rollup_func == "rrollup") {
                rollup <- self$rrollup
            }
            else {
                stop("Unknown rollup func: ", rollup_func)
            }
            
            prot_sdf <- self$protein_rollup(peptide_data, protein_data, as.matrix(sdf), rollup_func=rollup)
            write_tsv(prot_sdf, path=out_fp)
        },
        
        # Inspired by DANTE and pmartR
        protein_rollup = function(pep_ids, protein_ids, pep_mat, rollup_func=self$rrollup) {
            
            warning("protein_rollup: Use on your own risk, mostly untested")
            
            # pep_expression <- cbind(pep_ids, pep_mat)
            # prot_expression <- cbind(protein_ids, pep_mat)
            # protein_map <- cbind(protein_ids, pep_ids)
            
            unique_proteins <- unique(protein_ids)
            pep_counts <- list()

            pep_results <- list()
            for (protein in unique_proteins) {
                
                row_inds <- which(protein_ids == protein)
                current_peps_only <- pep_mat[row_inds, , drop=FALSE]
                
                scaled_peps <- rollup_func(current_peps_only)
                pep_results[[protein]] <- scaled_peps
                pep_counts[[protein]] <- nrow(current_peps_only)
            }

            cbind(Protein=unique_proteins, Peptides=unlist(pep_counts), data.frame(do.call("rbind", pep_results)))
        },
        
        # Based on InfernoRDN
        # https://github.com/PNNL-Comp-Mass-Spec/InfernoRDN/blob/master/Rscripts/Rollup/RRollup.R
        remove_outliers = function(pep_mat, minPs=5, pvalue_thres=0.05) {
            
            warning("Completely untested")
            
            xPeptideCount <- colSums(!is.na(pep_mat))
            
            for (sample_i in 1:dim(pep_mat)[2]) {
                if (xPeptideCount[sample_i] >= minPs) {
                    
                    # Repeat the outlier check until none are classified as outliers
                    repeat {
                        grubbs <- outliers::grubbs.test(pep_mat[, sample_i])
                        if ( (grubbs$p.value < pvalue_thres) && (!is.nan(grubbs$statistic[2])) && grubbs$statistic[2] != 0) {
                            # Fill outlier with median value
                            pep_mat[, sample_i] <- outliers::rm.outlier.1(pep_mat[, sample_i], fill=TRUE, median=TRUE)
                        }
                        else {
                            break
                        }
                    }
                }
            }
        },
        
        zrollup = function(peptides, combine_func=median) {
            warning("Relatively untested, and no Grubbs outlier test performed")
            
            num_peps <- nrow(peptides)
            #res <- matrix(NA, nrow=1, ncol=ncol(peptides))
            
            # Compute mean and sd of peptides
            mds <- matrixStats::rowMedians(peptides, na.rm=TRUE)
            sds <- matrixStats::rowSds(peptides, na.rm=TRUE)
            
            # Scale peptide data as pep_scaled = (pep - median) / sd
            medians_mat <- matrix(mds, nrow=num_peps, ncol=ncol(peptides), byrow=FALSE)
            standiv_mat <- matrix(sds, nrow=num_peps, ncol=ncol(peptides), byrow=FALSE)
            
            proteins_scaled <- apply((peptides - medians_mat) / standiv_mat, 2, combine_func, na.rm=TRUE)
            proteins_scaled
        },
        
        rrollup = function(peptides, combine_func=median) {
            # warning("Relatively untested, and no Grubbs outlier test performed")
            
            num_peps <- nrow(peptides)

            if (num_peps == 1) {
                protein_val <- unlist(peptides)
            }
            else {
                # 1. Select reference peptide
                na.cnt <- apply(is.na(peptides), 1, sum)
                least.na <- which(na.cnt == min(na.cnt))
                
                # If tied, select with highest median abundance
                if (length(least.na) > 1) {
                    mds <- matrixStats::rowMedians(peptides, na.rm=TRUE)[least.na]
                    least.na <- least.na[which(mds == max(mds))]
                }
                prot_val <- unlist(peptides[least.na, ])
                
                # 2. Ratio all peptides to the reference
                # Since data is on log scale this is the difference (?)
                ref_ratios <- matrix(prot_val, nrow=num_peps, ncol=ncol(peptides), byrow=TRUE) - peptides
                scaling_factor <- matrixStats::rowMedians(ref_ratios, na.rm=TRUE)
                
                # 3. Use median of ratio as scaling factor for each peptide
                x_scaled <- peptides + matrix(scaling_factor, nrow=num_peps, ncol=ncol(peptides))

                # 4. Set abundance as median peptide abundance
                protein_val <- apply(x_scaled, 2, combine_func, na.rm=TRUE)
            }
            
            protein_val
        },
        
        # qrollup_thres: 0 - 1 value, peptides above threshold used for rollup
        qrollup = function(peptides, qrollup_thres, combine_func=median) {
            
            warning("Relatively untested, and no Grubbs outlier test performed")
            num_peps <- nrow(peptides)
            if (num_peps == 1) {
                protein_val <- unlist(peptides)
                peps_used <- 1
            }
            else {
                # Step 1: Subset peptides with abundance >= qrollup_threshold
                means <- Matrix::rowMeans(peptides, na.rm=TRUE)
                quantil <- quantile(means, probs=qrollup_thres, na.rm=TRUE)
                
                quality_peps <- peptides[means >= quantil, ]
                peps_used <- nrow(quality_peps)
                
                # Step 1b: If only 1 peptide, set protein value to that
                if (nrow(quality_peps) == 1) {
                    protein_val <- unlist(quality_peps)
                }
                # Step 2: Set protein abundance to mean / median of peptide abundances
                else {
                    protein_val <- apply(quality_peps, 2, combine_func, na.rm=TRUE)
                }
                
                protein_val
            }
        }
    ),
    private = list(
        rollup.score = function(currPepSel, currProtSel, method) {
            
            N1 <- dim(currPepSel)[1]
            N2 <- dim(currPepSel)[2]
            pepCorr <- rep(numeric(0),N1)
            ws <- rep(numeric(0),N1)
            
            for(i in 1:N1)#finds correlation values between each peptide profile and calculated protein profile
            {
                ws[i] <- sum(!is.na(currPepSel[i,]))/N2
                if(method=="pearson")
                    pepCorr[i] <- cor(as.vector(currPepSel[i,]), as.vector(currProtSel), use="pairwise.complete.obs")
                else if(method=="kendall")
                    pepCorr[i] <- cor(as.vector(currPepSel[i,]), as.vector(currProtSel), use="pairwise.complete.obs",
                                      method="kendall")
                else
                    pepCorr[i] <- cor(as.vector(currPepSel[i,]), as.vector(currProtSel), use="pairwise.complete.obs",
                                      method="spearman")
            }
            #meanCorr <- mean(pepCorr, na.rm=TRUE) #mean correlation value for each protein
            meanCorr <- weighted.mean(pepCorr, ws, na.rm=TRUE) #mean correlation value for each protein
            Penalty1 <- 1-1/N1
            #Penalty2 <- sum(!is.na(currPepSel))/(N1*N2)
            #Score <- meanCorr * Penalty1 * Penalty2
            Score <- meanCorr * Penalty1
            
            out <- Score
            return(out)
        },
        
        dante_remove.outliers = function(Data, minPs=5, pvalue=0.05)
            # internal function used by normalize.proteins()
            # Calculates the protein abundances after removing outliers
            # Depends on package "outliers"
        {
            # library(outliers)
            # ColNames = colnames(Data)
            xPeptideCount <- colSums(!is.na(Data))
            for (i in 1:dim(Data)[2])
            {
                if (xPeptideCount[i] >= minPs)
                {
                    repeat
                    {
                        grubbs <- grubbs.test(Data[,i]) # Grubb's test
                        if ( (grubbs$p.value < pvalue) && (!is.nan(grubbs$statistic[2])) &&
                             (grubbs$statistic[2] != 0) ) # pass the p-value cutoff
                        {
                            Data[,i] <- rm.outlier(Data[,i], fill=TRUE, median=TRUE)
                            # fill the outlier with the median value
                        }
                        else { break }
                    }
                }
            }
            return(Data)
        },
        
        rm.outlier = function(data, fill, median) {
            
        }
    )
)


protroll <- ProteinRollup$new()
pr <- protroll
print("Loading module to 'protroll' and 'pr'")
