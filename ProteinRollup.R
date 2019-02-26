library(R6)

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
        
        # Inspired by DANTE and pmartR
        protein_rollup = function(pep_ids, protein_ids, pep_mat, rollup_func=self$zrollup) {
            
            pep_expression <- cbind(pep_ids, pep_mat)
            prot_expression <- cbind(protein_ids, pep_mat)
            protein_map <- cbind(pep_ids, protein_ids)
            
            unique_proteins <- unique(protein_ids)

            pep_results <- list()
            for (protein in unique_proteins) {
                
                row_inds <- which(protein_ids == protein)
                current_peps_only <- pep_mat[row_inds, ]
                
                scaled_peps <- rollup_func(current_peps_only)
                pep_results[[protein]] <- scaled_peps
            }
            
            do.call("rbind", pep_results)
        },
        
        zrollup = function(peptides, combine_func=median) {
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
            num_peps <- nrow(peptides)
            #res <- matrix(NA, nrow=1, ncol=ncol(peptides))
            
            if (num_peps == 1) {
                protein_val <- unlist(peptides)
            }
            else {
                # 1. Select reference peptide
                na.cnt <- apply(is.na(peptides), 1, sum)
                least.na <- which(na.cnt == min(na.cnt))
                
                # If tied, select with highest median abundance
                if (length(least.na) > 1) {
                    mds <- matrixStats::rowMedian(peptides, na.rm=TRUE)[least.na, ]
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
    private = list()
)


protroll <- ProteinRollup$new()
pr <- protroll
print("Loading module to 'protroll' and 'pr'")
