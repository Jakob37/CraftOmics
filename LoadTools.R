library(R6)

library(tidyverse)
library(SummarizedExperiment)

#' Utilities for loading mixed data/annotation matrices with separate design matrix 
#' 
#' load_summarized_experiment: Utility to load mixed data/annotation and design matrix
#' directly into a SummarizedExperiment object. If provided as a list it is automatically
#' loaded as a multiple of assays
LoadTools <- R6Class(
    
    public = list(
        
        load_summarized_experiment_from_path = function(ddf_fp, rdf_fp, sample_col, na_field="NA", comment_char="", silent=FALSE) {
            
            if (silent) {
                col_types <- cols()
            }
            else {
                col_types <- NULL
            }
            
            ddf <- read_tsv(ddf_fp, col_types=col_types)
            rdf <- read_tsv(rdf_fp, na=na_field, comment=comment_char, col_types=col_types)
            se <- self$load_summarized_experiment(ddf, rdf, sample_col)
            se
        },
        
        load_summarized_experiment = function(ddf, rdf, sample_col, extra_rdfs=NULL) {
            
            if (!(sample_col %in% colnames(ddf))) {
                stop("Target sample col: ", sample_col, " not found in design matrix headers: ", paste(ddf, collapse=", "))
            }
            
            sample_names <- as.character(ddf[[sample_col]])
            sdf <- rdf %>% select(sample_names)
            adf <- rdf %>% select(-one_of(sample_names))
            
            sdf_list <- list()
            sdf_list[["raw"]] <- as.matrix(sdf)
            
            if (!is.null(extra_rdfs)) {
                for (name in names(extra_rdfs)) {
                    extra_rdf <- extra_rdfs[[name]]
                    extra_sdf <- extra_rdf %>% select(sample_names)
                    sdf_list[[name]] <- as.matrix(extra_sdf)
                }
            }
            
            if (typeof(sdf_list[[1]]) != "double") {
                stop("Expected format for assay is 'double', found: ", typeof(sdf_list[[1]]))
            }
            
            SummarizedExperiment::SummarizedExperiment(
                assays=sdf_list,
                colData=ddf,
                rowData=adf
            )
        },
        
        split_dataset = function(data_fp, design_fp, split1_data_fp, split1_design_fp, 
                       split2_data_fp, split2_design_fp, sample_col, split_col) {
            
            rdf <- readr::read_tsv(data_fp, col_types=readr::cols())
            ddf <- readr::read_tsv(design_fp, col_types=readr::cols())
            sdf <- rdf %>% dplyr::select(dplyr::one_of(ddf[[sample_col]]))
            adf <- rdf %>% dplyr::select(-dplyr::one_of(ddf[[sample_col]]))
            
            split_levels <- sort(unique(ddf[[split_col]]))

            if (length(split_levels) != 2) {
                stop("This function is designed for two split levels, found: ", paste(split_levels, collapse=", "))
            }
            
            s1_ddf <- ddf %>% filter(UQ(as.name(split_col)) == split_levels[1]) %>% data.frame()
            s2_ddf <- ddf %>% filter(UQ(as.name(split_col)) == split_levels[2]) %>% data.frame()
            
            s1_sdf <- sdf %>% select(s1_ddf[[sample_col]])
            s2_sdf <- sdf %>% select(s2_ddf[[sample_col]])
            
            s1_rdf <- cbind(adf, s1_sdf)
            s2_rdf <- cbind(adf, s2_sdf)
            
            write_tsv(s1_rdf, path=split1_data_fp)
            write_tsv(s1_ddf, path=split1_design_fp)
            write_tsv(s2_rdf, path=split2_data_fp)
            write_tsv(s2_ddf, path=split2_design_fp)
        },
        
        write_derep_matrices = function(rdf, ddf, out_data_path, namecol, techrepcol, out_design_path=NULL, trim_na_col=NULL) {
            
            if (!is.null(trim_na_col)) {
                trimmed_rdf <- rdf %>% filter(!is.na(UQ(as.name(trim_na_col))))
                message("Trimmed before and after counts: ", nrow(rdf), ", ", nrow(trimmed_rdf))
                rdf <- trimmed_rdf
            }
            
            sdf <- rdf %>% select(as.character(ddf[[namecol]]))
            adf <- rdf %>% select(-one_of(as.character(ddf[[namecol]])))
            red_mats <- st$reduce_technical_replicates_for_matrices(designMat = ddf, dataMat = sdf, ddf[[techrepcol]])
            message("Writing ", nrow(red_mats$data), " to ", out_data_path)
            write_tsv(cbind(adf, data.frame(red_mats$data)), path=out_data_path)
            if (!is.null(out_design_path)) {
                message("Writing ", nrow(red_mats$design), " to ", out_design_path)
                write_tsv(red_mats$design, path=out_design_path)
            }
        },
        
        make_ses = function() {
            
        }
    )
)

lt <- LoadTools$new()
print("Loading module to 'lt'")
