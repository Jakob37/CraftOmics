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
        }
    )
)

lt <- LoadTools$new()
print("Loading module to 'lt'")
