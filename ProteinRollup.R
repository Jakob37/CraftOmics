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
        }
    ),
    private = list()
)


protroll <- ProteinRollup$new()
pr <- protroll
print("Loading module to 'protroll' and 'pr'")
