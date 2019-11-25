library(R6)
library(tidyverse)
library(xlsx)

library(xlsx)


ExcelTools <- R6(
    public = list(
            write_commented_excel = function(table, comment_df, out_fp) {
                
                wb <- xlsx::createWorkbook()
                sh <- xlsx::createSheet(wb, "Sheet1")
                addDataFrame(
                    comment_df,
                    sheet = sh,
                    startRow = 1,
                    row.names = FALSE,
                    col.names = FALSE
                )
                addDataFrame(
                    table,
                    sheet = sh,
                    startRow = nrow(comments) + 2,
                    row.names = FALSE
                )
                
                xlsx::saveWorkbook(wb, file=out_fp)
            },
                
            write_excel_with_comment_from_file = function(table, comment_file, out_fp) {
                
                comments_df <- read_tsv(comment_file, col_names = FALSE) %>% data.frame()
                self$write_commented_excel(table, comments_df, out_fp)
            }
    )
)