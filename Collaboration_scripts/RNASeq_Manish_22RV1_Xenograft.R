DEG.params  <- list(contrast    = c("SPB-Vehicle",
                                    "SPB.IR-Vehicle",
                                    "IR-Vehicle"),
                    design      = "Comparisons",
                    design.ref  = c("Condition:Vehicle"),
                    lfc.cutoff  = 0,
                    padj.cutoff = 0.1,
                    deseq2.batch.correct = FALSE,
                    proj        = "RNASeq_Manish_22RV1_Xenograft",
                    species     = "Homo sapiens")

heatmap.params <- list(anno.row       = NULL,        # NULL, c("Group")
                       anno.column    = c("Condition"),
                       row.split      = NA,     
                       col.split      = NA, #c("Treatment"),
                       row.cluster    = c("all"),           # c("alphabetical", "group", "all")
                       col.cluster    = c("all"),  # c("alphabetical", "group", "all")
                       discrete_panel = FALSE, 
                       log.transform  = TRUE,
                       scale          = TRUE,
                       border_color   = "white",
                       bar_width      = NA,              # NA , 5
                       bar_height     = NA,              # NA , 5
                       width          = 5,              # NA
                       height         = 5,              # NA 
                       matrix_color   = "rdbu",          # c("vrds", "rdbu")
                       expr_legend    = TRUE,  # set FALSE if it overlaps with annotation legends
                       file_format    = "tiff")
