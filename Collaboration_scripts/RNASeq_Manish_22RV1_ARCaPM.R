# In this case, the replicates for each group are clustered together.
# So, the within-group variability is low. Hence, we do not subset the data

DEG.params  <- list(contrast    = c("ARCaPM.NDRG1.4Gy-ARCaPM.NDRG1.0Gy",
                                    "ARCaPM.PRDX1.4Gy-ARCaPM.PRDX1.0Gy",
                                    "ARCaPM.WT.4Gy-ARCaPM.WT.0Gy",
                                    "RV1_22.NDRG1.4Gy-RV1_22.NDRG1.0Gy",
                                    "RV1_22.PRDX1.4Gy-RV1_22.PRDX1.0Gy",
                                    "RV1_22.WT.4Gy-RV1_22.WT.0Gy",
                                    "ARCaPM.NDRG1.4Gy-ARCaPM.WT.4Gy",
                                    "ARCaPM.PRDX1.4Gy-ARCaPM.WT.4Gy",
                                    "ARCaPM.NDRG1.0Gy-ARCaPM.WT.0Gy",
                                    "ARCaPM.PRDX1.0Gy-ARCaPM.WT.0Gy",
                                    "RV1_22.NDRG1.4Gy-RV1_22.WT.4Gy",
                                    "RV1_22.PRDX1.4Gy-RV1_22.WT.4Gy",
                                    "RV1_22.NDRG1.0Gy-RV1_22.WT.0Gy",
                                    "RV1_22.PRDX1.0Gy-RV1_22.WT.0Gy",
                                    "(ARCaPM.NDRG1.4Gy-ARCaPM.NDRG1.0Gy)-(ARCaPM.WT.4Gy-ARCaPM.WT.0Gy)",
                                    "(RV1_22.NDRG1.4Gy-RV1_22.NDRG1.0Gy)-(RV1_22.WT.4Gy-RV1_22.WT.0Gy)"),
                                    design      = "Comparisons",
                                    design.ref  = c("Cell.Line:ARCaPM", "Condition:WT", "Treatment:0Gy"),
                                    lfc.cutoff  = 0,
                                    padj.cutoff = 0.1,
                                    deseq2.batch.correct = FALSE,
                                    proj        = "RNASeq_Manish_22RV1_ARCaPM",
                                    species     = "Homo sapiens")

heatmap.params <- list(anno.row       = NULL,        # NULL, c("Group")
                       anno.column    = c("Treatment", "Cell.Line", "Condition"),
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
