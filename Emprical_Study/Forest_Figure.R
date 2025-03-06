
library(tidyverse)
library(forestplot)
library(cowplot)

custom_forestplot <- function(method_names, sample_size, point_estimates, variances, p_values, exposure, outcome) {
 
  if (!(length(method_names) == length(sample_size) &&
        length(sample_size) == length(point_estimates) &&
        length(point_estimates) == length(variances) &&
        length(variances) == length(p_values))) {
    stop("All input vectors must be of the same length")
  }
  
  
  z_value <- qnorm(0.975)  
  lower_ci <- point_estimates - z_value * variances
  upper_ci <- point_estimates + z_value * variances
  
  
  coef_ci <- c("Mean difference (95% CI)", paste(format(round(point_estimates, 2), nsmall = 2), " (",
                                                 format(round(lower_ci, 2), nsmall = 2), " ~ ",
                                                 format(round(upper_ci, 2), nsmall = 2), ")", sep = ""))
  N <- c("N", sample_size)
  pval <- c("P-value", formatC(p_values, format = "f", digits = 3))
  
  
  all_results <- tibble(
    mean = c(NA, point_estimates),
    lower = c(NA, lower_ci),
    upper = c(NA, upper_ci)
  )
  
 
  method <- c(paste("Effect of", exposure[1], "on", outcome[1]), method_names)
  tabletext <- cbind(
    method,
    N,
    coef_ci,
    pval
  )
  
  A <- forestplot(
    tabletext,
    all_results,
    graph.pos = 4, 
    new_page = TRUE,
    hrzl_lines = gpar(col="#000000"),
    is.summary=c(TRUE,rep(FALSE,7)),
    xlog = FALSE,
    col = fpColors(box = "black", line = "red"),
    boxsize = 0.1, 
    xlab = "Mean difference (95% CI)",
    txt_gp = fpTxtGp(
      label = gpar(fontfamily = "sans", fontsize = 6, cex = 1.2),  
      ticks = gpar(fontfamily = "sans", fontsize = 6),
      xlab = gpar(fontfamily = "sans", fontsize = 6),
    ),
    xticks = seq(-0.4,0.4,by = 0.1)
  )
  return(A)
  
}





custom_combined_forestplot <- function(A, B, C, group_name) {
  if (dev.cur() > 1) dev.off()
  
  A_plot <- custom_forestplot(A$method, A$nsnp, A$b, A$se, A$pval, "BMI", group_name)
  if (dev.cur() > 1) dev.off()
  B_plot <- custom_forestplot(B$method, B$nsnp, B$b, B$se, B$pval, "BMI", group_name)
  if (dev.cur() > 1) dev.off()
  C_plot <- custom_forestplot(C$method, C$nsnp, C$b, C$se, C$pval, "BMI", group_name)
  if (dev.cur() > 1) dev.off()
  A_grob <- grid.grabExpr(print(A_plot))
  B_grob <- grid.grabExpr(print(B_plot))
  C_grob <- grid.grabExpr(print(C_plot))
  
  combined_plot <- arrangeGrob(A_grob, B_grob, C_grob, ncol = 1)
  
  output_filename <- paste0("combined_forestplot_", group_name, ".png")
  ggsave(output_filename, plot = combined_plot, width = 15, height = 18, dpi = 600)
  return(combined_plot)
}



library(data.table)
BMI_lm        <- fread("BMI_lm.csv")
BMI_lm_parent <- fread("BMI_lm_parent.csv")
BMI_SEM       <- fread("BMI_base.csv")



BMI_lm$GWAS_method        <- "IVW lm"
BMI_lm_parent$GWAS_method <- "IVW lm_parent"
BMI_SEM$GWAS_method       <- "IVW FT-SEM"
BMI_lm$sample.size        <- "683,458"
BMI_lm_parent$sample.size <- "101,185"
BMI_SEM$sample.size       <- "101,185"





IVW_result <- rbind(BMI_SEM[BMI_SEM$method == "Inverse variance weighted",],
                    BMI_SEM[BMI_lm_parent$method == "Inverse variance weighted",],
                    BMI_lm[BMI_lm$method == "Inverse variance weighted",])

IVW_result_1 <- rbind(BMI_SEM[BMI_SEM$method == "Inverse variance weighted",][1,],
                      BMI_lm_parent[BMI_lm_parent$method == "Inverse variance weighted",][1,],
                      BMI_lm[BMI_lm$method == "Inverse variance weighted",][1,])

IVW_result_2 <- rbind(BMI_SEM[BMI_SEM$method == "Inverse variance weighted",][2,],
                      BMI_lm_parent[BMI_lm_parent$method == "Inverse variance weighted",][2,],
                      BMI_lm[BMI_lm$method == "Inverse variance weighted",][2,])

IVW_result_3 <- rbind(BMI_SEM[BMI_SEM$method == "Inverse variance weighted",][3,],
                      BMI_lm_parent[BMI_lm_parent$method == "Inverse variance weighted",][3,],
                      BMI_lm[BMI_lm$method == "Inverse variance weighted",][3,])

IVW_result_4 <- rbind(BMI_SEM[BMI_SEM$method == "Inverse variance weighted",][4,],
                      BMI_lm_parent[BMI_lm_parent$method == "Inverse variance weighted",][4,],
                      BMI_lm[BMI_lm$method == "Inverse variance weighted",][4,])

IVW_result_5 <- rbind(BMI_SEM[BMI_SEM$method == "Inverse variance weighted",][5,],
                      BMI_lm_parent[BMI_lm_parent$method == "Inverse variance weighted",][5,],
                      BMI_lm[BMI_lm$method == "Inverse variance weighted",][5,])

IVW_result_6 <- rbind(BMI_SEM[BMI_SEM$method == "Inverse variance weighted",][6,],
                      BMI_lm_parent[BMI_lm_parent$method == "Inverse variance weighted",][6,],
                      BMI_lm[BMI_lm$method == "Inverse variance weighted",][6,])



IVW_result_1$outcome <- "NIC "
IVW_result_2$outcome <- "CON "
IVW_result_3$outcome <- "DEP "
IVW_result_4$outcome <- "DRG "
IVW_result_5$outcome <- "BD  "




method_names    <- IVW_result_1$GWAS_method
sample_size     <- IVW_result_1$sample.size
point_estimates <- IVW_result_1$b
variances       <- IVW_result_1$se
p_values        <- IVW_result_1$pval
exposure        <- IVW_result_1$exposure
outcome         <- IVW_result_1$outcome




NIC <- custom_forestplot(IVW_result_1$GWAS_method,
                         IVW_result_1$sample.size,
                         IVW_result_1$b,
                         IVW_result_1$se,
                         IVW_result_1$pval,
                         IVW_result_1$exposure,
                         IVW_result_1$outcome)
CON <- custom_forestplot(IVW_result_2$GWAS_method,
                         IVW_result_2$sample.size,
                         IVW_result_2$b,
                         IVW_result_2$se,
                         IVW_result_2$pval,
                         IVW_result_2$exposure,
                         IVW_result_2$outcome)
DEP <- custom_forestplot(IVW_result_3$GWAS_method,
                         IVW_result_3$sample.size,
                         IVW_result_3$b,
                         IVW_result_3$se,
                         IVW_result_3$pval,
                         IVW_result_3$exposure,
                         IVW_result_3$outcome)
DRG <- custom_forestplot(IVW_result_4$GWAS_method,
                         IVW_result_4$sample.size,
                         IVW_result_4$b,
                         IVW_result_4$se,
                         IVW_result_4$pval,
                         IVW_result_4$exposure,
                         IVW_result_4$outcome)
BD <- custom_forestplot(IVW_result_5$GWAS_method,
                        IVW_result_5$sample.size,
                        IVW_result_5$b,
                        IVW_result_5$se,
                        IVW_result_5$pval,
                        IVW_result_5$exposure,
                        IVW_result_5$outcome)
FSIQ <- custom_forestplot(IVW_result_6$GWAS_method,
                          IVW_result_6$sample.size,
                          IVW_result_6$b,
                          IVW_result_6$se,
                          IVW_result_6$pval,
                          IVW_result_6$exposure,
                          IVW_result_6$outcome)






library(grid)
library(gridExtra)


NIC_grob <- grid.grabExpr(print(NIC))
CON_grob <- grid.grabExpr(print(CON))
DEP_grob <- grid.grabExpr(print(DEP))
DRG_grob <- grid.grabExpr(print(DRG))
BD_grob  <- grid.grabExpr(print(BD))
FSIQ_grob <- grid.grabExpr(print(FSIQ))

combined_plot <- plot_grid(NIC_grob, CON_grob, DEP_grob, DRG_grob, BD_grob, 
                           ncol = 1, labels = c("a", "b", "c", "d", "e"),
                           label_x = 0.01,
                           label_y = 1.05,
                           hjust = 0.5,
                           label_size = 14)  # label_x 控制标签的水平位置，0 为最左，1 为最右

ggsave(filename = "./TwoSampleMR.pdf",height = 180, width = 200, units = "mm")
ggsave(filename = "./TwoSampleMR.png",height = 180, width = 200, units = "mm", dpi = 300)




library(tidyverse)
library(forestplot)
library(cowplot)

custom_forestplot <- function(method_names, sample_size, point_estimates, variances, p_values, exposure, outcome) {
  
  if (!(length(method_names) == length(sample_size) &&
        length(sample_size) == length(point_estimates) &&
        length(point_estimates) == length(variances) &&
        length(variances) == length(p_values))) {
    stop("All input vectors must be of the same length")
  }
  
  
  z_value <- qnorm(0.975)  
  lower_ci <- point_estimates - z_value * variances
  upper_ci <- point_estimates + z_value * variances
  
 
  coef_ci <- c("Mean difference (95% CI)", paste(format(round(point_estimates, 2), nsmall = 2), " (",
                                                 format(round(lower_ci, 2), nsmall = 2), " ~ ",
                                                 format(round(upper_ci, 2), nsmall = 2), ")", sep = ""))
  N <- c("N", sample_size)
  pval <- c("P-value", formatC(p_values, format = "f", digits = 3))
  
  
  all_results <- tibble(
    mean = c(NA, point_estimates),
    lower = c(NA, lower_ci),
    upper = c(NA, upper_ci)
  )
  
 
  method <- c(paste("Effect of", exposure[1], "on", outcome[1]), method_names)
  tabletext <- cbind(
    method,
    coef_ci,
    pval
  )
  
  A <- forestplot(
    tabletext,
    all_results,
    graph.pos = 3, 
    new_page = TRUE,
    hrzl_lines = gpar(col="#000000"),
    is.summary=c(TRUE,rep(FALSE,7)),
    xlog = FALSE,
    col = fpColors(box = "black", line = "red"),
    boxsize = 0.1, 
    xlab = "Mean difference (95% CI)",
    txt_gp = fpTxtGp(
      label = gpar(fontfamily = "sans", fontsize = 6, cex = 1.2),  
      ticks = gpar(fontfamily = "sans", fontsize = 6),
      xlab = gpar(fontfamily = "sans", fontsize = 6),
    ),
    xticks = seq(-0.4,0.4,by = 0.1)
  )
  return(A)
  
}



BMI_SEM[BMI_SEM$method == "Inverse variance weighted",]$method <- "IVW"
BMI_lm_parent[BMI_lm_parent$method == "Inverse variance weighted",]$method <- "IVW"
BMI_lm[BMI_lm$method == "Inverse variance weighted",]$method <- "IVW"

BMI_SEM[BMI_SEM$method == "Weighted median",]$method <- "Weighted median (FT-SEM)"
BMI_lm_parent[BMI_lm_parent$method == "Weighted median",]$method <- "Weighted median (FT-SEM)"
BMI_lm[BMI_lm$method == "Weighted median",]$method <- "Weighted median (FT-SEM)"

BMI_SEM[BMI_SEM$method == "MR Egger",]$method <- "MR-Egger (FT-SEM)"
BMI_lm_parent[BMI_lm_parent$method == "MR Egger",]$method <- "MR-Egger (FT-SEM)"
BMI_lm[BMI_lm$method == "MR Egger",]$method <- "MR-Egger (FT-SEM)"

BMI_SEM[BMI_SEM$method == "Weighted mode",]$method <- "Weighted mode (FT-SEM)"
BMI_lm_parent[BMI_lm_parent$method == "Weighted mode",]$method <- "Weighted mode (FT-SEM)"
BMI_lm[BMI_lm$method == "Weighted mode",]$method <- "Weighted mode (FT-SEM)"




result_1 <- BMI_SEM[BMI_SEM$method != "Simple mode",][c(2,4,1),]
result_2 <- BMI_SEM[BMI_SEM$method != "Simple mode",][4+c(2,4,1),]
result_3 <- BMI_SEM[BMI_SEM$method != "Simple mode",][8+c(2,4,1),]
result_4 <- BMI_SEM[BMI_SEM$method != "Simple mode",][12+c(2,4,1),]
result_5 <- BMI_SEM[BMI_SEM$method != "Simple mode",][16+c(2,4,1),]
result_6 <- BMI_SEM[BMI_SEM$method != "Simple mode",][20+c(2,4,1),]





NIC <- custom_forestplot(result_1$method,
                         result_1$sample.size,
                         result_1$b,
                         result_1$se,
                         result_1$pval,
                         result_1$exposure,
                         result_1$outcome)
CON <- custom_forestplot(result_2$method,
                         result_2$sample.size,
                         result_2$b,
                         result_2$se,
                         result_2$pval,
                         result_2$exposure,
                         result_2$outcome)
DEP <- custom_forestplot(result_3$method,
                         result_3$sample.size,
                         result_3$b,
                         result_3$se,
                         result_3$pval,
                         result_3$exposure,
                         result_3$outcome)
DRG <- custom_forestplot(result_4$method,
                         result_4$sample.size,
                         result_4$b,
                         result_4$se,
                         result_4$pval,
                         result_4$exposure,
                         result_4$outcome)
BD <- custom_forestplot(result_5$method,
                        result_5$sample.size,
                        result_5$b,
                        result_5$se,
                        result_5$pval,
                        result_5$exposure,
                        result_5$outcome)
FSIQ <- custom_forestplot(result_6$method,
                          result_6$sample.size,
                          result_6$b,
                          result_6$se,
                          result_6$pval,
                          result_6$exposure,
                          result_6$outcome)




library(cowplot)

NIC_grob <- grid.grabExpr(print(NIC))
CON_grob <- grid.grabExpr(print(CON))
DEP_grob <- grid.grabExpr(print(DEP))
DRG_grob <- grid.grabExpr(print(DRG))
BD_grob  <- grid.grabExpr(print(BD))
FSIQ_grob <- grid.grabExpr(print(FSIQ))

combined_plot <- plot_grid(NIC_grob, CON_grob, DEP_grob, DRG_grob, BD_grob, 
                           ncol = 1, labels = c("a", "b", "c", "d", "e"),
                           label_x = 0.01,
                           label_y = 1.05,
                           hjust = 0.5,
                           label_size = 14)  

ggsave(filename = "./FT-SEM.pdf",height = 180, width = 200, units = "mm")
ggsave(filename = "./FT-SEM.png",height = 180, width = 200, units = "mm", dpi = 300)



BMI_SEM[BMI_SEM$method == "Weighted median (FT-SEM)",]$method <- "Weighted median (lm_parent)"
BMI_lm_parent[BMI_lm_parent$method == "Weighted median (FT-SEM)",]$method <- "Weighted median (lm_parent)"
BMI_lm[BMI_lm$method == "Weighted median (FT-SEM)",]$method <- "Weighted median (lm_parent)"

BMI_SEM[BMI_SEM$method == "MR-Egger (FT-SEM)",]$method <- "MR-Egger (lm_parent)"
BMI_lm_parent[BMI_lm_parent$method == "MR-Egger (FT-SEM)",]$method <- "MR-Egger (lm_parent)"
BMI_lm[BMI_lm$method == "MR-Egger (FT-SEM)",]$method <- "MR-Egger (lm_parent)"

BMI_SEM[BMI_SEM$method == "Weighted mode (FT-SEM)",]$method <- "Weighted mode (lm_parent)"
BMI_lm_parent[BMI_lm_parent$method == "Weighted mode (FT-SEM)",]$method <- "Weighted mode (lm_parent)"
BMI_lm[BMI_lm$method == "Weighted mode (FT-SEM)",]$method <- "Weighted mode (lm_parent)"




result_1 <- BMI_lm_parent[BMI_lm_parent$method != "Simple mode",][c(2,4,1),]
result_2 <- BMI_lm_parent[BMI_lm_parent$method != "Simple mode",][4+c(2,4,1),]
result_3 <- BMI_lm_parent[BMI_lm_parent$method != "Simple mode",][8+c(2,4,1),]
result_4 <- BMI_lm_parent[BMI_lm_parent$method != "Simple mode",][12+c(2,4,1),]
result_5 <- BMI_lm_parent[BMI_lm_parent$method != "Simple mode",][16+c(2,4,1),]
result_6 <- BMI_lm_parent[BMI_lm_parent$method != "Simple mode",][20+c(2,4,1),]



NIC <- custom_forestplot(result_1$method,
                         result_1$sample.size,
                         result_1$b,
                         result_1$se,
                         result_1$pval,
                         result_1$exposure,
                         result_1$outcome)
CON <- custom_forestplot(result_2$method,
                         result_2$sample.size,
                         result_2$b,
                         result_2$se,
                         result_2$pval,
                         result_2$exposure,
                         result_2$outcome)
DEP <- custom_forestplot(result_3$method,
                         result_3$sample.size,
                         result_3$b,
                         result_3$se,
                         result_3$pval,
                         result_3$exposure,
                         result_3$outcome)
DRG <- custom_forestplot(result_4$method,
                         result_4$sample.size,
                         result_4$b,
                         result_4$se,
                         result_4$pval,
                         result_4$exposure,
                         result_4$outcome)
BD <- custom_forestplot(result_5$method,
                        result_5$sample.size,
                        result_5$b,
                        result_5$se,
                        result_5$pval,
                        result_5$exposure,
                        result_5$outcome)
FSIQ <- custom_forestplot(result_6$method,
                          result_6$sample.size,
                          result_6$b,
                          result_6$se,
                          result_6$pval,
                          result_6$exposure,
                          result_6$outcome)




library(cowplot)

NIC_grob <- grid.grabExpr(print(NIC))
CON_grob <- grid.grabExpr(print(CON))
DEP_grob <- grid.grabExpr(print(DEP))
DRG_grob <- grid.grabExpr(print(DRG))
BD_grob  <- grid.grabExpr(print(BD))
FSIQ_grob <- grid.grabExpr(print(FSIQ))

combined_plot <- plot_grid(NIC_grob, CON_grob, DEP_grob, DRG_grob, BD_grob, 
                           ncol = 1, labels = c("a", "b", "c", "d", "e"),
                           label_x = 0.01,
                           label_y = 1.05,
                           hjust = 0.5,
                           label_size = 14)  # label_x 控制标签的水平位置，0 为最左，1 为最右

ggsave(filename = "./lm_parent.pdf",height = 180, width = 200, units = "mm")
ggsave(filename = "./lm_parent.png",height = 180, width = 200, units = "mm", dpi = 300)



BMI_SEM[BMI_SEM$method == "Weighted median (lm_parent)",]$method <- "Weighted median (lm)"
BMI_lm_parent[BMI_lm_parent$method == "Weighted median (lm_parent)",]$method <- "Weighted median (lm)"
BMI_lm[BMI_lm$method == "Weighted median (lm_parent)",]$method <- "Weighted median (lm)"

BMI_SEM[BMI_SEM$method == "MR-Egger (lm_parent)",]$method <- "MR-Egger (lm)"
BMI_lm_parent[BMI_lm_parent$method == "MR-Egger (lm_parent)",]$method <- "MR-Egger (lm)"
BMI_lm[BMI_lm$method == "MR-Egger (lm_parent)",]$method <- "MR-Egger (lm)"

BMI_SEM[BMI_SEM$method == "Weighted mode (lm_parent)",]$method <- "Weighted mode (lm)"
BMI_lm_parent[BMI_lm_parent$method == "Weighted mode (lm_parent)",]$method <- "Weighted mode (lm)"
BMI_lm[BMI_lm$method == "Weighted mode (lm_parent)",]$method <- "Weighted mode (lm)"




result_1 <- BMI_lm[BMI_lm$method != "Simple mode",][c(2,4,1),]
result_2 <- BMI_lm[BMI_lm$method != "Simple mode",][4+c(2,4,1),]
result_3 <- BMI_lm[BMI_lm$method != "Simple mode",][8+c(2,4,1),]
result_4 <- BMI_lm[BMI_lm$method != "Simple mode",][12+c(2,4,1),]
result_5 <- BMI_lm[BMI_lm$method != "Simple mode",][16+c(2,4,1),]
result_6 <- BMI_lm[BMI_lm$method != "Simple mode",][20+c(2,4,1),]



NIC <- custom_forestplot(result_1$method,
                         result_1$sample.size,
                         result_1$b,
                         result_1$se,
                         result_1$pval,
                         result_1$exposure,
                         result_1$outcome)
CON <- custom_forestplot(result_2$method,
                         result_2$sample.size,
                         result_2$b,
                         result_2$se,
                         result_2$pval,
                         result_2$exposure,
                         result_2$outcome)
DEP <- custom_forestplot(result_3$method,
                         result_3$sample.size,
                         result_3$b,
                         result_3$se,
                         result_3$pval,
                         result_3$exposure,
                         result_3$outcome)
DRG <- custom_forestplot(result_4$method,
                         result_4$sample.size,
                         result_4$b,
                         result_4$se,
                         result_4$pval,
                         result_4$exposure,
                         result_4$outcome)
BD <- custom_forestplot(result_5$method,
                        result_5$sample.size,
                        result_5$b,
                        result_5$se,
                        result_5$pval,
                        result_5$exposure,
                        result_5$outcome)
FSIQ <- custom_forestplot(result_6$method,
                          result_6$sample.size,
                          result_6$b,
                          result_6$se,
                          result_6$pval,
                          result_6$exposure,
                          result_6$outcome)




library(cowplot)

NIC_grob <- grid.grabExpr(print(NIC))
CON_grob <- grid.grabExpr(print(CON))
DEP_grob <- grid.grabExpr(print(DEP))
DRG_grob <- grid.grabExpr(print(DRG))
BD_grob  <- grid.grabExpr(print(BD))
FSIQ_grob <- grid.grabExpr(print(FSIQ))

combined_plot <- plot_grid(NIC_grob, CON_grob, DEP_grob, DRG_grob, BD_grob, 
                           ncol = 1, labels = c("a", "b", "c", "d", "e"),
                           label_x = 0.01,
                           label_y = 1.05,
                           hjust = 0.5,
                           label_size = 14)  

ggsave(filename = "./lm.pdf",height = 180, width = 200, units = "mm")
ggsave(filename = "./lm.png",height = 180, width = 200, units = "mm", dpi = 300)
