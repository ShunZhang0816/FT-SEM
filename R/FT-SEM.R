#' @title The main function for family trios structural equation model (FT-SEM)
#' @description A robust and powerful GWAS method for family trios
#' @param data_sub A data.frame object for FT-SEM contains six columns: the father's genotype, the mother's genotype, the offspring's genotype, the father's phenotype, the mother's phenotype, and the offspring's phenotype. These columns are named f_snp, m_snp, o_snp, Exposure_f, Exposure_m, and Exposure_o, respectively.
#' @import OpenMx
#' @export

FT_SEM <- function(data_sub = NULL, snp_name = NA){
  # Set up model parameters
  manifests <- names(data_sub)
  latents <- c("e_f","e_m","e_o","var_o_snp","var_f_snp","var_m_snp","gf_snp_1","gm_snp_1","gf_snp_2","gm_snp_2")
  IVmodel_base <- mxModel(model = "IV model", type = "RAM", mxData(observed = data_sub, type = "raw"),
                          mxPath(from = "gf_snp_1", arrows = 2, free = TRUE, values = 1, labels = "theta_snp"),
                          mxPath(from = "gm_snp_1", arrows = 2, free = TRUE, values = 1, labels = "theta_snp"),
                          mxPath(from = "gf_snp_2", arrows = 2, free = TRUE, values = 1, labels = "theta_snp"),
                          mxPath(from = "gm_snp_2", arrows = 2, free = TRUE, values = 1, labels = "theta_snp"),
                          mxPath(from = "var_f_snp", arrows = 2, free = TRUE, values = 1, labels = "theta_snp"),
                          mxPath(from = "var_m_snp", arrows = 2, free = TRUE, values = 1, labels = "theta_snp"),
                          mxPath(from = "var_o_snp", arrows = 2, free = TRUE, values = 1, labels = "theta_snp"),
                          
                          mxPath(from = "gf_snp_1", to = "f_snp", arrows = 1, free = FALSE,
                                 values = 0.5, labels = "snp_gf1_to_f"),
                          mxPath(from = "gm_snp_1", to = "f_snp", arrows = 1, free = FALSE,
                                 values = 0.5, labels = "snp_gm1_to_f"),
                          mxPath(from = "gf_snp_2", to = "m_snp", arrows = 1, free = FALSE,
                                 values = 0.5, labels = "snp_gf2_to_m"),
                          mxPath(from = "gm_snp_2", to = "m_snp", arrows = 1, free = FALSE,
                                 values = 0.5, labels = "snp_gm2_to_m"),
                          mxPath(from = "f_snp", to = "o_snp", arrows = 1, free = FALSE,
                                 values = 0.5, labels = "snp_f_to_o"),
                          mxPath(from = "m_snp", to = "o_snp", arrows = 1, free = FALSE,
                                 values = 0.5, labels = "snp_m_to_o"),
                          mxPath(from = "var_f_snp", to = "f_snp", arrows = 1, free = FALSE,
                                 values = sqrt(0.5), labels = "var_f_f"),
                          mxPath(from = "var_m_snp", to = "m_snp", arrows = 1, free = FALSE,
                                 values = sqrt(0.5), labels = "var_m_m"),
                          mxPath(from = "var_o_snp", to = "o_snp", arrows = 1, free = FALSE,
                                 values = sqrt(0.5), labels = "var_o_o"),
                          
                          mxPath(from="e_f", to="e_o", arrows=2, free=TRUE, values=0.2, labels=c("ρ")),
                          mxPath(from="e_m", to="e_o", arrows=2, free=TRUE, values=0.2, labels=c("ρ")),
                          mxPath(from="e_f", to="e_m", arrows=2, free=TRUE, values=0.2, labels=c("ρ")),
                          
                          mxPath(from = "e_f", arrows = 2, free = TRUE, values = 1, labels = "var_e_f"),
                          mxPath(from = "e_m", arrows = 2, free = TRUE, values = 1, labels = "var_e_m"),
                          mxPath(from = "e_o", arrows = 2, free = TRUE, values = 1, labels = "var_e_o"),
                          mxPath(from = "e_f", to = "Exposure_f", free = FALSE, values = 1,
                                 labels = "effect_e_f_E"),
                          mxPath(from = "e_m", to = "Exposure_m", free = FALSE, values = 1,
                                 labels = "effect_e_m_E"),
                          mxPath(from = "e_o", to = "Exposure_o", free = FALSE, values = 1,
                                 labels = "effect_e_o_E"),
                          
                          mxPath(from = "f_snp", to = "Exposure_o", free = TRUE,
                                 arrows = 1, labels = "effect_f_e"),
                          mxPath(from = "m_snp", to = "Exposure_o", free = TRUE,
                                 arrows = 1, labels = "effect_m_e"),
                          mxPath(from = "o_snp", to = "Exposure_o", free = TRUE,
                                 arrows = 1, labels = "effect_o_e"),
                          mxPath(from = "gf_snp_1", to = "Exposure_f", free = TRUE,
                                 arrows = 1, labels = "effect_f_e"),
                          mxPath(from = "gm_snp_1", to = "Exposure_f", free = TRUE,
                                 arrows = 1, labels = "effect_m_e"),
                          mxPath(from = "f_snp", to = "Exposure_f", free = TRUE,
                                 arrows = 1, labels = "effect_o_e"),
                          mxPath(from = "gf_snp_2", to = "Exposure_m", free = TRUE,
                                 arrows = 1, labels = "effect_f_e"),
                          mxPath(from = "gm_snp_2", to = "Exposure_m", free = TRUE,
                                 arrows = 1, labels = "effect_m_e"),
                          mxPath(from = "m_snp", to = "Exposure_m", free = TRUE,
                                 arrows = 1, labels = "effect_o_e"),
                          mxPath( from="one", to=manifests, arrows=1, free= c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
                                  values=rep(0,length(manifests)), labels=c("mean_snp_f","mean_snp_m","mean_snp_o","mean_ef","mean_em","mean_eo")),  
                          manifestVars = manifests,latentVars = latents
  )
  
  
  IVFit_base  <- mxRun(IVmodel_base)
  Zscore_base <- summary(IVFit_base)$parameters[,5]/summary(IVFit_base)$parameters[,6]
  Pval_base   <- 2*(1-pnorm(abs(Zscore_base),0,1))
  
  base_summary            <- summary(IVFit_base)
  
  
  base <- rbind(c(snp_name,"FT-SEM",
                  base_summary$parameters[1,5],
                  base_summary$parameters[1,6],Pval_base[1],
                  (base_summary$parameters[1,5]-1.96*base_summary$parameters[1,6]),
                  (base_summary$parameters[1,5]+1.96*base_summary$parameters[1,6]),
                  base_summary$parameters[2,5],
                  base_summary$parameters[2,6],Pval_base[2],
                  (base_summary$parameters[2,5]-1.96*base_summary$parameters[2,6]),
                  (base_summary$parameters[2,5]+1.96*base_summary$parameters[2,6]),
                  base_summary$parameters[3,5],
                  base_summary$parameters[3,6],Pval_base[3],
                  (base_summary$parameters[3,5]-1.96*base_summary$parameters[3,6]),
                  (base_summary$parameters[3,5]+1.96*base_summary$parameters[3,6])
                  )
  )
  
  base <- as.data.frame(base)
  names(base) <- c("SNP","Method",
                   "Beta_o_e", "SE_o_e", "P_wald_o_e","CI_lower_o_e","CI_upper_o_e",
                   "Beta_f_e", "SE_f_e", "P_wald_f_e","CI_lower_f_e","CI_upper_f_e",
                   "Beta_m_e", "SE_m_e", "P_wald_m_e","CI_lower_m_e","CI_upper_m_e"
  )
  return(base)
}

