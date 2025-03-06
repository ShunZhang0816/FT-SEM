library(foreach)
library(doParallel)
library(data.table)
library(dplyr)
library(OpenMx)
library(RNOmni)
library(qqman)
library(biomaRt)
library(tidyverse)

# Upload data before run the data
# Genetic data as "mydata"
# Phenotype data as "mydata_phe"
# Snp information as "snp_info
# PCA results as "eigvec"

names(snp_info) <- c("chr","snp_name","location1","location2","minor_allele","major_allele")
names(eigvec) = c("Fid","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
eigvec        <- eigvec[,2:ncol(eigvec)]


bind_data_1      <- cbind(merge(mydata[,c(1:6,7),with = F],mydata_phe,by.x = "IID", by.y = "SUBJECT_ID"))
data_family      <- left_join(bind_data_1,eigvec,by = c("IID"="IID"))
data_independent <- left_join(bind_data_1,eigvec,by = c("IID"="IID"))

snp_name      <- names(data_family[,7])

children <- data_family[!(data_family$PAT == 0 & data_family$MAT == 0), ]
children_name <- names(children)

parents <- data_family[data_family$IID %in% unique(c(children$PAT, children$MAT)), ]

names(children) <- paste(names(children),"_o",sep = "")


children_with_fathers <- inner_join(children, parents, by = c("PAT_o" = "IID"))
names(children_with_fathers) <- c(names(children),paste(children_name,"_f",sep = ""))[-(ncol(children)+1)]
family_data <- inner_join(children_with_fathers, parents, by = c("MAT_o" = "IID"))
names(family_data) <- c(names(children_with_fathers),paste(children_name,"_m",sep = ""))[-(ncol(children_with_fathers)+1)]


set.seed(2024)
family_unique <- family_data %>% group_by(FID_o) %>% sample_n(1)


independent_data <- data_independent %>% group_by(FID) %>% sample_n(1)


selected_IID_family      <- family_unique[,1]
selected_IID_independent <- independent_data[,1]


mydata_name      <- names(mydata) 
nCircle          <- 5000 
nSum             <- ncol(mydata) 
# nSum             <- 1050     

result_name      <- c("NIC_parents","CON_parents","DEP_parents","DRG_parents","BD_parents","FSIQ")


cl <- makeCluster(32)
registerDoParallel(cl)

time1 <- Sys.time()
for( k in 1:6){
  
  base_f_result      <- data.frame()
  lm_parent_f_result <- data.frame()
  base_m_result      <- data.frame()
  lm_parent_m_result <- data.frame()
  base_o_result      <- data.frame()
  lm_parent_o_result <- data.frame()
  lm_o_result        <- data.frame()
  
  end_result         <- data.frame()
  result_2           <- data.frame()
  for( i in 1:(nSum%/%nCircle + 1)){
    if (i == 1){
      
      data_2 <- mydata[,c(1:6,((i-1)*nCircle + 7):(i*nCircle)),with = F]
      nsnp   <- nCircle
    } else if (i == ((nSum %/% nCircle) + 1)){
      
      data_2 <- mydata[,c(1:6,((i-1)*nCircle + 1):(nSum)),with = F]
      nsnp   <- nSum - ((i-1)*nCircle + 1) + 1 + 6
    } else {
      
      data_2 <- mydata[,c(1:6,((i-1)*nCircle + 1):(i*nCircle)),with = F]
      nsnp   <- nCircle + 6
    }
    
    time11 <- Sys.time()
    
    result <- foreach (j = 7:nsnp, .combine = rbind,
                       .packages = c("data.table","dplyr","OpenMx","RNOmni"),
                       .errorhandling = "pass") %dopar%{
                         bind_data_1      <- cbind(merge(data_2[,c(1:6,j),with = F],mydata_phe,by.x = "IID", by.y = "SUBJECT_ID"))
                         data_family      <- left_join(bind_data_1,eigvec,by = c("IID"="IID"))
                         data_independent <- left_join(bind_data_1,eigvec,by = c("IID"="IID"))
                         
                         snp_name      <- names(data_family[,7])
                         
                         children <- data_family[!(data_family$PAT == 0 & data_family$MAT == 0), ]
                         children_name <- names(children)
                         
                         parents <- data_family[data_family$IID %in% unique(c(children$PAT, children$MAT)), ]
                         
                         names(children) <- paste(names(children),"_o",sep = "")
                         
                        
                         children_with_fathers <- inner_join(children, parents, by = c("PAT_o" = "IID"))
                         names(children_with_fathers) <- c(names(children),paste(children_name,"_f",sep = ""))[-(ncol(children)+1)]
                         family_data <- inner_join(children_with_fathers, parents, by = c("MAT_o" = "IID"))
                         names(family_data) <- c(names(children_with_fathers),paste(children_name,"_m",sep = ""))[-(ncol(children_with_fathers)+1)]
                         
                         
                         family_unique    <- inner_join(family_data,selected_IID_family,by = c("IID_o" = "IID_o"))
                         independent_data <- inner_join(data_independent,selected_IID_independent,by = c("IID" = "IID"))
                         
                         # data_sub          <- family_unique[,c(7,10:13,13 + k,19:28,
                         #                                       34,37:40,40 + k, 
                         #                                       61,64:67,67 + k 
                         #                                       ),with = F]
                         
                         
                         data_sub          <- family_unique[,c(7,10:13,13 + k,20:29,
                                                               35,38:41,41 + k, 
                                                               63,66:69,69 + k 
                         ),with = F]
                         
                         data_sub_2        <- independent_data[,c(7,10:13,13 + k,20:29),with = F]
                         names(data_sub)   <- c("o_snp","SEX_o","AGE_o","Birthyear_o","Generation_o","Exposure_o",
                                                paste0("PC", 1:10,"_o"),
                                                "f_snp","SEX_f","AGE_f","Birthyear_f","Generation_f","Exposure_f",
                                                "m_snp","SEX_m","AGE_m","Birthyear_m","Generation_m","Exposure_m"
                         )
                         names(data_sub_2) <- c("o_snp","SEX_o","AGE_o","Birthyear_o","Generation_o","Exposure_o",
                                                paste0("PC", 1:10,"_o"))
                         
                         data_sub          <- sapply(data_sub, as.numeric)
                         data_sub_2        <- sapply(data_sub_2, as.numeric)
                         data_sub          <- as.data.frame(data_sub)
                         data_sub_2        <- as.data.frame(data_sub_2)
                         
                         
                         data_sub   <- na.omit(data_sub)
                         data_sub_2 <- na.omit(data_sub_2)
                         
                         if (nrow(data_sub) == 0 || all(is.na(data_sub))) {
                           
                           empty_row <- matrix(NA,1,35)
                           return(empty_row)
                         }
                         
                         
                         custom_normalize <- function(x) {
                           (x - mean(x)) / sd(x)
                         }
                         data_sub          <- as.data.frame(lapply(data_sub, RankNorm))
                         # normalized_data   <- as.data.frame(lapply(data_sub, custom_normalize))
                         data_sub_2        <- as.data.frame(lapply(data_sub_2, RankNorm))
                         # normalized_data_2 <- as.data.frame(lapply(data_sub_2, custom_normalize))
                         # Set up model parameters
                         
                         manifests <- c("o_snp","SEX_o","AGE_o","Birthyear_o","Exposure_o",
                                        "f_snp","AGE_f","Birthyear_f","Exposure_f",
                                        "m_snp","AGE_m","Birthyear_m","Exposure_m"
                         )
                         latents <- c("e_f","e_m","e_o","var_o_snp","var_f_snp","var_m_snp",
                                      "gf_snp_1","gm_snp_1","gf_snp_2","gm_snp_2",
                                      "var_coeffect_o_e_1",
                                      "var_coeffect_o_e_2",
                                      "var_coeffect_o_e_3")
                         
                         IVmodel_base <- mxModel(model = "IV model", type = "RAM", mxData(observed = data_sub, type = "raw"),
                                                 
                                                 mxPath(from = "gf_snp_1", arrows = 2, free = TRUE, values = 0.2, labels = "theta_snp"),
                                                 mxPath(from = "gm_snp_1", arrows = 2, free = TRUE, values = 0.2, labels = "theta_snp"),
                                                 mxPath(from = "gf_snp_2", arrows = 2, free = TRUE, values = 0.2, labels = "theta_snp"),
                                                 mxPath(from = "gm_snp_2", arrows = 2, free = TRUE, values = 0.2, labels = "theta_snp"),
                                                 mxPath(from = "var_f_snp", arrows = 2, free = TRUE, values = 0.2, labels = "theta_snp"),
                                                 mxPath(from = "var_m_snp", arrows = 2, free = TRUE, values = 0.2, labels = "theta_snp"),
                                                 mxPath(from = "var_o_snp", arrows = 2, free = TRUE, values = 0.2, labels = "theta_snp"),
                                                 
                                                 mxPath(from = "gf_snp_1", to = "f_snp", arrows = 1, free = FALSE, values = 0.5, labels = "snp_gf1_to_f"),
                                                 mxPath(from = "gm_snp_1", to = "f_snp", arrows = 1, free = FALSE, values = 0.5, labels = "snp_gm1_to_f"),
                                                 mxPath(from = "gf_snp_2", to = "m_snp", arrows = 1, free = FALSE, values = 0.5, labels = "snp_gf2_to_m"),
                                                 mxPath(from = "gm_snp_2", to = "m_snp", arrows = 1, free = FALSE, values = 0.5, labels = "snp_gm2_to_m"),
                                                 mxPath(from = "f_snp", to = "o_snp", arrows = 1, free = FALSE, values = 0.5, labels = "snp_f_to_o"),
                                                 mxPath(from = "m_snp", to = "o_snp", arrows = 1, free = FALSE, values = 0.5, labels = "snp_m_to_o"),
                                                 mxPath(from = "var_f_snp", to = "f_snp", arrows = 1, free = FALSE, values = sqrt(0.5), labels = "var_f_f"),
                                                 mxPath(from = "var_m_snp", to = "m_snp", arrows = 1, free = FALSE, values = sqrt(0.5), labels = "var_m_m"),
                                                 mxPath(from = "var_o_snp", to = "o_snp", arrows = 1, free = FALSE, values = sqrt(0.5), labels = "var_o_o"),
                                                 
                                                 mxPath(from = "e_f", arrows = 2, free = TRUE, values = 1, labels = "var_e_f"),
                                                 mxPath(from = "e_m", arrows = 2, free = TRUE, values = 1, labels = "var_e_m"),
                                                 mxPath(from = "e_o", arrows = 2, free = TRUE, values = 1, labels = "var_e_o"),
                                                 mxPath(from = "e_f", to = "Exposure_f", free = FALSE, values = 1, labels = "effect_e_f_E"),
                                                 mxPath(from = "e_m", to = "Exposure_m", free = FALSE, values = 1, labels = "effect_e_m_E"),
                                                 mxPath(from = "e_o", to = "Exposure_o", free = FALSE, values = 1, labels = "effect_e_o_E"),
                                                 
                                                 mxPath(from="e_f", to="e_o", arrows=2, free=TRUE, values=0.2, labels=c("ρ")),
                                                 mxPath(from="e_m", to="e_o", arrows=2, free=TRUE, values=0.2, labels=c("ρ")),
                                                 mxPath(from="e_f", to="e_m", arrows=2, free=TRUE, values=0.2, labels=c("ρ")),
                                                 
                                                 mxPath(from = "f_snp", to = "Exposure_o", free = TRUE, arrows = 1, labels = "effect_f_e"),
                                                 mxPath(from = "m_snp", to = "Exposure_o", free = TRUE, arrows = 1, labels = "effect_m_e"),
                                                 mxPath(from = "o_snp", to = "Exposure_o", free = TRUE, arrows = 1, labels = "effect_o_e"),
                                                 mxPath(from = "gf_snp_1", to = "Exposure_f", free = TRUE, arrows = 1, labels = "effect_f_e"),
                                                 mxPath(from = "gm_snp_1", to = "Exposure_f", free = TRUE, arrows = 1, labels = "effect_m_e"),
                                                 mxPath(from = "f_snp", to = "Exposure_f", free = TRUE, arrows = 1, labels = "effect_o_e"),
                                                 mxPath(from = "gf_snp_2", to = "Exposure_m", free = TRUE, arrows = 1, labels = "effect_f_e"),
                                                 mxPath(from = "gm_snp_2", to = "Exposure_m", free = TRUE, arrows = 1, labels = "effect_m_e"),
                                                 mxPath(from = "m_snp", to = "Exposure_m", free = TRUE, arrows = 1, labels = "effect_o_e"),
                                                
                                                 mxPath(from = "SEX_o", to = "Exposure_o", free = TRUE, arrows = 1, labels = "coeffect_o_e_1"),
                                                 mxPath(from = "AGE_o", to = "Exposure_o", free = TRUE, arrows = 1, labels = "coeffect_o_e_2"),
                                                 mxPath(from = "AGE_f", to = "Exposure_f", free = TRUE, arrows = 1, labels = "coeffect_o_e_2"),
                                                 mxPath(from = "AGE_m", to = "Exposure_m", free = TRUE, arrows = 1, labels = "coeffect_o_e_2"),
                                                 mxPath(from = "Birthyear_o", to = "Exposure_o", free = TRUE, arrows = 1, labels = "coeffect_o_e_3"),
                                                 mxPath(from = "Birthyear_f", to = "Exposure_f", free = TRUE, arrows = 1, labels = "coeffect_o_e_3"),
                                                 mxPath(from = "Birthyear_m", to = "Exposure_m", free = TRUE, arrows = 1, labels = "coeffect_o_e_3"),
                                                 
                                                 mxPath(from = "SEX_o", free = TRUE, arrows = 2, values = 1, labels = "var_coeffect_o_e_1"),
                                                 mxPath(from = "AGE_o", free = TRUE, arrows = 2, values = 1, labels = "var_coeffect_o_e_2_o"),
                                                 mxPath(from = "AGE_f", free = TRUE, arrows = 2, values = 1, labels = "var_coeffect_o_e_2_f"),
                                                 mxPath(from = "AGE_m", free = TRUE, arrows = 2, values = 1, labels = "var_coeffect_o_e_2_m"),
                                                 mxPath(from = "Birthyear_o", free = TRUE, arrows = 2, values = 1, labels = "var_coeffect_o_e_3_o"),
                                                 mxPath(from = "Birthyear_f", free = TRUE, arrows = 2, values = 1, labels = "var_coeffect_o_e_3_f"),
                                                 mxPath(from = "Birthyear_m", free = TRUE, arrows = 2, values = 1, labels = "var_coeffect_o_e_3_m"),
                                                 mxPath( from="one", to=manifests, arrows=1, free= c(TRUE),
                                                         values=rep(0,length(manifests)),
                                                         labels=paste0("mean_", manifests)),
                                                 manifestVars = manifests,latentVars = latents
                         )
                         
                         IVFit_base <- mxRun(IVmodel_base)
                         Zscore_base <- summary(IVFit_base)$parameters[,5]/summary(IVFit_base)$parameters[,6]
                         Pval_base <- 2*(1-pnorm(abs(Zscore_base),0,1))
                         
                         IVFit_result     <- summary(IVFit_base)
                         
                         # Exposure_f_std   <- sd(data_sub$Exposure_o)
                         # f_snp_std        <- sd(data_sub$f_snp)
                         # Exposure_m_std   <- sd(data_sub$Exposure_o)
                         # m_snp_std        <- sd(data_sub$m_snp)
                         # Exposure_o_std   <- sd(data_sub$Exposure_o)
                         # o_snp_std        <- sd(data_sub$o_snp)
                         # 
                         # 
                         # tran_f_1         <- Exposure_f_std / f_snp_std
                         # tran_m_1         <- Exposure_m_std / m_snp_std
                         # tran_o_1         <- Exposure_o_std / o_snp_std
                         
                         lm_parent   <- lm(Exposure_o ~ o_snp + f_snp + m_snp + SEX_o +  Birthyear_o + AGE_o,
                                           data_sub)
                         lm_o        <- lm(Exposure_o ~ o_snp +
                                             SEX_o * Generation_o +
                                             Birthyear_o * Generation_o +
                                             AGE_o * Generation_o +
                                             PC1_o + PC2_o + PC3_o + PC4_o + PC5_o + PC6_o + PC7_o + PC8_o + PC9_o + PC10_o, data_sub_2)
                         
                         lm_parent_summary <- summary(lm_parent)
                         lm_summary        <- summary(lm_o)
                         
                         base_o_result <- rbind(c(
                           IVFit_result$parameters[1,5] ,
                           IVFit_result$parameters[1,6] ,
                           Pval_base[1],
                           (IVFit_result$parameters[1,5]-1.96*IVFit_result$parameters[1,6]),
                           (IVFit_result$parameters[1,5]+1.96*IVFit_result$parameters[1,6])
                         )
                         )
                         
                         
                         lm_parent_o_result  <- rbind(c(
                           lm_parent_summary$coefficients[2,1],
                           lm_parent_summary$coefficients[2,2],
                           lm_parent_summary$coefficients[2,4],
                           (lm_parent_summary$coefficients[2,1]-1.96*lm_parent_summary$coefficients[2,2]),
                           (lm_parent_summary$coefficients[2,1]+1.96*lm_parent_summary$coefficients[2,2])
                         )
                         )
                         
                         lm_o_result        <- rbind(c(
                           lm_summary$coefficients[2,1],
                           lm_summary$coefficients[2,2],
                           lm_summary$coefficients[2,4],
                           (lm_summary$coefficients[2,1]-1.96*lm_summary$coefficients[2,2]),
                           (lm_summary$coefficients[2,1]+1.96*lm_summary$coefficients[2,2])
                         )
                         )
                         
                         base_f_result <- rbind(c(
                           IVFit_result$parameters[5,5],
                           IVFit_result$parameters[5,6],
                           Pval_base[5],
                           (IVFit_result$parameters[5,5]-1.96*IVFit_result$parameters[5,6]),
                           (IVFit_result$parameters[5,5]+1.96*IVFit_result$parameters[5,6])
                         )
                         )
                         
                         
                         lm_parent_f_result  <- rbind(c(
                           lm_parent_summary$coefficients[3,1],
                           lm_parent_summary$coefficients[3,2],
                           lm_parent_summary$coefficients[3,4],
                           (lm_parent_summary$coefficients[3,1]-1.96*lm_parent_summary$coefficients[3,2]),
                           (lm_parent_summary$coefficients[3,1]+1.96*lm_parent_summary$coefficients[3,2])
                         )
                         )
                         
                         base_m_result <- rbind(c(
                           IVFit_result$parameters[6,5],
                           IVFit_result$parameters[6,6],
                           Pval_base[6],
                           (IVFit_result$parameters[6,5]-1.96*IVFit_result$parameters[6,6]),
                           (IVFit_result$parameters[6,5]+1.96*IVFit_result$parameters[6,6])
                         )
                         )
                         
                         
                         lm_parent_m_result  <- rbind(c(
                           lm_parent_summary$coefficients[4,1],
                           lm_parent_summary$coefficients[4,2],
                           lm_parent_summary$coefficients[4,4],
                           (lm_parent_summary$coefficients[4,1]-1.96*lm_parent_summary$coefficients[4,2]),
                           (lm_parent_summary$coefficients[4,1]+1.96*lm_parent_summary$coefficients[4,2])
                         )
                         )
                         
                         # names(base_result) <- c(
                         #                    "Beta_base", "SE_base", "P_wald_base","CI_lower_base","CI_upper_base"
                         #                    )
                         # names(lm_result) <- c(
                         #                    "Beta_lm", "SE_lm", "P_wald_lm","CI_lower_lm","CI_upper_lm"
                         #                    )
                         # names(lm_parent_result) <- c(
                         #                    "Beta_lm_parent","SE_lm_parent","P_wald_lm_parent","CI_lower_lm_parent","CI_upper_lm_parent"
                         #                    )
                         bind_result <- cbind(base_f_result,lm_parent_f_result,
                                              base_m_result,lm_parent_m_result,
                                              base_o_result,lm_parent_o_result,lm_o_result)
                         bind_result
                       }
    result_2 <- rbind(result_2,result)
    
    time22 <- Sys.time()
    print(i)
    print(time22-time11)
    gc() 
  }
  end_result <- as.data.frame(cbind(mydata_name[7:nSum],result_2))
  names(end_result) <- c("snp_name","Beta_f_base", "SE_f_base", "P_wald_f_base","CI_lower_base_f","CI_upper_base_f",
                         "Beta_lm_parent_f","SE_lm_parent_f","P_wald_lm_parent_f","CI_lower_lm_parent_f","CI_upper_lm_parent_f",
                         "Beta_m_base", "SE_m_base", "P_wald_m_base","CI_lower_base_m","CI_upper_base_m",
                         "Beta_lm_parent_m","SE_lm_parent_m","P_wald_lm_parent_m","CI_lower_lm_parent_m","CI_upper_lm_parent_m",
                         "Beta_o_base", "SE_o_base", "P_wald_o_base","CI_lower_base_o","CI_upper_base_o",
                         "Beta_lm_parent_o","SE_lm_parent_o","P_wald_lm_parent_o","CI_lower_lm_parent_o","CI_upper_lm_parent_o",
                         "Beta_lm_o","SE_lm_o","P_wald_lm_o","CI_lower_lm_o","CI_upper_lm_o"
  )
  
  # names(end_result) <- c("snp_name",
  #                        "Beta_lm_o","SE_lm_o","P_wald_lm_o","CI_lower_lm_o","CI_upper_lm_o"
  #                        )
  
  file_path <- paste0("G:/phs000620/MCTFR_clean/Output_end_20241121/", result_name[k], ".csv")
  fwrite(end_result,file_path)
  show_data <- end_result
  
  
  
  snp_tibble <- tibble(SNP_Name = show_data$snp_name)
  
  
  snp_tibble <- snp_tibble %>%
    separate(SNP_Name, into = c("snp_name", "allele"), sep = "_", extra = "merge")
  
  
  show_data[,1]    <- snp_tibble[,1]
  show_data$allele <- snp_tibble$allele
  
  
  show_data_2         <- inner_join(show_data, snp_info, by = c("snp_name" = "snp_name"))
  show_data_3         <- inner_join(show_data_2, snp_info_2[,c(2,5)],by = c("snp_name" = "SNP"))
  
  file_path_2 <- paste0(result_name[k], "_full.csv")
  fwrite(show_data_3,file_path_2)
}

time2 <- Sys.time()
print(time2-time1)

stopImplicitCluster()
stopCluster(cl)
