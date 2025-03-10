library(FTSEM)
library(data.table)


data_1 <- read.table("./Example_2_Exposure.txt", header = T)
data_1 <- as.data.table(data_1)

data_2 <- read.table("./Example_2_Outcome.txt", header = T)
data_2 <- as.data.table(data_2)

# Exposure part

for(i in 1:(ncol(data_1) - 6)){
  snp_name  <- names(data_1[,-(1:6)])
  data_sub  <- data_1[,c(1:6,(6+i)),with = F]
  trio_data <- process_family_data(data_sub, seed = 0)
  
  
  if (i == 1){
    Exposure_result <- FT_SEM(trio_data,snp_name[i])
  } else {
    Exposure_result <- as.data.frame(rbind(Exposure_result,FT_SEM(trio_data,snp_name[i])))
  }
}

# Outcome part

for(i in 1:(ncol(data_2) - 6)){
  snp_name  <- names(data_2[,-(1:6)])
  data_sub  <- data_2[,c(1:6,(6+i)),with = F]
  trio_data <- process_family_data(data_sub, seed = 0)
  
  
  if (i == 1){
    Outcome_result <- FT_SEM(trio_data,snp_name[i])
  } else {
    Outcome_result <- as.data.frame(rbind(Outcome_result,FT_SEM(trio_data,snp_name[i])))
  }
}

# IVW

merged_data <- merge(Exposure_result, Outcome_result, by = "SNP")

bon_p_value <- 0.05 / (nrow(merged_data) - 6)

filtered_data <- subset(merged_data, P_wald_o_e.x < bon_p_value)

beta_exp <- as.numeric(filtered_data$Beta_o_e.x)  
se_exp <- as.numeric(filtered_data$SE_o_e.x)      
beta_out <- as.numeric(filtered_data$Beta_o_e.y)  
se_out <- as.numeric(filtered_data$SE_o_e.y)         


weights <- 1 / se_out^2

beta_ivw <- sum(weights * beta_out / beta_exp) / sum(weights)

se_ivw <- sqrt(1 / sum(weights))

ci_lower <- beta_ivw - 1.96 * se_ivw
ci_upper <- beta_ivw + 1.96 * se_ivw

z_value <- beta_ivw / se_ivw
p_value <- 2 * (1 - pnorm(abs(z_value)))

ivw_result <- data.frame(
  Beta_IVW = beta_ivw,
  SE_IVW = se_ivw,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper,
  P_Value = p_value
)

print(ivw_result)
