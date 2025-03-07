library(MASS)
library(dplyr)
library(OpenMx)
library(foreach)
library(doParallel)
Simulate_IV <- function(Nrep = 100, N=1000, p1=0.35,p2=0.25, beta_o_e=0, beta_f_e=0, beta_m_e=0, beta_u_e=0.5,sd_e = 1,beta_l_e = 0, rou = 0.5){

  
  q1 <- 1 - p1                      
  q2 <- 1 - p2
  base             <- data.frame() 
  lm_result        <- data.frame() 
  lm_parent_result <- data.frame() 
  
  cl <- makeCluster(24)
  registerDoParallel(cl)
  time1 <- Sys.time()
  result <- foreach (j = 1:Nrep,
                     .combine = rbind,
                     .packages = c("data.table","dplyr","OpenMx","MASS"),
                     .errorhandling = "pass") %dopar% {
                       ### SIMULATE DATA ###
                       
                       gf_snp_1_1 <- sample(x = c(0,1,2), size = N/2, replace = TRUE, prob = c(q1^2, 2*p1*q1, p1^2))
                       gm_snp_1_1 <- sample(x = c(0,1,2), size = N/2, replace = TRUE, prob = c(q1^2, 2*p1*q1, p1^2))
                       gf_snp_1_2 <- sample(x = c(0,1,2), size = N/2, replace = TRUE, prob = c(q2^2, 2*p2*q2, p2^2))
                       gm_snp_1_2 <- sample(x = c(0,1,2), size = N/2, replace = TRUE, prob = c(q2^2, 2*p2*q2, p2^2))
                       gf_snp_2_1 <- sample(x = c(0,1,2), size = N/2, replace = TRUE, prob = c(q1^2, 2*p1*q1, p1^2))
                       gm_snp_2_1 <- sample(x = c(0,1,2), size = N/2, replace = TRUE, prob = c(q1^2, 2*p1*q1, p1^2))
                       gf_snp_2_2 <- sample(x = c(0,1,2), size = N/2, replace = TRUE, prob = c(q2^2, 2*p2*q2, p2^2))
                       gm_snp_2_2 <- sample(x = c(0,1,2), size = N/2, replace = TRUE, prob = c(q2^2, 2*p2*q2, p2^2))
                       
                       snp_inherit <- function(f_snp,m_snp,N){
                         o_snp <- vector(length=N)
                         ro    <- runif(N)
                         for (i in 1:N) { 
                           if((m_snp[i]==0) && (f_snp[i]==0)) {o_snp[i] = 0}
                           if((m_snp[i]==0) && (f_snp[i]==1))  {if(ro[i] <= 0.5) {o_snp[i] = 0} else {o_snp[i] = 1}}
                           if((m_snp[i]==0) && (f_snp[i]==2))  {o_snp[i] = 1}
                           if((m_snp[i]==1) &&  (f_snp[i]==0)) {if(ro[i] <= 0.5) {o_snp[i] = 0} else {o_snp[i] = 1}}
                           if((m_snp[i]==1) &&  (f_snp[i]==1))  {
                             if(ro[i] <= 0.25) {o_snp[i] = 0}
                             if(ro[i] > 0.25 && ro[i] <= 0.75) {o_snp[i] = 1}
                             if(ro[i] > 0.75) {o_snp[i] = 2}
                           }
                           if((m_snp[i]==1) &&  (f_snp[i]==2))  {if(ro[i] <= 0.5) {o_snp[i] = 2} else {o_snp[i] = 1}}
                           if((m_snp[i]==2) &&  (f_snp[i]==0)) {o_snp[i] = 1}
                           if((m_snp[i]==2) &&  (f_snp[i]==1))  {if(ro[i] <= 0.5) {o_snp[i] = 2} else {o_snp[i] = 1}}
                           if((m_snp[i]==2) &&  (f_snp[i]==2))  {o_snp[i] = 2}
                         }
                         return(o_snp)
                       }
                       f_snp_1 <- snp_inherit(gf_snp_1_1,gm_snp_1_1,N/2)
                       f_snp_2 <- snp_inherit(gf_snp_1_2,gm_snp_1_2,N/2)
                       m_snp_1 <- snp_inherit(gf_snp_2_1,gm_snp_2_1,N/2)
                       m_snp_2 <- snp_inherit(gf_snp_2_2,gm_snp_2_2,N/2)
                       o_snp_1 <- snp_inherit(f_snp_1,m_snp_1,N/2)
                       o_snp_2 <- snp_inherit(f_snp_2,m_snp_2,N/2)
                       
                       
                       U_f_1 <- rnorm(N/2, mean=0, sd=1)
                       U_f_2 <- rnorm(N/2, mean=0, sd=1)
                       U_m_1 <- rnorm(N/2, mean=0, sd=1)
                       U_m_2 <- rnorm(N/2, mean=0, sd=1)
                       U_o_1 <- rnorm(N/2, mean=0, sd=1)
                       U_o_2 <- rnorm(N/2, mean=0, sd=1)
                       
                       sigma_1 <- matrix(sd_e^2*c(1,rou,rou,rou,1,rou,rou,rou,1),3,3)
                       error_1 <- mvrnorm(N/2,c(0,0,0),sigma_1)
                       
                       sigma_2 <- matrix(sd_e^2*c(1,rou,rou,rou,1,rou,rou,rou,1),3,3)
                       error_2 <- mvrnorm(N/2,c(0,0,0),sigma_2)
                       
                       Exposure_f_1 <- beta_o_e*f_snp_1 +
                         beta_f_e*gf_snp_1_1 +
                         beta_m_e*gm_snp_1_1 +
                         beta_u_e*U_f_1 +
                         beta_l_e*1 +
                         error_1[,1]
                       
                       Exposure_m_1 <- beta_o_e*m_snp_1 +
                         beta_f_e*gf_snp_2_1 +
                         beta_m_e*gm_snp_2_1 +
                         beta_u_e*U_m_1 +
                         beta_l_e*1 +
                         error_1[,2]
                       
                       Exposure_o_1 <- beta_o_e*o_snp_1+
                         beta_f_e*f_snp_1 +
                         beta_m_e*m_snp_1 +
                         beta_u_e*U_o_1 +
                         beta_l_e*1 +
                         error_1[,3]
                       
                       Exposure_f_2 <- beta_o_e*f_snp_2 +
                         beta_f_e*gf_snp_1_2 +
                         beta_m_e*gm_snp_1_2 +
                         beta_u_e*U_f_2 +
                         beta_l_e*2 +
                         error_2[,1]
                       
                       Exposure_m_2 <- beta_o_e*m_snp_2 +
                         beta_f_e*gf_snp_2_2 +
                         beta_m_e*gm_snp_2_2 +
                         beta_u_e*U_m_2 +
                         beta_l_e*2 +
                         error_2[,2]
                       
                       Exposure_o_2 <- beta_o_e*o_snp_2+
                         beta_f_e*f_snp_2 +
                         beta_m_e*m_snp_2 +
                         beta_u_e*U_o_2 +
                         beta_l_e*2 +
                         error_2[,3]
                       
                       data_sub <- as.data.frame(cbind(c(f_snp_1,f_snp_2),
                                                       c(m_snp_1,m_snp_2),
                                                       c(o_snp_1,o_snp_2),
                                                       c(Exposure_f_1,Exposure_f_2),
                                                       c(Exposure_m_1,Exposure_m_2),
                                                       c(Exposure_o_1,Exposure_o_2),
                                                       rep(1:2,each=N/2)))
                       names(data_sub) <- c("f_snp","m_snp","o_snp",
                                            "Exposure_f","Exposure_m","Exposure_o",
                                            "location")
                       
                       
                       
                       # Set up model parameters
                       manifests <- c("f_snp","m_snp","o_snp",
                                      "Exposure_f","Exposure_m","Exposure_o")
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
                                               
                                               mxPath(from = "e_f", arrows = 2, free = TRUE, values = 1, labels = "var_e_f"),
                                               mxPath(from = "e_m", arrows = 2, free = TRUE, values = 1, labels = "var_e_m"),
                                               mxPath(from = "e_o", arrows = 2, free = TRUE, values = 1, labels = "var_e_o"),
                                               mxPath(from = "e_f", to = "Exposure_f", free = FALSE, values = 1,
                                                      labels = "effect_e_f_E"),
                                               mxPath(from = "e_m", to = "Exposure_m", free = FALSE, values = 1,
                                                      labels = "effect_e_m_E"),
                                               mxPath(from = "e_o", to = "Exposure_o", free = FALSE, values = 1,
                                                      labels = "effect_e_o_E"),
                                               
                                               mxPath(from="e_f", to="e_o", arrows=2, free=TRUE, values=0.2, labels=c("ρ")),
                                               mxPath(from="e_m", to="e_o", arrows=2, free=TRUE, values=0.2, labels=c("ρ")),
                                               mxPath(from="e_f", to="e_m", arrows=2, free=TRUE, values=0.2, labels=c("ρ")),
                                               
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
                       lm_parent   <- lm(Exposure_o ~ o_snp + f_snp + m_snp, data_sub)
                       lm_o        <- lm(Exposure_o ~ o_snp, data_sub)
                       lm_f        <- lm(Exposure_o ~ f_snp, data_sub)
                       lm_m        <- lm(Exposure_o ~ m_snp, data_sub)
                       
                       
                       base_summary      <- summary(IVFit_base)
                       lm_parent_summary <- summary(lm_parent)
                       lm_f_summary      <- summary(lm_f)
                       lm_m_summary      <- summary(lm_m)
                       lm_o_summary      <- summary(lm_o)
                       
                       # Save results
                       base <- rbind(c("base",
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
                       
                       lm_parent_result <- rbind(
                         c("lm_parent",
                           lm_parent_summary$coefficients[2,1],
                           lm_parent_summary$coefficients[2,2],lm_parent_summary$coefficients[2,4],
                           (lm_parent_summary$coefficients[2,1]-1.96*lm_parent_summary$coefficients[2,2]),
                           (lm_parent_summary$coefficients[2,1]+1.96*lm_parent_summary$coefficients[2,2]),
                           lm_parent_summary$coefficients[3,1],
                           lm_parent_summary$coefficients[3,2],lm_parent_summary$coefficients[3,4],
                           (lm_parent_summary$coefficients[3,1]-1.96*lm_parent_summary$coefficients[3,2]),
                           (lm_parent_summary$coefficients[3,1]+1.96*lm_parent_summary$coefficients[3,2]),
                           lm_parent_summary$coefficients[4,1],
                           lm_parent_summary$coefficients[4,2],lm_parent_summary$coefficients[4,4],
                           (lm_parent_summary$coefficients[4,1]-1.96*lm_parent_summary$coefficients[4,2]),
                           (lm_parent_summary$coefficients[4,1]+1.96*lm_parent_summary$coefficients[4,2])
                         )
                       )
                       lm_result <- rbind(
                         c("lm",
                           lm_o_summary$coefficients[2,1],
                           lm_o_summary$coefficients[2,2],lm_o_summary$coefficients[2,4],
                           (lm_o_summary$coefficients[2,1]-1.96*lm_o_summary$coefficients[2,2]),
                           (lm_o_summary$coefficients[2,1]+1.96*lm_o_summary$coefficients[2,2]),
                           lm_f_summary$coefficients[2,1],
                           lm_f_summary$coefficients[2,2],lm_f_summary$coefficients[2,4],
                           (lm_f_summary$coefficients[2,1]-1.96*lm_f_summary$coefficients[2,2]),
                           (lm_f_summary$coefficients[2,1]+1.96*lm_f_summary$coefficients[2,2]),
                           lm_m_summary$coefficients[2,1],
                           lm_m_summary$coefficients[2,2],lm_m_summary$coefficients[2,4],
                           (lm_m_summary$coefficients[2,1]-1.96*lm_m_summary$coefficients[2,2]),
                           (lm_m_summary$coefficients[2,1]+1.96*lm_m_summary$coefficients[2,2])
                         )
                       )
                       names(base) <- c("Method",
                                        "Beta_o_e", "SE_o_e", "P_wald_o_e","CI_lower_o_e","CI_upper_o_e",
                                        "Beta_f_e", "SE_f_e", "P_wald_f_e","CI_lower_f_e","CI_upper_f_e",
                                        "Beta_m_e", "SE_m_e", "P_wald_m_e","CI_lower_m_e","CI_upper_m_e"
                       )
                       names(lm_parent_result) <- c("Method",
                                                    "Beta_o_e", "SE_o_e", "P_wald_o_e","CI_lower_o_e","CI_upper_o_e",
                                                    "Beta_f_e", "SE_f_e", "P_wald_f_e","CI_lower_f_e","CI_upper_f_e",
                                                    "Beta_m_e", "SE_m_e", "P_wald_m_e","CI_lower_m_e","CI_upper_m_e"
                       )
                       names(lm_result) <- c("Method",
                                             "Beta_o_e", "SE_o_e", "P_wald_o_e","CI_lower_o_e","CI_upper_o_e",
                                             "Beta_f_e", "SE_f_e", "P_wald_f_e","CI_lower_f_e","CI_upper_f_e",
                                             "Beta_m_e", "SE_m_e", "P_wald_m_e","CI_lower_m_e","CI_upper_m_e"
                       )
                       bind_result <- rbind(base,lm_parent_result,lm_result)
                       bind_result
                     }
  result        <- as.data.frame(result)
  names(result) <- c("Method",
                     "Beta_o_e", "SE_o_e", "P_wald_o_e","CI_lower_o_e","CI_upper_o_e",
                     "Beta_f_e", "SE_f_e", "P_wald_f_e","CI_lower_f_e","CI_upper_f_e",
                     "Beta_m_e", "SE_m_e", "P_wald_m_e","CI_lower_m_e","CI_upper_m_e"
  )
  gc()
  stopCluster(cl)
  time2 <- Sys.time()
  print(time2-time1)
  return(result)
}

Nrep             <- 10000  
MAF              <- 0.3 
rou              <- c(0.3,0.6)
p1               <- 0.25
p2               <- 0.35
Nsample          <- c(1000,2000,3000)  
base_result      <- data.frame() 
lm_parent_result <- data.frame()
lm_result        <- data.frame()


PVE_snp_sum <- 0    
PVE_u       <- 0.3  
PVE_l       <- 0.05  
sd_e        <- sqrt(1 - PVE_snp_sum - PVE_u) 
beta_f_e    <- 0
beta_m_e    <- 0
set.seed(0628)
for (i in 1:length(rou)){
  for (j in 1:length(Nsample)){
    beta_o_e    <- sqrt(0 / (2 * MAF * (1 - MAF)))
    beta_u_e    <- sqrt(PVE_u)
    beta_l_e    <- sqrt(PVE_l / 0.25)
    raw_result       <- Simulate_IV(Nrep = Nrep, N = Nsample[j], p1 = p1, p2 = p2,
                                    beta_o_e = beta_o_e, beta_u_e = beta_u_e, beta_l_e = beta_l_e, beta_f_e = beta_f_e, beta_m_e = beta_m_e,
                                    sd_e = sd_e, rou = rou[i]
    )
    raw_result       <- as.data.frame(raw_result)
    raw_result[, -1] <- sapply(raw_result[, -1], as.numeric)
    
    base_bias_o_e      <- mean(raw_result[raw_result[,1]=="base",2] - beta_o_e)
    lm_parent_bias_o_e <- mean(raw_result[raw_result[,1]=="lm_parent",2] - beta_o_e)
    lm_bias_o_e        <- mean(raw_result[raw_result[,1]=="lm",2] - beta_o_e)
    
    base_bias_f_e      <- mean(raw_result[raw_result[,1]=="base",7] - beta_f_e)
    lm_parent_bias_f_e <- mean(raw_result[raw_result[,1]=="lm_parent",7] - beta_f_e)
    lm_bias_f_e        <- mean(raw_result[raw_result[,1]=="lm",7] - beta_f_e)
    
    base_bias_m_e      <- mean(raw_result[raw_result[,1]=="base",12] - beta_m_e)
    lm_parent_bias_m_e <- mean(raw_result[raw_result[,1]=="lm_parent",12] - beta_m_e)
    lm_bias_m_e        <- mean(raw_result[raw_result[,1]=="lm",12] - beta_m_e)
    
    ## RMSE
    base_RMSE_o_e      <- sqrt(mean((raw_result[raw_result[,1]=="base",2] - beta_o_e)^2))
    lm_parent_RMSE_o_e <- sqrt(mean((raw_result[raw_result[,1]=="lm_parent",2] - beta_o_e)^2))
    lm_RMSE_o_e        <- sqrt(mean((raw_result[raw_result[,1]=="lm",2] - beta_o_e)^2))
    
    base_RMSE_f_e      <- sqrt(mean((raw_result[raw_result[,1]=="base",7] - beta_f_e)^2))
    lm_parent_RMSE_f_e <- sqrt(mean((raw_result[raw_result[,1]=="lm_parent",7] - beta_f_e)^2))
    lm_RMSE_f_e        <- sqrt(mean((raw_result[raw_result[,1]=="lm",7] - beta_f_e)^2))
    
    base_RMSE_m_e      <- sqrt(mean((raw_result[raw_result[,1]=="base",12] - beta_m_e)^2))
    lm_parent_RMSE_m_e <- sqrt(mean((raw_result[raw_result[,1]=="lm_parent",12] - beta_m_e)^2))
    lm_RMSE_m_e        <- sqrt(mean((raw_result[raw_result[,1]=="lm",12] - beta_m_e)^2))
    
    base_coverage_o_e      <- mean(
      (raw_result[raw_result[,1]=="base",5] <= beta_o_e) & (beta_o_e < raw_result[raw_result[,1]=="base",6])
    )
    lm_parent_coverage_o_e <- mean(
      (raw_result[raw_result[,1]=="lm_parent",5] <= beta_o_e) & (beta_o_e < raw_result[raw_result[,1]=="lm_parent",6])
    )
    lm_coverage_o_e <- mean(
      (raw_result[raw_result[,1]=="lm",5] <= beta_o_e) & (beta_o_e < raw_result[raw_result[,1]=="lm",6])
    )
    
    base_coverage_f_e      <- mean(
      (raw_result[raw_result[,1]=="base",10] <= beta_f_e) & (beta_f_e < raw_result[raw_result[,1]=="base",11])
    )
    lm_parent_coverage_f_e <- mean(
      (raw_result[raw_result[,1]=="lm_parent",10] <= beta_f_e) & (beta_f_e < raw_result[raw_result[,1]=="lm_parent",11])
    )
    lm_coverage_f_e <- mean(
      (raw_result[raw_result[,1]=="lm",10] <= beta_f_e) & (beta_f_e < raw_result[raw_result[,1]=="lm",11])
    )
    
    base_coverage_m_e      <- mean(
      (raw_result[raw_result[,1]=="base",15] <= beta_m_e) & (beta_m_e < raw_result[raw_result[,1]=="base",16])
    )
    lm_parent_coverage_m_e <- mean(
      (raw_result[raw_result[,1]=="lm_parent",15] <= beta_m_e) & (beta_m_e < raw_result[raw_result[,1]=="lm_parent",16])
    )
    lm_coverage_m_e <- mean(
      (raw_result[raw_result[,1]=="lm",15] <= beta_m_e) & (beta_m_e < raw_result[raw_result[,1]=="lm",16])
    )
    
    base_CIlength_o_e      <- mean(
      raw_result[raw_result[,1]=="base",6] - raw_result[raw_result[,1]=="base",5]
    )
    lm_parent_CIlength_o_e <- mean(
      raw_result[raw_result[,1]=="lm_parent",6] - raw_result[raw_result[,1]=="lm_parent",5]
    )
    lm_CIlength_o_e <- mean(
      raw_result[raw_result[,1]=="lm",6] - raw_result[raw_result[,1]=="lm",5]
    )
    
    base_CIlength_f_e      <- mean(
      raw_result[raw_result[,1]=="base",11] - raw_result[raw_result[,1]=="base",10]
    )
    lm_parent_CIlength_f_e <- mean(
      raw_result[raw_result[,1]=="lm_parent",11] - raw_result[raw_result[,1]=="lm_parent",10]
    )
    lm_CIlength_f_e <- mean(
      raw_result[raw_result[,1]=="lm",11] - raw_result[raw_result[,1]=="lm",10]
    )
    
    base_CIlength_m_e      <- mean(
      raw_result[raw_result[,1]=="base",16] - raw_result[raw_result[,1]=="base",15]
    )
    lm_parent_CIlength_m_e <- mean(
      raw_result[raw_result[,1]=="lm_parent",16] - raw_result[raw_result[,1]=="lm_parent",15]
    )
    lm_CIlength_m_e <- mean(
      raw_result[raw_result[,1]=="lm",16] - raw_result[raw_result[,1]=="lm",15]
    )
    
    base_a_o_e      <- mean(raw_result[raw_result[,1]=="base",4]<0.05)
    lm_parent_a_o_e <- mean(raw_result[raw_result[,1]=="lm_parent",4]<0.05)
    lm_a_o_e        <- mean(raw_result[raw_result[,1]=="lm",4]<0.05)
    
    base_a_f_e      <- mean(raw_result[raw_result[,1]=="base",9]<0.05)
    lm_parent_a_f_e <- mean(raw_result[raw_result[,1]=="lm_parent",9]<0.05)
    lm_a_f_e        <- mean(raw_result[raw_result[,1]=="lm",9]<0.05)
    
    base_a_m_e      <- mean(raw_result[raw_result[,1]=="base",14]<0.05)
    lm_parent_a_m_e <- mean(raw_result[raw_result[,1]=="lm_parent",14]<0.05)
    lm_a_m_e        <- mean(raw_result[raw_result[,1]=="lm",14]<0.05)
    
    base_result <- rbind(base_result,
                         cbind(Nrep,"FT-SEM",p1,p2,rou[i],Nsample[j],
                               base_bias_o_e,base_RMSE_o_e,base_coverage_o_e,base_CIlength_o_e,base_a_o_e,
                               base_bias_f_e,base_RMSE_f_e,base_coverage_f_e,base_CIlength_f_e,base_a_f_e,
                               base_bias_m_e,base_RMSE_m_e,base_coverage_m_e,base_CIlength_m_e,base_a_m_e
                         ))
    lm_parent_result <- rbind(lm_parent_result,
                              cbind(Nrep,"lm_parent",p1,p2,rou[i],Nsample[j],
                                    lm_parent_bias_o_e,lm_parent_RMSE_o_e,lm_parent_coverage_o_e,
                                    lm_parent_CIlength_o_e,lm_parent_a_o_e,
                                    lm_parent_bias_f_e,lm_parent_RMSE_f_e,lm_parent_coverage_f_e,
                                    lm_parent_CIlength_f_e,lm_parent_a_f_e,
                                    lm_parent_bias_m_e,lm_parent_RMSE_m_e,lm_parent_coverage_m_e,
                                    lm_parent_CIlength_m_e,lm_parent_a_m_e
                              ))
    lm_result <- rbind(lm_result,
                       cbind(Nrep,"lm",p1,p2,rou[i],Nsample[j],
                             lm_bias_o_e,lm_RMSE_o_e,lm_coverage_o_e,lm_CIlength_o_e,lm_a_o_e,
                             lm_bias_f_e,lm_RMSE_f_e,lm_coverage_f_e,lm_CIlength_f_e,lm_a_f_e,
                             lm_bias_m_e,lm_RMSE_m_e,lm_coverage_m_e,lm_CIlength_m_e,lm_a_m_e
                       ))
  }
}
names(base_result)   <- c("Simulation","Method","MAF_1","MAF_2","rou","Sample_size",
                          "Offspring_effect_Bias","Offspring_effect_RMSE"," Offspring_effect_Coverage","Offspring_effect_CILength","Offspring_effect_Type_one_error_rates",
                          "Paternal_effect_Bias","Paternal_effect_RMSE","Paternal_effect_Coverage","Paternal_effect_CILength","Paternal_effect_Type_one_error_rates",
                          "Maternal_effect_Bias","Maternal_effect_RMSE","Maternal_effect_Coverage","Maternal_effect_CILength","Maternal_effect_Type_one_error_rates"
)
names(lm_parent_result)   <- c("Simulation","Method","MAF_1","MAF_2","rou","Sample_size",
                               "Offspring Bias","Offspring_effect_RMSE","Offspring_effect_Coverage","Offspring_effect_CILength","Offspring_effect_Type_one_error_rates",
                               "Paternal effectBias","Paternal effectRMSE","Paternal effectCoverage","Paternal effectCILength","Paternal_effect_Type_one_error_rates",
                               "Maternal_effect_Bias","Maternal_effect_RMSE","Maternal_effect_Coverage","Maternal_effect_CILength","Maternal_effect_Type_one_error_rates"
)
names(lm_result)   <- c("Simulation","Method","MAF_1","MAF_2","rou","Sample_size",
                        "Offspring Bias","Offspring_effect_RMSE","Offspring_effect_Coverage","Offspring_effect_CILength","Offspring_effect_Type_one_error_rates",
                        "Paternal_effect_Bias","Paternal_effect_RMSE","Paternal_effect_Coverage","Paternal_effect_CILength","Paternal_effect_Type_one_error_rates",
                        "Maternal_effect_Bias","Maternal_effect_RMSE","Maternal_effect_Coverage","Maternal_effect_CILength","Maternal_effect_Type_one_error_rates"
)
result               <- rbind(base_result,lm_parent_result,lm_result)
