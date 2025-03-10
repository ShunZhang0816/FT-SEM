library(MASS)
library(dplyr)
library(OpenMx)
library(foreach)
library(doParallel)
Simulate_IV <- function(Nrep = 100, N=1000, p1=0.25, p2=0.35, beta_o_e=0, beta_f_e=0, beta_m_e=0, beta_u_e=0.5,sd_e = 1,beta_l_e = 0,
                        beta_o_O = 0, beta_f_O = 0, beta_m_O = 0, beta_u_O = 0.4, beta_e_O = 0, beta_l_O = 0, sd_o = 1, rou = 0.6){
  
  
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
                       
                       sigma_e_1 <- matrix(sd_e^2*c(1,rou,rou,rou,1,rou,rou,rou,1),3,3)
                       error_e_1 <- mvrnorm(N/2,c(0,0,0),sigma_e_1)
                       sigma_e_2 <- matrix(sd_e^2*c(1,rou,rou,rou,1,rou,rou,rou,1),3,3)
                       error_e_2 <- mvrnorm(N/2,c(0,0,0),sigma_e_2)
                       

                       Exposure_f_1 <- beta_o_e*f_snp_1 +
                         beta_f_e*gf_snp_1_1 +
                         beta_m_e*gm_snp_1_1 +
                         beta_u_e*U_f_1 +
                         beta_l_e*1 +
                         error_e_1[,1]
                       
                       Exposure_m_1 <- beta_o_e*m_snp_1 +
                         beta_f_e*gf_snp_2_1 +
                         beta_m_e*gm_snp_2_1 +
                         beta_u_e*U_m_1 +
                         beta_l_e*1 +
                         error_e_1[,2]
                       
                       Exposure_o_1 <- beta_o_e*o_snp_1+
                         beta_f_e*f_snp_1 +
                         beta_m_e*m_snp_1 +
                         beta_u_e*U_o_1 +
                         beta_l_e*1 +
                         error_e_1[,3]
                       
                       Exposure_f_2 <- beta_o_e*f_snp_2 +
                         beta_f_e*gf_snp_1_2 +
                         beta_m_e*gm_snp_1_2 +
                         beta_u_e*U_f_2 +
                         beta_l_e*2 +
                         error_e_2[,1]
                       
                       Exposure_m_2 <- beta_o_e*m_snp_2 +
                         beta_f_e*gf_snp_2_2 +
                         beta_m_e*gm_snp_2_2 +
                         beta_u_e*U_m_2 +
                         beta_l_e*2 +
                         error_e_2[,2]
                       
                       Exposure_o_2 <- beta_o_e*o_snp_2+
                         beta_f_e*f_snp_2 +
                         beta_m_e*m_snp_2 +
                         beta_u_e*U_o_2 +
                         beta_l_e*2 +
                         error_e_2[,3]
                       
                       sigma_o_1 <- matrix(sd_o^2*c(1,rou,rou,rou,1,rou,rou,rou,1),3,3)
                       error_o_1 <- mvrnorm(N/2,c(0,0,0),sigma_o_1)
                       sigma_o_2 <- matrix(sd_o^2*c(1,rou,rou,rou,1,rou,rou,rou,1),3,3)
                       error_o_2 <- mvrnorm(N/2,c(0,0,0),sigma_o_2)
                       
                       Outcome_f_1 <- beta_o_O*f_snp_1 +
                         beta_f_O * gf_snp_1_1 +
                         beta_m_O * gm_snp_1_1 +
                         beta_e_O * Exposure_f_1 +
                         beta_u_O * U_f_1 +
                         beta_l_O * 1 +
                         error_o_1[,1]
                       
                       Outcome_m_1 <- beta_o_O*m_snp_1 +
                         beta_f_O * gf_snp_2_1 +
                         beta_m_O * gm_snp_2_1 +
                         beta_e_O * Exposure_m_1 +
                         beta_u_O * U_m_1 +
                         beta_l_O * 1 +
                         error_o_1[,2]
                       
                       Outcome_o_1 <- beta_o_O*o_snp_1 +
                         beta_f_O * f_snp_1 +
                         beta_m_O * m_snp_1 +
                         beta_e_O * Exposure_o_1 +
                         beta_u_O * U_o_1 +
                         beta_l_O * 1 +
                         error_o_1[,3]
                       
                       Outcome_f_2 <- beta_o_O*f_snp_2 +
                         beta_f_O * gf_snp_1_2 +
                         beta_m_O * gm_snp_1_2 +
                         beta_e_O * Exposure_f_2 +
                         beta_u_O * U_f_2 +
                         beta_l_O * 2 +
                         error_o_2[,1]
                       
                       Outcome_m_2 <- beta_o_O*m_snp_2 +
                         beta_f_O * gf_snp_2_2 +
                         beta_m_O * gm_snp_2_2 +
                         beta_e_O * Exposure_m_2 +
                         beta_u_O * U_m_2 +
                         beta_l_O * 2 +
                         error_o_2[,2]
                       
                       Outcome_o_2 <- beta_o_O*o_snp_2 +
                         beta_f_O * f_snp_2 +
                         beta_m_O * m_snp_2 +
                         beta_e_O * Exposure_o_2 +
                         beta_u_O * U_o_2 +
                         beta_l_O * 2 +
                         error_o_2[,3]
                       
                       data_sub <- as.data.frame(cbind(c(f_snp_1,f_snp_2),
                                                       c(m_snp_1,m_snp_2),
                                                       c(o_snp_1,o_snp_2),
                                                       c(Exposure_f_1,Exposure_f_2),
                                                       c(Exposure_m_1,Exposure_m_2),
                                                       c(Exposure_o_1,Exposure_o_2),
                                                       c(Outcome_f_1,Outcome_f_2),
                                                       c(Outcome_m_1,Outcome_m_2),
                                                       c(Outcome_o_1,Outcome_o_2),
                                                       rep(1:2,each=N/2)))
                       names(data_sub) <- c("f_snp","m_snp","o_snp",
                                            "Exposure_f","Exposure_m","Exposure_o",
                                            "Outcome_f","Outcome_m","Outcome_o",
                                            "location")
                       
                       random_index <- sample(nrow(data_sub))
                       data_sub_1 <- data_sub[random_index[1:(N/2)],1:6]
                       data_sub_2 <- data_sub[random_index[(N/2+1):(N)],c(1:3,7:9)]
                       
                       
                       
                       # Set up model parameters
                       manifests <- names(data_sub_1)
                       latents <- c("e_f","e_m","e_o","var_o_snp","var_f_snp","var_m_snp","gf_snp_1","gm_snp_1","gf_snp_2","gm_snp_2")
                       IVmodel_base <- mxModel(model = "IV model", type = "RAM", mxData(observed = data_sub_1, type = "raw"),
                                               
                                               mxPath(from = "gf_snp_1", arrows = 2, free = TRUE, values = 0.2, labels = "theta_snp"),
                                               mxPath(from = "gm_snp_1", arrows = 2, free = TRUE, values = 0.2, labels = "theta_snp"),
                                               mxPath(from = "gf_snp_2", arrows = 2, free = TRUE, values = 0.2, labels = "theta_snp"),
                                               mxPath(from = "gm_snp_2", arrows = 2, free = TRUE, values = 0.2, labels = "theta_snp"),
                                               mxPath(from = "var_f_snp", arrows = 2, free = TRUE, values = 0.2, labels = "theta_snp"),
                                               mxPath(from = "var_m_snp", arrows = 2, free = TRUE, values = 0.2, labels = "theta_snp"),
                                               mxPath(from = "var_o_snp", arrows = 2, free = TRUE, values = 0.2, labels = "theta_snp"),
                                               
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
                       
                       IVFit_base_1  <- mxRun(IVmodel_base)
                       lm_parent_1 <- lm(Exposure_o ~ o_snp + f_snp + m_snp, data_sub_1)
                       lm_1      <- lm(Exposure_o ~ o_snp, data_sub_1)
                       
                       manifests <- names(data_sub_2)
                       latents <- c("e_f","e_m","e_o","var_o_snp","var_f_snp","var_m_snp","gf_snp_1","gm_snp_1","gf_snp_2","gm_snp_2")
                       IVmodel_base <- mxModel(model = "IVmodel", type = "RAM", mxData(observed = data_sub_2, type = "raw"),
                                               
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
                                               mxPath(from = "e_f", to = "Outcome_f", free = FALSE, values = 1, labels = "effect_e_f_O"),
                                               mxPath(from = "e_m", to = "Outcome_m", free = FALSE, values = 1, labels = "effect_e_m_O"),
                                               mxPath(from = "e_o", to = "Outcome_o", free = FALSE, values = 1, labels = "effect_e_o_O"),
                                               
                                               mxPath(from="e_f", to="e_o", arrows=2, free=TRUE, values=0.2, labels=c("ρ")),
                                               mxPath(from="e_m", to="e_o", arrows=2, free=TRUE, values=0.2, labels=c("ρ")),
                                               mxPath(from="e_f", to="e_m", arrows=2, free=TRUE, values=0.2, labels=c("ρ")),
                                               
                                               mxPath(from = "f_snp", to = "Outcome_o", free = TRUE, arrows = 1, labels = "effect_f_O"),
                                               mxPath(from = "m_snp", to = "Outcome_o", free = TRUE, arrows = 1, labels = "effect_m_O"),
                                               mxPath(from = "o_snp", to = "Outcome_o", free = TRUE, arrows = 1, labels = "effect_o_O"),
                                               mxPath(from = "gf_snp_1", to = "Outcome_f", free = TRUE, arrows = 1, labels = "effect_f_O"),
                                               mxPath(from = "gm_snp_1", to = "Outcome_f", free = TRUE, arrows = 1, labels = "effect_m_O"),
                                               mxPath(from = "f_snp", to = "Outcome_f", free = TRUE, arrows = 1, labels = "effect_o_O"),
                                               mxPath(from = "gf_snp_2", to = "Outcome_m", free = TRUE, arrows = 1, labels = "effect_f_O"),
                                               mxPath(from = "gm_snp_2", to = "Outcome_m", free = TRUE, arrows = 1, labels = "effect_m_O"),
                                               mxPath(from = "m_snp", to = "Outcome_m", free = TRUE, arrows = 1, labels = "effect_o_O"),
                                               mxPath( from="one", to=manifests, arrows=1, free= c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
                                                       values=rep(0,length(manifests)), labels=c("mean_snp_f","mean_snp_m","mean_snp_o","mean_ef","mean_em","mean_eo")),  
                                               manifestVars = manifests,latentVars = latents
                       )
                       IVFit_base_2  <- mxRun(IVmodel_base)
                       lm_parent_2   <- lm(Outcome_o ~ o_snp + f_snp + m_snp, data_sub_2)
                       lm_2          <- lm(Outcome_o ~ o_snp, data_sub_2)
                       
                       
                       effect_base          <- (summary(IVFit_base_2)$parameters[1,5]) / (summary(IVFit_base_1)$parameters[1,5])
                       effect_base_std      <- abs((summary(IVFit_base_2)$parameters[1,6]) / (summary(IVFit_base_1)$parameters[1,5]))
                       effect_lm_parent     <- (summary(lm_parent_2)$coefficients[2,1]) / (summary(lm_parent_1)$coefficients[2,1])
                       effect_lm_parent_std <- abs((summary(lm_parent_2)$coefficients[2,2]) / (summary(lm_parent_1)$coefficients[2,1]))
                       effect_lm            <- (summary(lm_2)$coefficients[2,1]) / (summary(lm_1)$coefficients[2,1])
                       effect_lm_std        <- abs((summary(lm_2)$coefficients[2,2]) / (summary(lm_1)$coefficients[2,1]))
                       
                       
                       Zscore_base      <- effect_base/effect_base_std
                       Zscore_lm_parent <- effect_lm_parent/effect_lm_parent_std
                       Zscore_lm        <- effect_lm/effect_lm_std
                       Pval_base        <- 2*(1-pnorm(abs(Zscore_base),0,1))
                       Pval_lm_parent   <- 2*(1-pnorm(abs(Zscore_lm_parent),0,1))
                       Pval_lm          <- 2*(1-pnorm(abs(Zscore_lm),0,1))
                       
                       base      <- rbind(
                         cbind("base",
                               effect_base,
                               effect_base_std,
                               Pval_base,
                               effect_base - 1.96*effect_base_std,
                               effect_base + 1.96*effect_base_std
                         ))
                       lm_parent <- rbind(
                         cbind("lm_parent",
                               effect_lm_parent,
                               effect_lm_parent_std,
                               Pval_lm_parent,
                               effect_lm_parent - 1.96*effect_lm_parent_std,
                               effect_lm_parent + 1.96*effect_lm_parent_std
                         ))
                       lm_result <- rbind(
                         cbind("lm",
                               effect_lm,
                               effect_lm_std,
                               Pval_lm,
                               effect_lm - 1.96*effect_lm_std,
                               effect_lm + 1.96*effect_lm_std
                         ))
                       names(base)      <- c("Method","Estimates","SE","P","CI_lower","CI_upper")
                       names(lm_parent) <- c("Method","Estimates","SE","P","CI_lower","CI_upper")
                       names(lm_result) <- c("Method","Estimates","SE","P","CI_lower","CI_upper")
                       bind_result <- rbind(base,lm_parent,lm_result)
                       bind_result
                     }
  result             <- as.data.frame(result)
  names(result)      <- c("Method","Esitimates","SE","P","CI_lower","CI_upper")
  gc()
  stopCluster(cl)
  time2 <- Sys.time()
  print(time2-time1)
  return(result)
}

Nrep            <- 100
rou             <- c(0.3,0.6)
MAF             <- 0.3  
Nsample         <- c(2000,4000,6000) 
p1               <- 0.25
p2               <- 0.35
result           <- data.frame()
base_result      <- data.frame()
lm_parent_result <- data.frame()
lm_result        <- data.frame()

PVE_zx    <- 0.1                       
PVE_zy    <- 0                        
PVE_ux    <- 0.2                       
PVE_uy    <- 0.3                       
PVE_lx    <- 0.05                      
PVE_ly    <- 0.1                       
sd_e      <- sqrt(1 - PVE_zx - PVE_ux) 
sd_o      <- sqrt(1 - PVE_zy - PVE_uy) 
beta_e_O  <- sqrt(0 / 0.1)             

set.seed(0628)
for (i in 1:length(rou)){
  for (j in 1:length(Nsample)){
    beta_o_e    <- sqrt(0.1 / (2 * MAF * (1 - MAF)))
    beta_f_e    <- sqrt(0 / (2 * MAF * (1 - MAF)))
    beta_m_e    <- sqrt(0 / (2 * MAF * (1 - MAF)))
    beta_u_e    <- sqrt(PVE_ux)
    beta_l_e    <- sqrt(PVE_lx)
    beta_u_O    <- sqrt(PVE_uy)
    beta_f_O    <- sqrt(0 / (2 * MAF * (1 - MAF)))
    beta_m_O    <- sqrt(0 / (2 * MAF * (1 - MAF)))
    beta_l_O    <- sqrt(PVE_ly)
    
    raw_result       <- Simulate_IV(Nrep = Nrep,p1 = p1, p2 = p2, N = Nsample[j],
                                    beta_e_O = beta_e_O,
                                    beta_o_e = beta_o_e, beta_f_e = beta_f_e, beta_m_e = beta_m_e, beta_u_e = beta_u_e, sd_e = sd_e,
                                    beta_l_e = beta_l_e,
                                    beta_o_O = 0, beta_f_O = beta_f_O, beta_m_O = beta_m_O, beta_u_O = beta_u_O, sd_o = sd_o, rou = rou[i],
                                    beta_l_O = beta_l_O)
    raw_result       <- as.data.frame(raw_result)
    raw_result[, -1] <- sapply(raw_result[, -1], as.numeric)
    
    ## Bias
    base_bias      <- mean(raw_result[raw_result[,1]=="base",2] - beta_e_O)
    lm_parent_bias <- mean(raw_result[raw_result[,1]=="lm_parent",2] - beta_e_O)
    lm_bias       <- mean(raw_result[raw_result[,1]=="lm",2] - beta_e_O)
    
    ## RMSE
    base_RMSE      <- sqrt(mean((raw_result[raw_result[,1]=="base",2] - beta_e_O)^2))
    lm_parent_RMSE <- sqrt(mean((raw_result[raw_result[,1]=="lm_parent",2] - beta_e_O)^2))
    lm_RMSE        <- sqrt(mean((raw_result[raw_result[,1]=="lm",2] - beta_e_O)^2))
    
    base_coverage    <- mean(
      (raw_result[raw_result[,1]=="base",5] <= beta_e_O) & (beta_e_O < raw_result[raw_result[,1]=="base",6])
    )
    lm_parent_coverage <- mean(
      (raw_result[raw_result[,1]=="lm_parent",5] <= beta_e_O) & (beta_e_O < raw_result[raw_result[,1]=="lm_parent",6])
    )
    lm_coverage <- mean(
      (raw_result[raw_result[,1]=="lm",5] <= beta_e_O) & (beta_e_O < raw_result[raw_result[,1]=="lm",6])
    )
    
    base_CIlength      <- mean(
      raw_result[raw_result[,1]=="base",6] - raw_result[raw_result[,1]=="base",5]
    )
    lm_parent_CIlength <- mean(
      raw_result[raw_result[,1]=="lm_parent",6] - raw_result[raw_result[,1]=="lm_parent",5]
    )
    lm_CIlength <- mean(
      raw_result[raw_result[,1]=="lm",6] - raw_result[raw_result[,1]=="lm",5]
    )
    
    base_a      <- mean(raw_result[raw_result[,1]=="base",4]<0.05)
    lm_parent_a <- mean(raw_result[raw_result[,1]=="lm_parent",4]<0.05)
    lm_a        <- mean(raw_result[raw_result[,1]=="lm",4]<0.05)
    
    base_result <- rbind(base_result,
                         cbind(Nrep,p1,p2,rou[i],Nsample[j],"FT-SEM",
                               base_bias,base_RMSE,base_coverage,base_CIlength,base_a
                         ))
    lm_parent_result <- rbind(lm_parent_result,
                              cbind(Nrep,p1,p2,rou[i],Nsample[j],"lm_parent",
                                    lm_parent_bias,lm_parent_RMSE,lm_parent_coverage,
                                    lm_parent_CIlength,lm_parent_a
                              ))
    lm_result <- rbind(lm_result,
                       cbind(Nrep,p1,p2,rou[i],Nsample[j],"lm",
                             lm_bias,lm_RMSE,lm_coverage,lm_CIlength,lm_a
                       ))
  }
}
names(base_result)   <- c("Simulation","MAF_1","MAF_2","rou","Sample_size","Method",
                          "Bias","RMSE","Coverage","CILength","Type_one_error_rates"
)
names(lm_parent_result)   <- c("Simulation","MAF_1","MAF_2","rou","Sample_size","Method",
                               "Bias","RMSE","Coverage","CILength","Type_one_error_rates"
)
names(lm_result)   <- c("Simulation","MAF_1","MAF_2","rou","Sample_size","Method",
                        "Bias","RMSE","Coverage","CILength","Type_one_error_rates"
)
result               <- rbind(base_result,lm_parent_result,lm_result)
result