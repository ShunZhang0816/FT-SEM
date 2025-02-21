library(data.table)
library(TwoSampleMR)
library(dplyr)

# Modify the address of summary results about phs000620 before run the 

dataset <- "ieu-b-40"
# 综合分析的文件名
result_name      <- c("NIC_parents","CON_parents","DEP_parents","DRG_parents","BD_parents","FSIQ")
bind_result      <- data.frame()
for (k in 1:length(result_name)){
  instrument <- extract_instruments(
    dataset,                     # 指定访问数据集
    p1 = 5e-08,                       # 过滤 SNP 的 p 值阈值
    clump = FALSE                     # 进行 SNP clumping
  )
  
  Outcome <- read_outcome_data(
    filename = paste0("G:/phs000620/MCTFR_clean/Output_end_20241121/", result_name[k], "_full.csv"),
    snps = instrument$SNP,
    sep = ",",
    snp_col = "snp_name",
    beta_col = "Beta_lm_o",
    se_col = "SE_lm_o",
    pval_col = "P_wald_lm_o",
    chr_col = "chr",
    pos_col = "location2",
    effect_allele_col = "minor_allele",
    other_allele_col = "major_allele",
    eaf_col = "MAF"
  )
  
  instrument         <- instrument[instrument$SNP %in% Outcome$SNP,]
  
  instrument_clump   <- clump_data(instrument,clump_r2=0.001 ,pop = "EUR")
  instrument         <- instrument %>% filter(SNP %in% instrument_clump$SNP)
  
  # 筛选工具变量（去除弱工具变量偏倚）
  Ffilter = 10
  dat <- instrument
  dat$eaf.exposure <- Outcome[Outcome$SNP %in% instrument$SNP,]$eaf.outcome
  N=dat[1,"samplesize.exposure"]
  dat=transform(dat,R2=2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))
  dat=transform(dat,F=(N-2)*R2/(1-R2))
  outTab=dat[dat$F>Ffilter,]
  
  instrument <- instrument[instrument$SNP %in% outTab$SNP,]
  
  # 对数据进行预处理，使其效应等位与效应量保持统一
  dat <- harmonise_data(
    exposure_dat = instrument, 
    outcome_dat = Outcome
  )
  
  res         <- mr(dat)
  ple         <- mr_pleiotropy_test(dat)
  result      <- res[,c(4,3,5:9)]
  result[,1]  <- "BMI"
  result[,2]  <- result_name[k]
  
  result$egger_intercept <- ple$egger_intercept
  result$egger_se        <- ple$se
  result$egger_pval      <- ple$pval
  if (k == 1){
    bind_result <- result
  } else {
    bind_result <- rbind(bind_result,
                         result
    )
  }
}

# 筛选工具变量（去除弱工具变量偏倚）
Ffilter = 10
dat <- instrument
dat$eaf.exposure <- Outcome[Outcome$SNP %in% instrument$SNP,]$eaf.outcome
N=dat[1,"samplesize.exposure"]
dat=transform(dat,R2=2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))
dat=transform(dat,F=(N-2)*R2/(1-R2))
outTab=dat[dat$F>Ffilter,]

dataset <- "ieu-b-4816"
# 综合分析的文件名
result_name      <- c("NIC_parents","CON_parents","DEP_parents","DRG_parents","BD_parents","FSIQ")
bind_result      <- data.frame()
for (k in 1:length(result_name)){
  instrument <- extract_instruments(
    dataset,                     # 指定访问数据集
    p1 = 5e-08,                       # 过滤 SNP 的 p 值阈值
    clump = FALSE                     # 进行 SNP clumping
  )
  
  Outcome <- read_outcome_data(
    filename = paste0("G:/phs000620/MCTFR_clean/Output_end_20241121/", result_name[k], "_full.csv"),
    snps = instrument$SNP,
    sep = ",",
    snp_col = "snp_name",
    beta_col = "Beta_o_base",
    se_col = "SE_o_base",
    pval_col = "P_wald_o_base",
    chr_col = "chr",
    pos_col = "location2",
    effect_allele_col = "minor_allele",
    other_allele_col = "major_allele",
    eaf_col = "MAF"
  )
  
  instrument         <- instrument[instrument$SNP %in% Outcome$SNP,]
  
  instrument_clump   <- clump_data(instrument,clump_r2=0.001 ,pop = "EUR")
  instrument         <- instrument %>% filter(SNP %in% instrument_clump$SNP)
  
  # 筛选工具变量（去除弱工具变量偏倚）
  Ffilter = 10
  dat <- instrument
  dat$eaf.exposure <- Outcome[Outcome$SNP %in% instrument$SNP,]$eaf.outcome
  N=dat[1,"samplesize.exposure"]
  dat=transform(dat,R2=2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))
  dat=transform(dat,F=(N-2)*R2/(1-R2))
  outTab=dat[dat$F>Ffilter,]
  
  instrument <- instrument[instrument$SNP %in% outTab$SNP,]
  
  # 对数据进行预处理，使其效应等位与效应量保持统一
  dat <- harmonise_data(
    exposure_dat = instrument, 
    outcome_dat = Outcome
  )
  
  res         <- mr(dat)
  ple         <- mr_pleiotropy_test(dat)
  result      <- res[,c(4,3,5:9)]
  result[,1]  <- "BMI"
  result[,2]  <- result_name[k]
  
  result$egger_intercept <- ple$egger_intercept
  result$egger_se        <- ple$se
  result$egger_pval      <- ple$pval
  if (k == 1){
    bind_result <- result
  } else {
    bind_result <- rbind(bind_result,
                         result
    )
  }
}

dataset <- "ieu-b-4816"
# 综合分析的文件名
result_name      <- c("NIC_parents","CON_parents","DEP_parents","DRG_parents","BD_parents","FSIQ")
bind_result      <- data.frame()
for (k in 1:length(result_name)){
  instrument <- extract_instruments(
    dataset,                     # 指定访问数据集
    p1 = 5e-08,                       # 过滤 SNP 的 p 值阈值
    clump = FALSE                     # 进行 SNP clumping
  )
  
  Outcome <- read_outcome_data(
    filename = paste0("G:/phs000620/MCTFR_clean/Output_end_20241121/", result_name[k], "_full.csv"),
    snps = instrument$SNP,
    sep = ",",
    snp_col = "snp_name",
    beta_col = "Beta_lm_parent_o",
    se_col = "SE_lm_parent_o",
    pval_col = "P_wald_lm_parent_o",
    chr_col = "chr",
    pos_col = "location2",
    effect_allele_col = "minor_allele",
    other_allele_col = "major_allele",
    eaf_col = "MAF"
  )
  
  instrument         <- instrument[instrument$SNP %in% Outcome$SNP,]
  
  instrument_clump   <- clump_data(instrument,clump_r2=0.001 ,pop = "EUR")
  instrument         <- instrument %>% filter(SNP %in% instrument_clump$SNP)
  
  # 筛选工具变量（去除弱工具变量偏倚）
  Ffilter = 10
  dat <- instrument
  dat$eaf.exposure <- Outcome[Outcome$SNP %in% instrument$SNP,]$eaf.outcome
  N=dat[1,"samplesize.exposure"]
  dat=transform(dat,R2=2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))
  dat=transform(dat,F=(N-2)*R2/(1-R2))
  outTab=dat[dat$F>Ffilter,]
  
  instrument <- instrument[instrument$SNP %in% outTab$SNP,]
  
  # 对数据进行预处理，使其效应等位与效应量保持统一
  dat <- harmonise_data(
    exposure_dat = instrument, 
    outcome_dat = Outcome
  )
  
  res         <- mr(dat)
  ple         <- mr_pleiotropy_test(dat)
  result      <- res[,c(4,3,5:9)]
  result[,1]  <- "BMI"
  result[,2]  <- result_name[k]
  
  result$egger_intercept <- ple$egger_intercept
  result$egger_se        <- ple$se
  result$egger_pval      <- ple$pval
  if (k == 1){
    bind_result <- result
  } else {
    bind_result <- rbind(bind_result,
                         result
    )
  }
}

Ffilter = 10
dat <- instrument
dat$eaf.exposure <- Outcome[Outcome$SNP %in% instrument$SNP,]$eaf.outcome
N=dat[1,"samplesize.exposure"]
dat=transform(dat,R2=2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))
dat=transform(dat,F=(N-2)*R2/(1-R2))
outTab=dat[dat$F>Ffilter,]