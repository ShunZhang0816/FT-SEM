#' @title The main function for family trios structural equation model (FT-SEM)
#' @description A robust and powerful GWAS method for family trios
#' @param data A data.frame about genotype from plink1.9, containing FID,IID,PAT,MAT,SEX, phnotype and genotype (only 1 coloum) 
#' @param seed random seed for select offspring from a family with multiple offsprings
#' @export


process_family_data <- function(data, seed = 123) {
  library(dplyr)
  names(data)  <- c("FID","IID","PAT","MAT","SEX","y","SNP")
  
  # 筛选子代个体（PAT 和 MAT 不等于 0）
  offspring_data <- data[data$PAT != 0 & data$MAT != 0, ]
  
  # 随机选择每个家系 (FID) 内的一个子代
  set.seed(seed)  # 设置随机种子，保证结果可复现
  offspring_sampled <- offspring_data %>%
    group_by(FID) %>%
    sample_n(1) %>%
    ungroup()
  
  # 识别父母个体
  father_data <- data[data$IID %in% offspring_sampled$PAT, ]
  mother_data <- data[data$IID %in% offspring_sampled$MAT, ]
  
  # 按 FID 和 IID 匹配家庭成员
  data_sub <- merge(offspring_sampled, father_data, by.x = "PAT", by.y = "IID", suffixes = c("_o", "_f"))
  data_sub <- merge(data_sub, mother_data, by.x = "MAT_o", by.y = "IID", suffixes = c("", "_m"))
  
  # 选择所需列，并重命名以符合 FT-SEM 需求
  data_sub <- data_sub[, c("SNP_f", "SNP", "SNP_o",
                           "y_f", "y", "y_o")]
  
  # 统一列名
  colnames(data_sub) <- c("f_snp", "m_snp", "o_snp", "Exposure_f", "Exposure_m", "Exposure_o")
  
  return(data_sub)
}

