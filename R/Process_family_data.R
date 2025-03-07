#' @title The main function for family trios structural equation model (FT-SEM)
#' @description A robust and powerful GWAS method for family trios
#' @param data A data.frame about genotype from plink1.9, containing FID,IID,PAT,MAT,SEX, phnotype and genotype (only 1 coloum) 
#' @param seed random seed for select offspring from a family with multiple offsprings
#' @import dplyr
#' @export


process_family_data <- function(data, seed = 123) {
  names(data)  <- c("FID","IID","PAT","MAT","SEX","y","SNP")
  
  offspring_data <- data[data$PAT != 0 & data$MAT != 0, ]
  
  set.seed(seed)  
  offspring_sampled <- offspring_data %>%
    group_by(FID) %>%
    sample_n(1) %>%
    ungroup()
  
  father_data <- data[data$IID %in% offspring_sampled$PAT, ]
  mother_data <- data[data$IID %in% offspring_sampled$MAT, ]
  
  data_sub <- merge(offspring_sampled, father_data, by.x = "PAT", by.y = "IID", suffixes = c("_o", "_f"))
  data_sub <- merge(data_sub, mother_data, by.x = "MAT_o", by.y = "IID", suffixes = c("", "_m"))
  
  data_sub <- data_sub[, c("SNP_f", "SNP", "SNP_o",
                           "y_f", "y", "y_o")]
  
  colnames(data_sub) <- c("f_snp", "m_snp", "o_snp", "Exposure_f", "Exposure_m", "Exposure_o")
  
  return(data_sub)
}

