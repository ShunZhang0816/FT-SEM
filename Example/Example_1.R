library(FTSEM)
# install.packages("data.table")
library(data.table)
data <- read.table("Example_1.txt")
data <- as.data.table(data)

for(i in 1:(ncol(data) - 6)){
  snp_name  <- names(data[,-(1:6)])
  data_sub  <- data[,c(1:6,(6+i)),with = F]
  trio_data <- process_family_data(data_sub, seed = 0)
  
  
  if (i == 1){
    result    <- FT_SEM(trio_data,snp_name[i])
  } else {
    result <- as.data.frame(rbind(result,FT_SEM(trio_data,snp_name[i])))
  }
}
result
