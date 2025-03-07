# Insert the data as dynastic_null, population_null and none_null before run the code

library(ggplot2)
library(gridExtra)
library(dplyr)
library(showtext)
library(cowplot)
library(tidyverse)

dynastic_null   <- dynastic_null[dynastic_null$rou == 0.3,]
population_null <- population_null[population_null$rou == 0.3,]
none_null       <- none_null[none_null$rou == 0.3,]

bias_plot <- function(data,title){
  bias <- ggplot(data, aes(x = Sample size, y = Offspring_effect_Bias, group = Method,
                           shape = factor(Method, levels = c("FT-SEM","lm_parent","lm")),
                           linetype = factor(Method, levels = c("FT-SEM","lm_parent","lm")))) +
    geom_point(size = 2) +  
    geom_line(linewidth = 0.8) +
    labs(x = "Sample size", y = "Bias",shape = "Method",linetype = "Method") +
    # theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), 
      panel.background = element_rect(fill = "white", color = NA),  
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(family = "sans",size = 9), 
      axis.text = element_text(family = "sans",size = 9),  
      strip.text = element_text(family = "sans",size = 9), 
      aspect.ratio = 1,  
      plot.title = element_text(hjust = 0.5, family = "sans", face = "bold", size = 9), 
    ) +
    scale_x_continuous(breaks = c(1000, 2000, 3000)) + 
    scale_y_continuous(limits = c(-0.005,0.15)) +
    ggtitle(title)
  bias
}
bias_null_dynastic   <- bias_plot(dynastic_null, "Dynastic effect")
bias_null_population <- bias_plot(population_null, "Residual population stratification")
bias_null_none       <- bias_plot(none_null, "No bias")

RMSE_plot <- function(data,title){
  RMSE <- ggplot(data, aes(x = Sample size, y = Offspring_effect_RMSE, group = Method,
                           shape = factor(Method, levels = c("FT-SEM","lm_parent","lm")),
                           linetype = factor(Method, levels = c("FT-SEM","lm_parent","lm")))) +
    geom_point(size = 2) +  
    geom_line(linewidth = 0.8) +
    labs(x = "Sample size", y = "RMSE",shape = "Method",linetype = "Method") +
    # theme_minimal() +
    theme(
      panel.grid = element_blank(), 
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  
      panel.background = element_rect(fill = "white", color = NA),  
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(family = "sans",size = 9), 
      axis.text = element_text(family = "sans",size = 9),  
      strip.text = element_text(family = "sans",size = 9), 
      aspect.ratio = 1,  
      plot.title = element_text(hjust = 0.5, family = "sans", face = "bold", size = 9), 
    ) +
    scale_x_continuous(breaks = c(1000, 2000, 3000)) + 
    scale_y_continuous(limits = c(-0.005,0.15)) +
    ggtitle(title)
  RMSE
}
RMSE_null_dynastic   <- RMSE_plot(dynastic_null, "Dynastic effect")
RMSE_null_population <- RMSE_plot(population_null, "Residual population stratification")
RMSE_null_none       <- RMSE_plot(none_null, "No bias")

Coverage_plot <- function(data,title){
  Coverage <- ggplot(data, aes(x = Sample size, y = Offspring_effect_Coverage, group = Method,
                               shape = factor(Method, levels = c("FT-SEM","lm_parent","lm")),
                               linetype = factor(Method, levels = c("FT-SEM","lm_parent","lm")))) +
    geom_point(size = 2) +  
    geom_line(linewidth = 0.8) +
    labs(x = "Sample size", y = "CP of CI",shape = "Method",linetype = "Method") +
    # theme_minimal() +
    theme(
      panel.grid = element_blank(), 
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  
      panel.background = element_rect(fill = "white", color = NA), 
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(family = "sans",size = 9),
      axis.text = element_text(family = "sans",size = 9),  
      strip.text = element_text(family = "sans",size = 9), 
      aspect.ratio = 1,  
      plot.title = element_text(hjust = 0.5, family = "sans", face = "bold", size = 9), 
    ) +
    scale_x_continuous(breaks = c(1000, 2000, 3000)) + 
    scale_y_continuous(limits = c(-0.005,1)) +
    ggtitle(title)
  Coverage
}
Coverage_null_dynastic   <- Coverage_plot(dynastic_null, "Dynastic effect")
Coverage_null_population <- Coverage_plot(population_null, "Residual population stratification")
Coverage_null_none       <- Coverage_plot(none_null, "No bias")

CI_plot <- function(data,title){
  CI <- ggplot(data, aes(x = Sample size, y = Offspring_effect_CILength, group = Method,
                         shape = factor(Method, levels = c("FT-SEM","lm_parent","lm")),
                         linetype = factor(Method, levels = c("FT-SEM","lm_parent","lm")))) +
    geom_point(size = 2) +  
    geom_line(linewidth = 0.8) +
    labs(x = "Sample size", y = "Width of CI",shape = "Method",linetype = "Method") +
    # theme_minimal() +
    theme(
      panel.grid = element_blank(), 
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  
      panel.background = element_rect(fill = "white", color = NA),  
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(family = "sans",size = 9), 
      axis.text = element_text(family = "sans",size = 9),  
      strip.text = element_text(family = "sans",size = 9), 
      aspect.ratio = 1,  
      plot.title = element_text(hjust = 0.5, family = "sans", face = "bold", size = 9),  
    ) +
    scale_x_continuous(breaks = c(1000, 2000, 3000)) + 
    scale_y_continuous(limits = c(-0.005,0.5)) +
    ggtitle(title)
  CI
}
CI_null_dynastic   <- CI_plot(dynastic_null, "Dynastic effect")
CI_null_population <- CI_plot(population_null, "Residual population stratification")
CI_null_none       <- CI_plot(none_null, "No bias")

Type_one_error_rate_plot <- function(data,title){
  Type_one_error_rate <- ggplot(data, aes(x = Sample size, y = Offspring_effect_Type I error rates, group = Method,
                                          shape = factor(Method, levels = c("FT-SEM","lm_parent","lm")),
                                          linetype = factor(Method, levels = c("FT-SEM","lm_parent","lm")))) +
    geom_point(size = 2) +  
    geom_line(linewidth = 0.8) +
    labs(x = "Sample size", y = "Type I error rate",shape = "Method",linetype = "Method") +
    # theme_minimal() +
    theme(
      panel.grid = element_blank(), 
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), 
      panel.background = element_rect(fill = "white", color = NA),  
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(family = "sans",size = 9), 
      axis.text = element_text(family = "sans",size = 9),  
      strip.text = element_text(family = "sans",size = 9), 
      aspect.ratio = 1,  
      plot.title = element_text(hjust = 0.5, family = "sans", face = "bold", size = 9),  # 设置标题居中，字体样式和大小
    ) +
    scale_x_continuous(breaks = c(1000, 2000, 3000)) + 
    scale_y_continuous(limits = c(-0.005,1)) +
    ggtitle(title)
  Type_one_error_rate
}
Type_one_error_rate_null_dynastic   <- Type_one_error_rate_plot(dynastic_null, "Dynastic effect")
Type_one_error_rate_null_population <- Type_one_error_rate_plot(population_null, "Residual population stratification")
Type_one_error_rate_null_none       <- Type_one_error_rate_plot(none_null, "No bias")

library(cowplot)


grid <- plot_grid(
  bias_null_dynastic + theme(legend.position="none"),
  bias_null_population + theme(legend.position="none"),
  bias_null_none + theme(legend.position="none"),
  RMSE_null_dynastic + theme(legend.position="none"),
  RMSE_null_population + theme(legend.position="none"),
  RMSE_null_none + theme(legend.position="none"),
  Coverage_null_dynastic + theme(legend.position="none"),
  Coverage_null_population + theme(legend.position="none"),
  Coverage_null_none + theme(legend.position="none"),
  CI_null_dynastic + theme(legend.position="none"),
  CI_null_population + theme(legend.position="none"),
  CI_null_none + theme(legend.position="none"),
  Type_one_error_rate_null_dynastic + theme(legend.position="none"),
  Type_one_error_rate_null_population + theme(legend.position="none"),
  Type_one_error_rate_null_none + theme(legend.position="none"),
  # labels = c("a","b","c",
  #            "d","e","f",
  #            "g","h","i",
  #            "j","k","l",
  #            "m","n","o"),
  labels = "auto",
  label_fontface = "bold",
  label_size = 14,
  hjust = 0.5,
  label_x = 0.08,
  ncol = 3
  
  
)

legend <- get_legend(bias_null_dynastic +
                       theme(legend.text = element_text(size = 9),
                             legend.title = element_text(size = 12)))

plot_grid(grid, legend, rel_widths = c(3, .4))
