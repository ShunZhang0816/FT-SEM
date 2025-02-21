# Insert the data as dynastic_null, population_null and none_null before run the code

library(ggplot2)
library(gridExtra)
library(dplyr)
library(showtext)
library(cowplot)
library(tidyverse)

dynastic_null   <- dynastic_null[dynastic_null$残差相关性 == 0.3,]
population_null <- population_null[population_null$残差相关性 == 0.3,]
none_null       <- none_null[none_null$残差相关性 == 0.3,]

bias_plot <- function(data,title){
  bias <- ggplot(data, aes(x = 样本量, y = 直接效应Bias, group = 方法,
                           shape = factor(方法, levels = c("FT-SEM","lm_parent","lm")),
                           linetype = factor(方法, levels = c("FT-SEM","lm_parent","lm")))) +
    geom_point(size = 2) +  
    geom_line(linewidth = 0.8) +
    labs(x = "Sample size", y = "Bias",shape = "Method",linetype = "Method") +
    # theme_minimal() +
    theme(
      panel.grid = element_blank(), # 去除背景的方框线
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # 在面板周围添加黑色实线边框
      panel.background = element_rect(fill = "white", color = NA),  # 设置面板背景为白色
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(family = "sans",size = 9), # 加粗并增大轴标题字体
      axis.text = element_text(family = "sans",size = 9),  # 增大坐标轴刻度字体
      strip.text = element_text(family = "sans",size = 9), # 增大 facet 标题的字体大小并加粗
      aspect.ratio = 1,  # 设置图形为正方形
      plot.title = element_text(hjust = 0.5, family = "sans", face = "bold", size = 9),  # 设置标题居中，字体样式和大小
    ) +
    scale_x_continuous(breaks = c(1000, 2000, 3000)) + # 设置横坐标刻度为1000、2000和3000
    scale_y_continuous(limits = c(-0.005,0.15)) +
    ggtitle(title)
  bias
}
bias_null_dynastic   <- bias_plot(dynastic_null, "Dynastic effect")
bias_null_population <- bias_plot(population_null, "Residual population stratification")
bias_null_none       <- bias_plot(none_null, "No bias")

RMSE_plot <- function(data,title){
  RMSE <- ggplot(data, aes(x = 样本量, y = 直接效应RMSE, group = 方法,
                           shape = factor(方法, levels = c("FT-SEM","lm_parent","lm")),
                           linetype = factor(方法, levels = c("FT-SEM","lm_parent","lm")))) +
    geom_point(size = 2) +  # 添加点
    geom_line(linewidth = 0.8) +
    labs(x = "Sample size", y = "RMSE",shape = "Method",linetype = "Method") +
    # theme_minimal() +
    theme(
      panel.grid = element_blank(), # 去除背景的方框线
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # 在面板周围添加黑色实线边框
      panel.background = element_rect(fill = "white", color = NA),  # 设置面板背景为白色
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(family = "sans",size = 9), # 加粗并增大轴标题字体
      axis.text = element_text(family = "sans",size = 9),  # 增大坐标轴刻度字体
      strip.text = element_text(family = "sans",size = 9), # 增大 facet 标题的字体大小并加粗
      aspect.ratio = 1,  # 设置图形为正方形
      plot.title = element_text(hjust = 0.5, family = "sans", face = "bold", size = 9),  # 设置标题居中，字体样式和大小
    ) +
    scale_x_continuous(breaks = c(1000, 2000, 3000)) + # 设置横坐标刻度为1000、2000和3000
    scale_y_continuous(limits = c(-0.005,0.15)) +
    ggtitle(title)
  RMSE
}
RMSE_null_dynastic   <- RMSE_plot(dynastic_null, "Dynastic effect")
RMSE_null_population <- RMSE_plot(population_null, "Residual population stratification")
RMSE_null_none       <- RMSE_plot(none_null, "No bias")

Coverage_plot <- function(data,title){
  Coverage <- ggplot(data, aes(x = 样本量, y = 直接效应Coverage, group = 方法,
                               shape = factor(方法, levels = c("FT-SEM","lm_parent","lm")),
                               linetype = factor(方法, levels = c("FT-SEM","lm_parent","lm")))) +
    geom_point(size = 2) +  # 添加点
    geom_line(linewidth = 0.8) +
    labs(x = "Sample size", y = "CP of CI",shape = "Method",linetype = "Method") +
    # theme_minimal() +
    theme(
      panel.grid = element_blank(), # 去除背景的方框线
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # 在面板周围添加黑色实线边框
      panel.background = element_rect(fill = "white", color = NA),  # 设置面板背景为白色
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(family = "sans",size = 9), # 加粗并增大轴标题字体
      axis.text = element_text(family = "sans",size = 9),  # 增大坐标轴刻度字体
      strip.text = element_text(family = "sans",size = 9), # 增大 facet 标题的字体大小并加粗
      aspect.ratio = 1,  # 设置图形为正方形
      plot.title = element_text(hjust = 0.5, family = "sans", face = "bold", size = 9),  # 设置标题居中，字体样式和大小
    ) +
    scale_x_continuous(breaks = c(1000, 2000, 3000)) + # 设置横坐标刻度为1000、2000和3000
    scale_y_continuous(limits = c(-0.005,1)) +
    ggtitle(title)
  Coverage
}
Coverage_null_dynastic   <- Coverage_plot(dynastic_null, "Dynastic effect")
Coverage_null_population <- Coverage_plot(population_null, "Residual population stratification")
Coverage_null_none       <- Coverage_plot(none_null, "No bias")

CI_plot <- function(data,title){
  CI <- ggplot(data, aes(x = 样本量, y = 直接效应CILength, group = 方法,
                         shape = factor(方法, levels = c("FT-SEM","lm_parent","lm")),
                         linetype = factor(方法, levels = c("FT-SEM","lm_parent","lm")))) +
    geom_point(size = 2) +  # 添加点
    geom_line(linewidth = 0.8) +
    labs(x = "Sample size", y = "Width of CI",shape = "Method",linetype = "Method") +
    # theme_minimal() +
    theme(
      panel.grid = element_blank(), # 去除背景的方框线
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # 在面板周围添加黑色实线边框
      panel.background = element_rect(fill = "white", color = NA),  # 设置面板背景为白色
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(family = "sans",size = 9), # 加粗并增大轴标题字体
      axis.text = element_text(family = "sans",size = 9),  # 增大坐标轴刻度字体
      strip.text = element_text(family = "sans",size = 9), # 增大 facet 标题的字体大小并加粗
      aspect.ratio = 1,  # 设置图形为正方形
      plot.title = element_text(hjust = 0.5, family = "sans", face = "bold", size = 9),  # 设置标题居中，字体样式和大小
    ) +
    scale_x_continuous(breaks = c(1000, 2000, 3000)) + # 设置横坐标刻度为1000、2000和3000
    scale_y_continuous(limits = c(-0.005,0.5)) +
    ggtitle(title)
  CI
}
CI_null_dynastic   <- CI_plot(dynastic_null, "Dynastic effect")
CI_null_population <- CI_plot(population_null, "Residual population stratification")
CI_null_none       <- CI_plot(none_null, "No bias")

Type_one_error_rate_plot <- function(data,title){
  Type_one_error_rate <- ggplot(data, aes(x = 样本量, y = 直接效应第一类错误率, group = 方法,
                                          shape = factor(方法, levels = c("FT-SEM","lm_parent","lm")),
                                          linetype = factor(方法, levels = c("FT-SEM","lm_parent","lm")))) +
    geom_point(size = 2) +  # 添加点
    geom_line(linewidth = 0.8) +
    labs(x = "Sample size", y = "Type I error rate",shape = "Method",linetype = "Method") +
    # theme_minimal() +
    theme(
      panel.grid = element_blank(), # 去除背景的方框线
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # 在面板周围添加黑色实线边框
      panel.background = element_rect(fill = "white", color = NA),  # 设置面板背景为白色
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(family = "sans",size = 9), # 加粗并增大轴标题字体
      axis.text = element_text(family = "sans",size = 9),  # 增大坐标轴刻度字体
      strip.text = element_text(family = "sans",size = 9), # 增大 facet 标题的字体大小并加粗
      aspect.ratio = 1,  # 设置图形为正方形
      plot.title = element_text(hjust = 0.5, family = "sans", face = "bold", size = 9),  # 设置标题居中，字体样式和大小
    ) +
    scale_x_continuous(breaks = c(1000, 2000, 3000)) + # 设置横坐标刻度为1000、2000和3000
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
# 将图例和之前的图形进行组合，并设置宽度比例
plot_grid(grid, legend, rel_widths = c(3, .4))