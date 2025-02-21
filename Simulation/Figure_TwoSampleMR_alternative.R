# Insert the data as dynastic_alt, population_alt and none_alt before run the code

library(ggplot2)
library(gridExtra)
library(dplyr)
library(showtext)
library(cowplot)
library(tidyverse)

dynastic_alt   <- dynastic_alt[dynastic_alt$残差相关性 == 0.3,]
population_alt <- population_alt[population_alt$残差相关性 == 0.3,]
none_alt       <- none_alt[none_alt$残差相关性 == 0.3,]

bias_plot <- function(data,title){
  bias <- ggplot(data, aes(x = 样本量, y = 直接效应Bias, group = 方法,
                           shape = factor(方法, levels = c("FT-SEM","lm_parent","lm")),
                           linetype = factor(方法, levels = c("FT-SEM","lm_parent","lm")))) +
    geom_point(size = 2) +  # 添加点
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
bias_alt_dynastic   <- bias_plot(dynastic_alt, "Dynastic effect")
bias_alt_population <- bias_plot(population_alt, "Residual population stratification")
bias_alt_none       <- bias_plot(none_alt, "No bias")

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
RMSE_alt_dynastic   <- RMSE_plot(dynastic_alt, "Dynastic effect")
RMSE_alt_population <- RMSE_plot(population_alt, "Residual population stratification")
RMSE_alt_none       <- RMSE_plot(none_alt, "No bias")

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
Coverage_alt_dynastic   <- Coverage_plot(dynastic_alt, "Dynastic effect")
Coverage_alt_population <- Coverage_plot(population_alt, "Residual population stratification")
Coverage_alt_none       <- Coverage_plot(none_alt, "No bias")

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
CI_alt_dynastic   <- CI_plot(dynastic_alt, "Dynastic effect")
CI_alt_population <- CI_plot(population_alt, "Residual population stratification")
CI_alt_none       <- CI_plot(none_alt, "No bias")

Power_plot <- function(data,title){
  Power <- ggplot(data, aes(x = 样本量, y = 直接效应Power, group = 方法,
                            shape = factor(方法, levels = c("FT-SEM","lm_parent","lm")),
                            linetype = factor(方法, levels = c("FT-SEM","lm_parent","lm")))) +
    geom_point(size = 2) +  # 添加点
    geom_line(linewidth = 0.8) +
    labs(x = "Sample size", y = "Power",shape = "Method",linetype = "Method") +
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
  Power
}
Power_alt_dynastic   <- Power_plot(dynastic_alt, "Dynastic effect")
Power_alt_population <- Power_plot(population_alt, "Residual population stratification")
Power_alt_none       <- Power_plot(none_alt, "No bias")

library(cowplot)


grid <- plot_grid(
  bias_alt_dynastic + theme(legend.position="none"),
  bias_alt_population + theme(legend.position="none"),
  bias_alt_none + theme(legend.position="none"),
  RMSE_alt_dynastic + theme(legend.position="none"),
  RMSE_alt_population + theme(legend.position="none"),
  RMSE_alt_none + theme(legend.position="none"),
  Coverage_alt_dynastic + theme(legend.position="none"),
  Coverage_alt_population + theme(legend.position="none"),
  Coverage_alt_none + theme(legend.position="none"),
  CI_alt_dynastic + theme(legend.position="none"),
  CI_alt_population + theme(legend.position="none"),
  CI_alt_none + theme(legend.position="none"),
  Power_alt_dynastic + theme(legend.position="none"),
  Power_alt_population + theme(legend.position="none"),
  Power_alt_none + theme(legend.position="none"),
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

legend <- get_legend(bias_alt_dynastic +
                       theme(legend.text = element_text(size = 9),
                             legend.title = element_text(size = 12)))
# 将图例和之前的图形进行组合，并设置宽度比例
plot_grid(grid, legend, rel_widths = c(3, .4))