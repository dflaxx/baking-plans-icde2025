## Packages --------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggallin)
library(scales)

# Working directory = dir of this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# -----------------------------------------------------------------------------
# Data Load

estimator_names <- c("Crude_base", "Crude_sel", "IA_L", "IA_M", "IA_S", "Join_base", "Join_sel")
data_filenames <- c("cp_crude_base_error.dat", "cp_crude_sel_error.dat", 
                    "cp_ia_l_error.dat", "cp_ia_m_error.dat", "cp_ia_s_error.dat", "cp_join_base_error.dat", "cp_join_sel_error.dat")
stopifnot(length(estimator_names) == length(data_filenames))

# Read all input files,
# make a list of dataframes (https://stackoverflow.com/questions/17499013/how-do-i-make-a-list-of-data-frames)
data_list <- list()
data_full <- data.frame()
for (i in 1:length(estimator_names)) {
  file_data <- read.table(file = data_filenames[i],
                          col.names = c("planclass_size", "ce_error"),
                          header = FALSE, sep = " ", stringsAsFactors = TRUE) %>%
    mutate(planclass_size = as.factor(planclass_size),
           estimator = estimator_names[i])
  data_list[[i]] <- file_data
  data_full <- bind_rows(data_full, file_data)
  
  print(data_list[[i]] %>% group_by(estimator) %>% summarise(err_min = min(ce_error), err_max = max(ce_error)))
}

stopifnot(length(estimator_names) == length(data_list))


# -----------------------------------------------------------------------------
# Plot Estimators

# Computation function for boxplot stats
quantile_probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
get_boxplot_stats <- function(x) {
  # ymin/ymax: lower/upper whisker
  # lower/upper: box start/end
  # middle: median line
  r <- quantile(x, probs = quantile_probs)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
get_outliers <- function(x) {
  subset(x, x < quantile(x, quantile_probs)[1] | quantile(x, quantile_probs)[5] < x)
}

make_boxplot <- function(aData, aBreaks, aEstName) {
  p <- ggplot(aData, aes(x = planclass_size, y = ce_error)) +
    stat_summary(fun.data = get_boxplot_stats, geom = "boxplot",
                 staplewidth = 0.4,
                 linewidth = 0.7)
    stat_summary(fun = get_outliers, geom = "point", size = 0.5,
                 position = position_jitter(width = 0)) +
    scale_y_continuous(trans = pseudolog10_trans,
                       breaks = aBreaks) +
    # labs(title = paste("CE", aEstName)) +
    xlab(label = "Size of plan class") +
    ylab(label = "Cardinality estimation error (p-error)") +
    theme_bw() +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5)
    theme(
      panel.grid.major.y = element_line(color = "grey75", size = 0.25),
      panel.grid.minor.y = element_line(color = "grey90", size = 0.15),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank() ,
      axis.text = element_text(color="black")
    ) +
    expand_limits(y = c(10, -10))
}

for (est in unique(data_full$estimator)) {
  print(est)
  cur_data <- data_full %>% filter(estimator == est)
  p <- make_boxplot(cur_data,
                    c(-1e12, -1e15, -1e9, -1e6, -1000, -100, -10,
                      0,
                      10, 100, 1000, 1e6, 1e9, 1e12, 1e15, 1e21, 1e27, 1e33, 1e40, 1e44, 1e47),
                    est)
  print(p)
  ggsave(filename = paste("plot_", est, ".pdf", sep = ""),
         unit = "mm",
         scale = 0.5)
}
