# summarize simulation 1
library(tidyverse)
library(stringr)

results_dir <- str_subset(dir(getwd()), 'results_true_model')

results_all <- map(results_dir, function(dir){
  tib <- read_rds(file.path(dir, "means_all_long.rds"))
  tib$true_model <- str_sub(dir, 20)
  tib
})

results_all <- bind_rows(results_all)

mse_all <- results_all %>% 
  group_by(true_model, n, `Fitted Model`) %>% 
  summarise(MSE = mean((mean-0)^2), sd_mean = sd(mean), bias = mean(mean), `# simulations` = n())
mse_all$n <- factor(mse_all$n, levels = c("25", "50", "100"))

ggplot(mse_all, aes(x = n, y = MSE, shape = `Fitted Model`, group = `Fitted Model`)) + geom_point(alpha = .5) + geom_line(alpha = .5) + facet_wrap(~true_model, nrow = 2, scales = 'free')

ggplot(mse_all, aes(x = n, y = sd_mean, col = `Fitted Model`, group = `Fitted Model`)) + geom_point(alpha = .5) + geom_line(alpha = .5) + facet_wrap(~true_model, nrow = 2, scales = 'free')


ggplot(mse_all, aes(x = n, y = bias^2 , col = `Fitted Model`, group = `Fitted Model`)) + geom_point(alpha = .5) + geom_line(alpha = .5) + facet_wrap(~true_model, nrow = 2, scales = 'free')

#plot relative MSE's

# make the levels of true_model match those of `Fitted Model`
mse_all$true_model <- factor(sapply(mse_all$true_model, function(x){
  if(x == 'mixture 1') {a <- 'Mixture 1'}
  if(x == 'mixture 2') {a <- 'Mixture 2'}
  if(x == 'normal') {a <- 'Normal'}
  if(x == 'Student-t') {a <- 'Student-t'}
  a
}), levels = c('Normal', 'Student-t', 'Mixture 1', 'Mixture 2'))


relative_mse <- mse_all %>% 
  group_by(true_model, n) %>% 
  mutate(indicator = true_model == `Fitted Model`, base_mse = MSE[indicator], relative_mse = MSE/base_mse)

relative_mse$`Fitted Model` <- factor(relative_mse$`Fitted Model`, levels = c("Normal", "Student-t", "Mixture 1", "Mixture 2",  "Huber","Tukey"))

ggplot(relative_mse, aes(x = n, y = relative_mse, shape = `Fitted Model`, group = `Fitted Model`)) + geom_point(alpha = .5) + geom_line(alpha = .5) + facet_wrap(~true_model, nrow = 2, scales = 'free')
