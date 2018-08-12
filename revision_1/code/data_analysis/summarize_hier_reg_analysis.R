# Summarize the output from pooled regression analysis.
library(tidyverse)
library(MASS)

ns <- c(1000,2000)

logical_ytibble <- function(y_open, colname){
  #works on y_open, y_type1 logical matrices
  #creates a tibble
  y_open <- t(y_open)
  rownames(y_open) <- 1:nrow(y_open)
  y_open <- as_tibble(y_open, rownames = 'holdout sample')
  
  y_open %>% 
    rename_at(vars(contains('V')), funs(sub('V', '', .))) %>% 
    gather('Repetition',!!colname, -"holdout sample") %>% 
    mutate(Repetition = as.numeric(Repetition), `holdout sample` = as.numeric(`holdout sample`))
}


get_marginals <- function(hier_results, n){
y_open <- logical_ytibble(hier_results$y_open, 'Open')
y_type1 <- logical_ytibble(hier_results$y_type1, 'Type1')

margs <- apply(hier_results$marginals, 1, function(marg) {
  marg <- t(marg)
  rownames(marg) <- 1:nrow(marg)
  marg <- as_tibble(marg, rownames = 'holdout sample')
  marg <- marg %>% 
    rename_at(vars(contains('V')), funs(sub('V', '', .))) %>% 
    gather('Repetition','Marginal', -"holdout sample") %>% 
    mutate(Repetition = as.numeric(Repetition), `holdout sample` = as.numeric(`holdout sample`))
  
full_join(marg, y_open, by = c("holdout sample", "Repetition")) %>% full_join(y_type1, by = c("holdout sample", "Repetition"))
  
})
margs_sd <- apply(hier_results$marginals_sd, 1, function(marg) {
  marg <- t(marg)
  rownames(marg) <- 1:nrow(marg)
  marg <- as_tibble(marg, rownames = 'holdout sample')
  marg <- marg %>% 
    rename_at(vars(contains('V')), funs(sub('V', '', .))) %>% 
    gather('Repetition','SD', -"holdout sample") %>% 
    mutate(Repetition = as.numeric(Repetition), `holdout sample` = as.numeric(`holdout sample`))
  
full_join(marg, y_open, by = c("holdout sample", "Repetition")) %>% full_join(y_type1, by = c("holdout sample", "Repetition"))
  
})


names(margs) <-names(margs_sd) <- c("OLS",                "Rlm - Tukey",       
                  "Rlm - Huber",        "Normal",            
                  "Restricted - Tukey", "Restricted - Huber",
                  "Student-t")
margs <- bind_rows(margs, .id = 'Model')
margs <- margs %>% mutate(n = factor(n, levels = ns))

margs_sd <- bind_rows(margs_sd, .id = 'Model')
margs_sd <- margs_sd %>% mutate(n = factor(n, levels = ns))
full_join(margs, margs_sd)
}

# hier_results <- readRDS(file.path(getwd(), paste0('hier_reg_n', 1000,".rds")))
# plot(hier_results$group_converge)

hier_marginals <- ns %>% map(.f = function(n){
  hier_results <- readRDS(file.path(getwd(), paste0('hier_reg_n', n,".rds")))
  get_marginals(hier_results,n)
}) %>% 
  bind_rows()

# 
# hier_results <- readRDS(file.path(getwd(), paste0('hier_reg_n', 1000,".rds")))

# Finding TLM with a specific base model -----

get_tlm_tibble <- function(marginals_tibble, base_model = 'Student-t', trimming_fraction = 0.3){
marginals_open_type1 <- marginals_tibble %>% 
  filter(Open == TRUE & Type1 == TRUE)

base_model_marg <-  marginals_open_type1 %>% 
  filter(Model == base_model) 

lower_tlm_vals <- base_model_marg %>% 
  group_by(Repetition, n) %>% 
  summarise(lower_tlm = quantile(Marginal, trimming_fraction)) 


#for each repetition, find the holdout sample indexes to include from the base model----

holdout_include <- left_join(base_model_marg,lower_tlm_vals, by = c('Repetition', 'n')) %>% 
  mutate(include = Marginal > lower_tlm) %>% 
  dplyr::select(`holdout sample`, Repetition, include, n) %>% 
  filter(include == TRUE)

include_by_rep <- holdout_include %>%
  split(list(.$Repetition, .$n)) %>% 
  map(function(x) x$`holdout sample`)

# map(marg_split, .f = function(yy){
# tst <- yy %>%
#   split(.$Model)
# sapply(tst, function(x) all.equal(x$`holdout sample`, tst[[1]]$`holdout sample`))
# }) %>% bind_rows() %>%
#   apply(., 2, function(a) any(a == FALSE))
#include_by_rep: has the holdout sample included for tlm calculation for each rep. 
#for each Repetition value, marginals_open_type_1 will have the sample houldout samples for each model. 
marg_split <- marginals_open_type1 %>% 
  split(list(.$Repetition, .$n)) %>% 
  map2(.x = ., .y = include_by_rep, .f = function(x, y) filter(x,  x$`holdout sample` %in% y)) %>% 
  bind_rows()

marg_split %>% 
  group_by(Repetition, Model, n) %>% 
  summarise(tlm_mean = mean(Marginal), tlm_sd = sd(Marginal), tlm_mean_sd = mean(SD)) %>% 
  ungroup() %>% 
  group_by(Model, n) %>% 
  summarise(mean = mean(tlm_mean), 
            sd = sd(tlm_mean),
            median = median(tlm_mean), 
            lower = quantile(tlm_mean, 0.1),
            upper = quantile(tlm_mean, 0.9),
            mean_sd = mean(tlm_mean_sd)) %>% 
  ungroup() %>% 
  mutate(Model = factor(Model, levels = c( 'Restricted - Huber','Rlm - Huber', 'Restricted - Tukey', 'Rlm - Tukey','Student-t', 'Normal', 'OLS')), n = factor(n))
}



summary_tibble <- get_tlm_tibble(hier_marginals, base_model = 'Student-t', trimming_fraction = 0.3)


theme_set(theme_bw(base_family = 'Times'))
ggplot(filter(summary_tibble, !Model %in% c('Normal' ,'OLS')), aes(x = n, y = mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + #geom_line(position = position_dodge(width = .5)) +
  geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), width = 0.05, position  = position_dodge(width = .5), linetype = 'dotted') 



ggplot(filter(summary_tibble, !Model %in% c('Normal' ,'OLS')), aes(x = n, y = mean_sd, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + #geom_line(position = position_dodge(width = .5)) +
  geom_errorbar(mapping = aes(ymin = mean_sd - sd, ymax = mean_sd + sd), width = 0.05, position  = position_dodge(width = .5), linetype = 'dotted') 



# ggplot(filter(summary_tibble, !Model %in% c('Normal', 'OLS')), aes(x = n, y = median, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  +# geom_line(position = position_dodge(width = .5)) +
#   geom_errorbar(mapping = aes(ymin = lower, ymax = upper), width = 0.05, position  = position_dodge(width = .5), linetype = 'dotted')

