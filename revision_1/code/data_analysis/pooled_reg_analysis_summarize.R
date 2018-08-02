# Summarize the output from pooled regression analysis.
library(tidyverse)
library(MASS)
ns <- c(25, 100) #, 1000)

# pooled_results <- readRDS("~/Dropbox/school/osu/dissertationResearch/snm/journal_papers/bayes_rest_like_methods/revision_1/code/data_analysis/pooled_reg_n25.rds")


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


single_marginal_tibble <- function(pooled_results, MargName, ModelName){
  # given pooled_results and the marginal name, creates tibble with holdoutsample, repetition number, marginal value, model name and indicator columns for open and type 1 agencies:
  #holdout sample is nested within repetition. 
  
tmp <- t(pooled_results[[MargName]])
rownames(tmp) <- 1:nrow(tmp)
tmp <- as_tibble(tmp, rownames = 'holdout sample')
tmp <- tmp %>% 
  rename_at(vars(contains('V')), funs(sub('V', '', .))) %>% 
  gather('Repetition', 'Marginal', -"holdout sample") %>% 
  mutate(Model = ModelName, Repetition = as.numeric(Repetition), `holdout sample` = as.numeric(`holdout sample`))

y_open <- logical_ytibble(pooled_results$y_open, colname = 'Open')
y_type1 <- logical_ytibble(pooled_results$y_type1, colname = 'Type1')

full_join(tmp, y_open) %>% full_join(., y_type1)

}


marginals_tibble <- function(pooled_results){
  #pooled results is the list of output for a set of runs from the pooled regression - there is one for each sample size. 
  # function takes marginals from each model and each run and applies single_marginal_tibble and combines the results
t1 <- single_marginal_tibble(pooled_results,"margOls",'OLS') 
t2 <- single_marginal_tibble(pooled_results, "margRlm", 'Rlm - Tukey')  
t3 <- single_marginal_tibble(pooled_results, "margRlmHuber", 'Rlm - Huber')  
t4 <- single_marginal_tibble(pooled_results, "margNTheory", 'Normal')   
t5 <- single_marginal_tibble(pooled_results, "margRest", 'Restricted - Tukey')  
t6 <- single_marginal_tibble(pooled_results, "margRestHuber", 'Restricted - Huber')  
t7 <- single_marginal_tibble(pooled_results, "margT", 'Student-t')  
bind_rows(t1,t2,t3,t4,t5,t6,t7)
}


pooled_marginals_all <- function(ns){
  #ns is a vector of sample size. this will combine the results of marginals_tibble for all sample sizes in ns. 
  tbl_list <- lapply(ns, function(n){
    pooled_results <- read_rds(paste0("pooled_reg_n", n, ".rds"))
    marginals_tibble(pooled_results) %>% 
      mutate(n = n)
})
bind_rows(tbl_list)
}

pooled_marginals <- pooled_marginals_all(ns)
anyNA(pooled_marginals)

# Trimmed log marginal distribution 
# tlm <- function(marginals, fun = mean, trimming_fraction = 0.3){
#   log_marg <- log(marginals)
#   include <- which(log_marg >= quantile(log_marg, trimming_fraction))
#   fun(log_marg[include])
# }
#tlm(pooled_results$margNTheory[1,])

# Finding TLM with a specific base model -----


base_model <- 'Rlm - Huber'
trimming_fraction <- 0.3

marginals_open_type1 <- pooled_marginals %>% 
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



# pooled2 <- left_join(marginals_open_type1 ,lower_tlm_vals, by = c('Repetition', 'n'))

summary_tibble <- marg_split %>% 
 group_by(Repetition, Model, n) %>% 
summarise(tlm_mean = mean(Marginal), tlm_sd = sd(Marginal)) %>% 
  ungroup() %>% 
  group_by(Model, n) %>% 
  summarise(mean = mean(tlm_mean), 
            sd = sd(tlm_mean),
            median = median(tlm_mean), 
            lower = quantile(tlm_mean, 0.1),
            upper = quantile(tlm_mean, 0.9)) %>% 
  ungroup() %>% 
  mutate(Model = factor(Model, levels = c( 'Restricted - Huber','Rlm - Huber', 'Restricted - Tukey', 'Rlm - Tukey','Student-t', 'Normal', 'OLS')), n = factor(n))







 theme_set(theme_bw(base_family = 'Times'))
# ggplot(filter(summary_tibble, !Model %in% c('Normal', 'OLS')), aes(x = n, y = median, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  +# geom_line(position = position_dodge(width = .5)) +
#   geom_errorbar(mapping = aes(ymin = lower, ymax = upper), width = 0.05, position  = position_dodge(width = .5), linetype = 'dotted') 



ggplot(filter(summary_tibble, !Model %in% c('Normal', 'OLS')), aes(x = n, y = mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + #geom_line(position = position_dodge(width = .5)) +
  geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), width = 0.05, position  = position_dodge(width = .5), linetype = 'dotted') 


