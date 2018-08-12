# Summarize the output from single regression per state analyses.
library(tidyverse)
library(MASS)


ns <- c(10, 25, 50) #sample size for training set
states <- c(2,3,27)
v_inflate <- c(5, 10, 50, 100)


logical_ytibble <- function(y_open, colname){
  #works on y_open, y_type1 logical matrices
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
  
  if(MargName %in% c("margNTheory","margRest", "margRestHuber", "margT")){
  tmp_array <- pooled_results[[MargName]]
  v_tibbles <- apply(tmp_array, 2, function(tmp){
  tmp <- as.matrix(t(tmp))
  rownames(tmp) <- 1:nrow(tmp)
  tmp <- as_tibble(tmp, rownames = 'holdout sample')
  tmp <- tmp %>% 
    rename_at(vars(contains('V')), funs(sub('V', '', .))) %>% 
    gather('Repetition', 'Marginal', -"holdout sample") %>% 
    mutate(Model = ModelName, Repetition = as.numeric(Repetition), `holdout sample` = as.numeric(`holdout sample`))
  
  y_open <- logical_ytibble(pooled_results$y_open, colname = 'Open')
  y_type1 <- logical_ytibble(pooled_results$y_type1, colname = 'Type1')
  
  full_join(tmp, y_open,  by = c("holdout sample", "Repetition")) %>% full_join(., y_type1,  by = c("holdout sample", "Repetition"))
  })
  names(v_tibbles) <- v_inflate
  bind_rows(v_tibbles, .id = 'var_inflate')
  } else {
    tmp <- t(pooled_results[[MargName]])
    rownames(tmp) <- 1:nrow(tmp)
    tmp <- as_tibble(tmp, rownames = 'holdout sample')
    tmp <- tmp %>% 
      rename_at(vars(contains('V')), funs(sub('V', '', .))) %>% 
      gather('Repetition', 'Marginal', -"holdout sample") %>% 
      mutate(Model = ModelName, Repetition = as.numeric(Repetition), `holdout sample` = as.numeric(`holdout sample`))
    
    y_open <- logical_ytibble(pooled_results$y_open, colname = 'Open')
    y_type1 <- logical_ytibble(pooled_results$y_type1, colname = 'Type1')
    
    full_join(tmp, y_open,  by = c("holdout sample", "Repetition")) %>% full_join(., y_type1,  by = c("holdout sample", "Repetition"))
  }
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


pooled_marginals_all <- function(ns, states){
    list_1 <- lapply(states, function(State_keep){  
      tbl_list <- lapply(ns, function(n){
        pooled_results <- read_rds(paste0('single_reg_state_', State_keep, '_n_', n, '.rds' ))
       marginals_tibble(pooled_results) %>% 
       mutate(n = n, State = State_keep)
  })
   bind_rows(tbl_list)
})
bind_rows(list_1)
}

pooled_marginals <- pooled_marginals_all(ns, states)
anyNA(pooled_marginals)
pooled_marginals <- pooled_marginals %>% 
  mutate(n = as.factor(n), State = as.factor(State))
# var_inflate = as.factor(var_inflate)
#add artificial levels of v-inflate for the classical models

p0 <- pooled_marginals %>% filter(is.na(var_inflate)) %>% dplyr::select(-var_inflate) %>% mutate(var_inflate = "5")
p1 <- pooled_marginals %>% filter(is.na(var_inflate)) %>% dplyr::select(-var_inflate) %>% mutate(var_inflate = "10")
p2 <- pooled_marginals %>% filter(is.na(var_inflate)) %>% dplyr::select(-var_inflate) %>% mutate(var_inflate = "50")
p3 <- pooled_marginals %>% filter(is.na(var_inflate)) %>% dplyr::select(-var_inflate) %>% mutate(var_inflate = "100")
pooled_marginals <- bind_rows(pooled_marginals %>% filter(!is.na(var_inflate)), p0, p1, p2, p3) 
# Trimmed log marginal distribution 
# tlm <- function(marginals, fun = mean, trimming_fraction = 0.3){
#   log_marg <- log(marginals)
#   include <- which(log_marg >= quantile(log_marg, trimming_fraction))
#   fun(log_marg[include])
# }
#tlm(pooled_results$margNTheory[1,])

# Finding TLM with a specific base model -----


base_model <- 'Student-t'
trimming_fraction <- 0.3

marginals_open_type1 <- pooled_marginals %>% 
  filter(Open == TRUE & Type1 == TRUE)

base_model_marg <-  marginals_open_type1 %>% 
  filter(Model == base_model) 

lower_tlm_vals <- base_model_marg %>% 
  group_by(Repetition, n, State, var_inflate) %>% 
  summarise(lower_tlm = quantile(Marginal, trimming_fraction)) 


#for each repetition, find the holdout sample indexes to include from the base model----

holdout_include <- left_join(base_model_marg,lower_tlm_vals, by = c('Repetition', 'n', "State","var_inflate")) %>% 
  mutate(include = Marginal > lower_tlm) %>% 
  dplyr::select(`holdout sample`, Repetition, include, n, State, var_inflate, lower_tlm) %>% 
  filter(include == TRUE)

include_by_rep <- holdout_include %>%
  split(list(.$Repetition, .$n, .$State, .$var_inflate)) %>% map(function(x) x$`holdout sample`)


#include_by_rep: has the holdout sample included for tlm calculation for each rep. 
#for each Repetition value, marginals_open_type_1 will have the sample houldout samples for each model. 
marg_split <- marginals_open_type1 %>% 
  split(list(.$Repetition, .$n, .$State, .$var_inflate)) %>% map2(.x = ., .y = include_by_rep, .f = function(x, y) filter(x,  x$`holdout sample` %in% y)) %>% 
  bind_rows()


summary_tibble <- marg_split %>% 
  group_by(Repetition, Model, n, State, var_inflate) %>% 
  summarise(tlm_mean = mean(Marginal), tlm_sd = sd(Marginal)) %>% 
  ungroup() %>% 
  group_by(Model, n,State, var_inflate) %>% 
  summarise(mean = mean(tlm_mean), 
            sd = sd(tlm_mean),
            mean_sd = mean(tlm_sd),
            sd_sd = sd(tlm_sd)) %>% 
  ungroup() %>% 
  mutate(Model = factor(Model, levels = c( 'Restricted - Huber','Rlm - Huber', 'Restricted - Tukey', 'Rlm - Tukey','Student-t', 'Normal', 'OLS')), n = factor(n), var_inflate = factor(var_inflate, levels = v_inflate))



theme_set(theme_bw(base_family = 'Times'))
ggplot(filter(summary_tibble, !Model %in% c('OLS', 'Normal')), aes(x = n, y = mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + #geom_line(position = position_dodge(width = .5)) +
  geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), width = 0.05, position  = position_dodge(width = .5), linetype = 'dotted') +
  facet_wrap(~State+var_inflate, drop = FALSE)



ggplot(filter(summary_tibble, !Model %in% c('OLS', 'Normal')), aes(x = n, y = mean_sd/mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + #geom_line(position = position_dodge(width = .5)) +
  geom_errorbar(mapping = aes(ymin = (mean_sd - sd_sd)/mean, ymax = (mean_sd + sd_sd)/mean), width = 0.05, position  = position_dodge(width = .5), linetype = 'dotted') +
  facet_wrap(~State+var_inflate, drop = FALSE, scales = 'free')

