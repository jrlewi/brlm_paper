# Summarize the output from single regression per state analyses.
library(tidyverse)
library(MASS)


ns <- c(25, 50) #sample size for training set
states <- c(2, 15, 27, 36)
v_inflate <- c(100)

# plots ---- 
analysis_data <- read_rds(file.path(here::here(), 'data', 'analysis_data.rds'))

analysis_data <- analysis_data %>% 
  mutate(sqrt_count_2010 = sqrt(Count_2010), sqrt_count_2012 = sqrt(Count_2012)) %>% filter(State %in% states) %>% 
  mutate(State = droplevels(State))

label_vals <- c('2' = 'State 2', '15' = 'State 15', '27' = 'State 27', '36' = 'State 36')
ggplot(analysis_data) + geom_point(aes(x = sqrt_count_2010, y = sqrt_count_2012, col = Type), size = 1, alpha = .5) +
  theme_bw() + facet_wrap(~State, labeller = labeller(State = label_vals)) + labs(x = 'square root of 2010 houshold count', y = 'square root of 2012 houshold count') +
  theme(text = element_text(family = 'Times'))
ggsave(file.path(getwd(), "..", "..", "figs", 'scatter_by_state.png'), width = 6, height = 4)



cnts <- analysis_data %>% group_by(State) %>% 
  summarise(n = n())
# 
# ggplot(analysis_data) + geom_point(aes(x = sqrt_count_2010, y = sqrt_count_2012, col = Type), size = 1, alpha = .25) +
#   facet_wrap(~State, labeller = label_bquote(State ~ .(State)))+ 
#   theme_bw() + 
#   labs(x = 'square root of 2010 houshold count', y = 'square root of 2012 houshold count') +
#   theme(text = element_text(family = 'Times'))
# ggsave(file.path(getwd(), "..", "..", "figs", 'scatter_by_state.png'), width = 6, height = 4)
# 



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
    
    full_join(tmp, y_open,  by = c("holdout sample", "Repetition")) %>% full_join(., y_type1,  by = c("holdout sample", "Repetition")) %>% 
      mutate(var_inflate = NA)
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
# anyNA(pooled_marginals)
pooled_marginals <- pooled_marginals %>% 
  mutate(n = as.factor(n), State = as.factor(State))
# var_inflate = as.factor(var_inflate)
#add artificial levels of v-inflate for the classical models



# p0 <- pooled_marginals %>% filter(is.na(var_inflate)) %>% dplyr::select(-var_inflate) %>% mutate(var_inflate = "5")
# p1 <- pooled_marginals %>% filter(is.na(var_inflate)) %>% dplyr::select(-var_inflate) %>% mutate(var_inflate = "10")
# p2 <- pooled_marginals %>% filter(is.na(var_inflate)) %>% dplyr::select(-var_inflate) %>% mutate(var_inflate = "50")
# p3 <- pooled_marginals %>% filter(is.na(var_inflate)) %>% dplyr::select(-var_inflate) %>% mutate(var_inflate = "100")


p0 <- map(as.character(v_inflate), .f = function(v){
  pooled_marginals %>% filter(is.na(var_inflate)) %>% dplyr::select(-var_inflate) %>% mutate(var_inflate = v)
})

pooled_marginals <- bind_rows(pooled_marginals %>% filter(!is.na(var_inflate)), p0) 
# Trimmed log marginal distribution 
# tlm <- function(marginals, fun = mean, trimming_fraction = 0.3){
#   log_marg <- log(marginals)
#   include <- which(log_marg >= quantile(log_marg, trimming_fraction))
#   fun(log_marg[include])
# }
#tlm(pooled_results$margNTheory[1,])

# Finding TLM with a specific base model -----


base_model <- 'Student-t' 
get_tlm_tibble <- function(marginals_tibble, base_model = base_model, trimming_fractions = c(0, .1, .2, .3)){

marginals_open_type1 <- marginals_tibble %>% 
  filter(Type1 == TRUE)

base_model_marg <-  marginals_open_type1 %>% 
  filter(Model == base_model) 

summary_tibble_by_trim <- purrr::map_dfr(trimming_fractions, .f = function(trimming_fraction){
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

  marg_split %>% 
  group_by(Repetition, Model, n, State, var_inflate) %>% 
  summarise(tlm_mean = mean(Marginal), tlm_sd = sd(Marginal)/tlm_mean) %>% 
  ungroup() %>% 
  group_by(Model, n,State, var_inflate) %>% 
  summarise(mean = mean(tlm_mean), 
            sd = sd(tlm_mean),
            mean_sd = mean(tlm_sd),
            sd_sd = sd(tlm_sd)) %>% 
  ungroup() %>% 
  mutate(Model = factor(Model, levels = c('Restricted - Huber','Rlm - Huber', 'Restricted - Tukey', 'Rlm - Tukey','Student-t', 'Normal', 'OLS')), n = factor(n), var_inflate = factor(var_inflate, levels = v_inflate), `Trimming Fraction` = factor(trimming_fraction))
})
summary_tibble_by_trim %>% 
  mutate(`Trimming Fraction` = as.factor(`Trimming Fraction`))
}

summary_tibble <- get_tlm_tibble(pooled_marginals, base_model = 'Student-t', trimming_fractions = c(0, .1, .2, .3))

theme_set(theme_bw(base_family = 'Times'))
ggplot(filter(summary_tibble, var_inflate == 100, `Trimming Fraction` == 0.3), aes(x = n, y = mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + #geom_line(position = position_dodge(width = .5)) +
  geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), width = 0.05, position  = position_dodge(width = .5), linetype = 1)  + 
  facet_wrap(~State, drop = FALSE, scales = 'free', labeller = labeller(State = label_vals)) + xlab('Percent of sample used in training set') + ylab('Average TLM')
ggsave(file.path(getwd(), "..", "..", "figs", paste0('tlm_base_',base_model, '.png')), width = 6, height = 4)


# ggplot(filter(summary_tibble, var_inflate == 100, `Trimming Fraction` == 0.3), aes(x = State, y = mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + #geom_line(position = position_dodge(width = .5)) +
#   geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), width = 0.05, position  = position_dodge(width = .5), linetype = 'dotted')  + 
#   facet_wrap(~n, drop = FALSE)
# 

# ggplot(filter(summary_tibble, !Model %in% c('OLS', 'Normal')), aes(x = n, y = mean_sd, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + #geom_line(position = position_dodge(width = .5)) +
#   geom_errorbar(mapping = aes(ymin = (mean_sd - sd_sd), ymax = (mean_sd + sd_sd)), width = 0.05, position  = position_dodge(width = .5), linetype = 'dotted') +
#   facet_wrap(~State, drop = FALSE, scales = 'free')


# ggplot(filter(summary_tibble, !Model %in% c('OLS', 'Normal'), `Trimming Fraction` == 0.3), aes(x = n, y = sd/mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .2)) + #geom_line(position = position_dodge(width = .2), lty = 2) +
#   #geom_errorbar(mapping = aes(ymin = (mean_sd - sd_sd), ymax = (mean_sd + sd_sd)), width = 0.05, position  = position_dodge(width = .5), linetype = 'dotted') +
#   facet_wrap(~State, drop = FALSE)


ggplot(filter(summary_tibble, !Model %in% c('OLS', 'Normal'), `Trimming Fraction` == 0.3), aes(x = State, y = sd/mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .2)) + geom_line(position = position_dodge(width = .2), lty = 2) +
  #geom_errorbar(mapping = aes(ymin = (mean_sd - sd_sd), ymax = (mean_sd + sd_sd)), width = 0.05, position  = position_dodge(width = .5), linetype = 'dotted') +
  facet_wrap(~n, drop = FALSE, labeller = label_both)



ggplot(filter(summary_tibble, n == 50), aes(x = as.factor(`Trimming Fraction`), y = mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .75)) + #geom_line(position = position_dodge(width = .5), lty = 2) +
  geom_errorbar(mapping = aes(ymin = (mean - sd), ymax = (mean + sd)), width = 0.05, position  = position_dodge(width = .75), linetype = 1) + facet_wrap(~State, labeller = labeller(State = label_vals), scales = 'free') + xlab(bquote(`Trimming Fraction`(alpha)))
ggsave(file.path(getwd(), "..", "..", "figs", paste0('tlm_base_',base_model, 'byTrimming.png')), width = 6, height = 4)



ggplot(analysis_data) + geom_point(aes(x = sqrt_count_2010, y = sqrt_count_2012, col = Type), size = 1, alpha = .25) +
  facet_wrap(~State, labeller = label_bquote(State ~ .(State)))+ 
  theme_bw() + 
  labs(x = 'square root of 2010 houshold count', y = 'square root of 2012 houshold count') +
  theme(text = element_text(family = 'Times'))


# acceptance rates -----

list_1 <- lapply(states, function(State_keep){  
  tbl_list <- lapply(ns, function(n){
    read_rds(paste0('single_reg_state_', State_keep, '_n_', n, '.rds' ))
    })
     })

range(sapply(list_1, function(x){
  sapply(x, function(y) range(c(y$acceptYHuber, y$acceptY)))
}))

