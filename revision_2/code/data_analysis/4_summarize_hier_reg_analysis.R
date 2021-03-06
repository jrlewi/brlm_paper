# Summarize the output from pooled regression analysis.
library(tidyverse)
library(MASS)
library(abind)

ns <- c(1547)
#process in batches? combine results first
num_batch <- 5
rds_path <- "~/Dropbox/school/osu/dissertationResearch/snm/journal_papers/brlm_paper/revision_2/code/data_analysis/hier_reg_n1547_sim_number_"
tmp1 <- readRDS(file.path(paste0(rds_path, 1, '.rds')))

#combine results into one
for(i in 2:num_batch){
  tmp <-  readRDS(file.path(paste0(rds_path, i, '.rds')))
  
  tmp1$y_hold <- rbind(tmp1$y_hold, tmp$y_hold)
  tmp1$y_open <- rbind(tmp1$y_open, tmp$y_open)
  tmp1$y_type1 <- rbind(tmp1$y_hold, tmp$y_type1)
  tmp1$holdIndices <- rbind(tmp1$holdIndices, tmp$holdIndices)
  tmp1$group_estimates <- abind(tmp1$group_estimates, tmp$group_estimates, along = 4)
  tmp1$group_estimates_sds <- abind(tmp1$group_estimates_sds, tmp$group_estimates_sds, along = 4)
  tmp1$estimates <- abind(tmp1$estimates, tmp$estimates, along = 3)
  tmp1$estimates_sd <- abind(tmp1$estimates_sd, tmp$estimates_sd, along = 3)
  tmp1$marginals <- abind(tmp1$marginals, tmp$marginals, along = 2)
  tmp1$marginals_sd <- abind(tmp1$marginals_sd, tmp$marginals_sd, along = 2)
  tmp1$predictions <- abind(tmp1$predictions, tmp$predictions, along = 2)
  tmp1$acceptY <- abind(tmp1$acceptY, tmp$acceptY, along = 3)
  tmp1$group_converge <- abind(tmp1$group_converge, tmp$group_converge, along = 4) 
  tmp1$converge <- abind(tmp1$converge, tmp$converge, along = 3) 
  
  }

#write combined results to rds.
write_rds(tmp1, file.path(file.path(getwd(), paste0('hier_reg_n', ns,".rds"))))


analysis_data <- read_rds(file.path(here::here(), 'data', 'analysis_data.rds'))
analysis_data <- analysis_data %>% 
  mutate(sqrt_count_2010 = sqrt(Count_2010), sqrt_count_2012 = sqrt(Count_2012)) %>%
  group_by(State) %>% 
  filter(n() >= 25) %>% ungroup() %>% 
  mutate(State = factor(State)) %>% #filter(!State %in% c(14, 30)) %>% 
  arrange(State)

analysis_data <- analysis_data %>% 
  filter(!State %in% c('14', '30')) %>% 
  mutate(State = droplevels(State))


# tmp <- readRDS("~/Dropbox/school/osu/dissertationResearch/snm/journal_papers/brlm_paper/revision_1/code/data_analysis/hier_reg_n1590.rds")
# holdIndices <- tmp$holdIndices

#creates a dataframe with holdout indices and states for each rep
holdout_tibble <- function(holdIndices){
  holdIndices <- t(holdIndices)
  rownames(holdIndices) <- 1:nrow(holdIndices)
  holdIndices <- as_tibble(holdIndices, rownames = 'holdout sample')
  
  holdIndices <- holdIndices %>% 
    rename_at(vars(contains('V')), funs(sub('V', '', .))) %>% 
    gather('Repetition','holdout index', -"holdout sample") %>% 
    mutate(Repetition = as.numeric(Repetition), `holdout sample` = as.numeric(`holdout sample`))
  
  #grab the state
  holdIndices %>%
    split(list(.$Repetition)) %>% 
    map(.f = function(x){
      holddata <- analysis_data[x$`holdout index`,]
      x %>% mutate(State = holddata$State, Type =holddata$Type, Open = holddata$Count_2012 > 0)
    }) %>% 
    bind_rows() 
}


get_marginals <- function(hier_results, n){
  
hold_index <- holdout_tibble(hier_results$holdIndices)  

margs <- apply(hier_results$marginals, 1, function(marg) {
  marg <- t(marg)
  rownames(marg) <- 1:nrow(marg)
  marg <- as_tibble(marg, rownames = 'holdout sample')
  marg <- marg %>% 
    rename_at(vars(contains('V')), funs(sub('V', '', .))) %>% 
    gather('Repetition','Marginal', -"holdout sample") %>% 
    mutate(Repetition = as.numeric(Repetition), `holdout sample` = as.numeric(`holdout sample`))
  
full_join(marg, hold_index, by = c("holdout sample", "Repetition")) 
})
margs_sd <- apply(hier_results$marginals_sd, 1, function(marg) {
  marg <- t(marg)
  rownames(marg) <- 1:nrow(marg)
  marg <- as_tibble(marg, rownames = 'holdout sample')
  marg <- marg %>% 
    rename_at(vars(contains('V')), funs(sub('V', '', .))) %>% 
    gather('Repetition','SD', -"holdout sample") %>% 
    mutate(Repetition = as.numeric(Repetition), `holdout sample` = as.numeric(`holdout sample`))
  
full_join(marg, hold_index, by = c("holdout sample", "Repetition")) 
})


names(margs) <-names(margs_sd) <- c("OLS", 
                                    "Rlm - Tukey",
                                    "Rlm - Huber",
                                    "Normal", 
                                    "Restricted - Tukey", 
                                    "Restricted - Huber",
                                    "Student-t")
margs <- bind_rows(margs, .id = 'Model') %>% 
  mutate(n = factor(n, levels = ns))

margs_sd <- bind_rows(margs_sd, .id = 'Model') %>% 
  mutate(n = factor(n, levels = ns))
full_join(margs, margs_sd,  by = c("Model", "holdout sample", "Repetition", "holdout index", "State", "Type", "Open", "n"))
}


# combine hier results in batch to one hier results



hier_marginals <- ns %>% map(.f = function(n){
  hier_results <- readRDS(file.path(getwd(), paste0('hier_reg_n', n,".rds")))
  get_marginals(hier_results,n)
}) %>% 
  bind_rows()


# Finding TLM with a specific base model, by state -----
#base_model <- 'Student-t'
get_tlm_tibble <- function(marginals_tibble, base_model = base_model , trimming_fractions = c(0, .1, .15, 0.2,.25, .3)){

marginals_open_type1 <- marginals_tibble %>% 
  filter(Type == '1')

base_model_marg <-  marginals_open_type1 %>% 
  filter(Model == base_model) 

summary_tibble_by_trim <- purrr::map_dfr(trimming_fractions, .f = function(trimming_fraction){
lower_tlm_vals <- base_model_marg %>% 
  group_by(Repetition, State, n) %>% 
  summarise(lower_tlm = quantile(Marginal, trimming_fraction)) 
#for each repetition, find the holdout sample indexes to include from the base model
holdout_include <- left_join(base_model_marg,lower_tlm_vals, by = c('Repetition', 'State', 'n')) %>% 
  mutate(include = Marginal > lower_tlm) %>% 
  dplyr::select(-Model, -Marginal, -SD)
# %>% 
#   filter(include == TRUE)

holdout_include_all <- c('Restricted - Huber','Rlm - Huber', 'Restricted - Tukey', 'Rlm - Tukey','Student-t', 'Normal', 'OLS') %>% map_dfr(.f = function(model){
    tmp <- marginals_open_type1 %>% 
      filter(Model == model) 
    full_join(holdout_include, tmp, by = c('holdout sample', 'holdout index', 'State', 'Repetition', 'Type', 'Open', 'n')
) }) %>% filter(include == TRUE)



  holdout_include_all %>%  
  group_by(Repetition, Model, State, n) %>% 
  mutate(Marginal = log(Marginal + 1e-300)) %>% #log, avoid -Inf
  summarise(tlm = mean(Marginal), tlm_sd = sd(Marginal)) %>% 
  ungroup() %>% 
  # group_by(Model, State, n) %>% 
  # summarise(mean = mean(tlm_mean), 
  #           sd = sd(tlm_mean)) %>% 
  # ungroup() %>% 
  mutate(Model = factor(Model, levels = c('Restricted - Huber','Rlm - Huber', 'Restricted - Tukey', 'Rlm - Tukey','Student-t', 'Normal', 'OLS')), n = factor(n), `Trimming Fraction` = trimming_fraction)
})
summary_tibble_by_trim %>% 
  mutate(`Trimming Fraction` = as.factor(`Trimming Fraction`))
}

summary_tibble <- get_tlm_tibble(hier_marginals, base_model = 'Student-t', trimming_fractions = c(0, .1, .15, 0.2,.25, .3))


# overall average ------
#by rep - average over states, then average over reps
overall_average <- summary_tibble %>% 
  group_by(Model, Repetition, n, `Trimming Fraction`) %>%
  summarise(mean_tlm = mean(tlm)) %>% #ave over states
  ungroup() %>% 
  group_by(Model, n,  `Trimming Fraction`) %>% #ave/sd over reps
  summarize(mean = mean(mean_tlm), sd = sd(mean_tlm))
  
#K <- length(unique(summary_tibble$Repetition))

theme_set(theme_bw(base_family = 'Times'))
ggplot(overall_average %>% filter(`Trimming Fraction` %in% c(0.1,.15,.2,.25,.3)) , aes(x = `Trimming Fraction` , y = mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), width = 0.05, position  = position_dodge(width = .5), linetype = 1) + labs(y = 'Average TLM') + xlab(bquote(`Trimming Fraction`(alpha))) 
ggsave(file.path(getwd(), "..", "..", "figs", 'hier_average_tlm.png'), width = 6, height = 4)


ggplot(overall_average %>%  filter(!Model %in% c('OLS', 'Normal', 'Student-t')), aes(x = `Trimming Fraction` , y = sd, col = Model, group = Model)) + geom_point(position = position_dodge(width = .0))  + geom_line(position  = position_dodge(width = .0), linetype = 1) + xlab(bquote(`Trimming Fraction`(alpha)))
ggsave(file.path(getwd(), "..", "..", "figs", 'hier_sd_tlm.png'), width = 6, height = 4)



#Averages by State

average_by_state <- summary_tibble %>% 
  group_by(Model, State, n, `Trimming Fraction`) %>%
  summarise(mean = mean(tlm), sd = sd(tlm)) %>% 
  ungroup() 

state_count <- analysis_data %>% group_by(State) %>% 
  summarize(state_count = n())

average_by_state <- left_join(average_by_state, state_count, by = 'State') %>% unite(col = 'State(n)', c('State', 'state_count'), sep = '(', remove = FALSE) %>% 
  mutate(')' = ')') %>% unite(col = 'State(n)', c('State(n)', ')'), sep = '', remove = FALSE) %>%  mutate(`State(n)` = reorder(as.factor(`State(n)`), state_count))


theme_set(theme_bw(base_family = 'Times'))
ggplot(average_by_state  %>% filter(`Trimming Fraction` == '0.3', Model %in% c('Restricted - Tukey', 'Rlm - Tukey'), State != '28'), aes(x = `State(n)` , y = mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), width = 0.05, position  = position_dodge(width = .5), linetype = 1) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "Average TLM") #+ scale_y_log10() + 
ggsave(file.path(getwd(), "..", "..", "figs", 'hier_ave_tlm_state.png'), width = 7.5, height = 5)


theme_set(theme_bw(base_family = 'Times'))
ggplot(average_by_state  %>% filter(`Trimming Fraction` == '0.3', Model %in% c('Restricted - Huber', 'Rlm - Huber'), State != '28'), aes(x = `State(n)` , y = mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), width = 0.05, position  = position_dodge(width = .5), linetype = 1) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "Average TLM") #+ scale_y_log10() + 
ggsave(file.path(getwd(), "..", "..", "figs", 'hier_ave_tlm_state_huber.png'), width = 7.5, height = 5)




# theme_set(theme_bw(base_family = 'Times'))
# ggplot(average_by_state  %>% filter(`Trimming Fraction` == '0.3', State %in% c(2,15,27,36)), aes(x = State , y = mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), width = 0.05, position  = position_dodge(width = .5), linetype = 1) 
# 
# ggplot(average_by_state  %>% filter(State %in% c(2,15,27,36)), aes(x = `Trimming Fraction` , y = mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), width = 0.05, position  = position_dodge(width = .5), linetype = 1) + facet_wrap(~State, scales = 'fixed')
# 
# 
# ggplot(average_by_state %>% filter(!Model %in% c('OLS', 'Normal')), aes(x = `Trimming Fraction` , y = mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), width = 0.05, position  = position_dodge(width = .5), linetype = 1) + facet_wrap(~State, scales = 'free')
# 
# 




# 
# label_vals <- c('2' = 'State 2', '15' = 'State 15', '27' = 'State 27', '36' = 'State 36')
# 
# ggplot( average_by_state  %>% filter(State %in% c(2,15,27,36), !Model %in% c('OLS', 'Normal', 'Student-t')), aes(x = `Trimming Fraction` , y = sd/mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = 0)) + geom_line(position = position_dodge(width = 0)) + facet_wrap(~State, scales = 'free', labeller = labeller(State = label_vals))
# 
# 
# # ggplot( average_by_state %>% filter(!Model %in% c('OLS', 'Normal', 'Student-t'), `Trimming Fraction` == '0.3'), aes(x =state_count , y = sd/mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = 0)) + geom_line(position = position_dodge(width = 0)) + scale_x_log10() 
# 
# 
# filter_df <- average_by_state  %>% filter(!Model %in% c('OLS', 'Normal', 'Student-t'))
# filter_df <- filter_df %>% mutate(State = as.factor(paste0('State ', State)))
# ggplot(filter_df, aes(x = `Trimming Fraction` , y = sd/mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = 0)) + geom_line(position = position_dodge(width = 0)) + facet_wrap(~State, scales = 'free')

theme_set(theme_bw(base_family = 'Times'))
sub_df <- average_by_state  %>% filter(`Trimming Fraction` == '0.3', Model %in% c('Restricted - Tukey', 'Rlm - Tukey'))
ggplot(sub_df, aes(x = `State(n)` , y = sd/abs(mean), col = Model, group = Model)) + geom_point(position = position_dodge(width = 0)) + geom_line() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+ scale_y_log10()
ggsave(file.path(getwd(), "..", "..", "figs", 'hier_sd_tlm_state.png'), width = 7.5, height = 5)

# ----- 




#relatations to state sizes ----
# state_count <- analysis_data %>% 
#   group_by(State) %>% 
#   summarise(state_count = n())
# n <- ns
# hier_results <- readRDS(file.path(getwd(), paste0('hier_reg_n', n,".rds")))
# names(hier_results )
# ols_ind <- 1; 
# rlm_ind <- 2; 
# rlm_huber_ind <- 3
# norm_ind <- 4
# rest_ind <- 5
# rest_hub_ind <- 6
# t_ind <- 7
# 
# 
# 
# state_ind <- which(state_count$State == 36)
# 
# dim(hier_results$group_estimates)
# plot(hier_results$group_estimates[rest_ind, 1, state_ind, ], hier_results$group_estimates[rlm_ind, 1, state_ind, ])
# abline(0,1)
# plot(hier_results$group_estimates[rest_ind, 2, state_ind, ], hier_results$group_estimates[rlm_ind, 2, state_ind, ])
# abline(0,1)
# # abline(h = b0/(a0-1), v = b0/(a0-1))
# 
# 
# 
# 
# state_count <- left_join(average_by_state,state_count, by = "State") 
# sub_df <- state_count  %>% filter(`Trimming Fraction.x` == '0.3', Model.x %in% c('Restricted - Tukey', 'Rlm - Tukey'))
# 
# ggplot(sub_df, aes(x = as.factor(state_count) , y = sd.x, col = Model.x, group = Model.x)) + geom_point(position = position_dodge(width = 0)) + geom_line() 
# 
# ggplot(sub_df, aes(x = state_count , y = mean, col = Model, group = Model))  + geom_point(position = position_dodge(width =  0))  + geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), width = 0.0, position  = position_dodge(width = 0), linetype = 1) + scale_x_log10()



#Convergence Diagnostics ----
# tmp <- readRDS("~/Dropbox/school/osu/dissertationResearch/snm/journal_papers/brlm_paper/revision_1/code/data_analysis/hier_reg_n1590.rds")
range(tmp1$acceptY)
summary(tmp1$acceptY)
sum(tmp1$acceptY<.1)
# names(tmp)
# lapply(tmp, dim)
# 
# 
# plot(tmp$marginals[3,,],tmp$marginals[6,,], cex = .2)
# abline(0,1, col = 2)
# 
# plot(as.numeric(tmp$marginals[3,1,])-as.numeric(tmp$marginals[6,1,]), cex = .2)
# 
# 
# plot(tmp$group_estimates[2,1,,],tmp$group_estimates[5,1,,])
# abline(0,1)
# abline(h = tmp$estimates[5,1,1])
# 
# 
# plot(tmp$converge, pch = 19, cex = .2)
# mean(abs(tmp$converge[7,,]) > 2, na.rm = TRUE)
# 
# plot(tmp$group_converge, pch = 19, cex = .2)
# mean(abs(tmp$group_converge[-7,1,,]) > 2, na.rm = TRUE)
# mean(abs(tmp$group_converge[-7,,,]) > 2, na.rm = TRUE)
# 
# mean(abs(tmp$group_converge[-7,,,]) > 2, na.rm = TRUE)
# 
# which(tmp$group_converge > 5)
# dim(tmp$group_converge)
# which(tmp$group_converge[5,1,,] > 5)/22
# which(tmp$group_converge[5,2,,] > 5)/22
# 
# which(tmp$group_converge[7,1,,] > 5)/ 22
# which(tmp$group_converge[7,2,,] > 5)/22
# 
# dim(tmp$group_converge)
# 
# 

# Evaluate predictions for specific states ----------
# states <- c(2, 15, 27, 36) # state we care about
# 
# hier_marginals_filter <- ns %>% map(.f = function(n){
#   hier_results <- readRDS(file.path(getwd(), paste0('hier_reg_n', n,".rds")))
#   get_marginals(hier_results,n)
# }) %>% 
#   bind_rows() %>% 
#   filter(State %in% states)
# 

# get_tlm_tibble_by_state <- function(marginals_tibble, base_model = base_model, trimming_fractions = c(0, .1, .2, .3), state_filter = c(2, 15, 27, 36)){
#   marginals_tibble <- marginals_tibble %>% filter(State %in% state_filter) %>%   mutate(State = factor(State))
#   
#   marginals_open_type1 <- marginals_tibble %>% 
#     filter(Type == '1')
#   
#   base_model_marg <-  marginals_open_type1 %>% 
#     filter(Model == base_model) 
#   
#   summary_tibble_by_trim <- purrr::map_dfr(trimming_fractions, .f = function(trimming_fraction){
#     lower_tlm_vals <- base_model_marg %>% 
#       group_by(Repetition, n, State) %>% 
#       summarise(lower_tlm = quantile(Marginal, trimming_fraction)) 
#     
#     
#     #for each repetition, find the holdout sample indexes to include from the base model----
#     holdout_include <- left_join(base_model_marg,lower_tlm_vals, by = c('Repetition', 'n', "State")) %>% 
#       mutate(include = Marginal > lower_tlm) %>% 
#       dplyr::select(`holdout sample`, Repetition, include, n, State, lower_tlm) %>% 
#       filter(include == TRUE)
#     
#     include_by_rep <- holdout_include %>%
#       split(list(.$Repetition, .$n, .$State)) %>% map(function(x) x$`holdout sample`)
#     
#     
#     #include_by_rep: has the holdout sample included for tlm calculation for each rep. 
#     #for each Repetition value, marginals_open_type_1 will have the sample houldout samples for each model. 
#     marg_split <- marginals_open_type1 %>% 
#       split(list(.$Repetition, .$n, .$State)) %>% map2(.x = ., .y = include_by_rep, .f = function(x, y) filter(x,  x$`holdout sample` %in% y)) %>% 
#       bind_rows()
#     
#     marg_split %>% 
#       group_by(Repetition, Model, n, State) %>% 
#       summarise(tlm_mean = mean(Marginal), tlm_sd = sd(Marginal)/tlm_mean) %>%
#       ungroup() %>% 
#       group_by(Model, n,State) %>% 
#       summarise(mean = mean(tlm_mean), 
#                 sd = sd(tlm_mean),
#                 mean_sd = mean(tlm_sd),
#                 sd_sd = sd(tlm_sd)) %>% 
#       ungroup() %>% 
#       mutate(Model = factor(Model, levels = c('Restricted - Huber','Rlm - Huber', 'Restricted - Tukey', 'Rlm - Tukey','Student-t', 'Normal', 'OLS')), n = factor(n), `Trimming Fraction` = factor(trimming_fraction))
#   })
#   summary_tibble_by_trim %>% 
#     mutate(`Trimming Fraction` = as.factor(`Trimming Fraction`))
# }
# 
# summary_tibble_by_state <- get_tlm_tibble_by_state(hier_marginals_filter , base_model = 'Student-t', trimming_fractions = c(0, .1, .2, .3), state_filter = c(2, 15, 27, 36))
# 
# theme_set(theme_bw(base_family = 'Times'))
# ggplot(filter(summary_tibble_by_state, `Trimming Fraction` == 0.3), aes(x = State, y = mean, col = Model, group = Model)) + geom_point(position = position_dodge(width = .5))  + #geom_line(position = position_dodge(width = .5)) +
#   geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), width = 0.05, position  = position_dodge(width = .5), linetype = 1)  + 
#   facet_wrap(~State, drop = FALSE, scales = 'free', labeller = labeller(State = label_vals)) + xlab('Percent of sample used in training set') + ylab('Average TLM')
# 


