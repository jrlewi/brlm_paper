library(MASS)
library(tidyverse)

load('/Users/john/Dropbox/school/osu/dissertationResearch/snm/linearRegression/nwPaper1upDate/workSpaces/nwdataPrepPaper1DataCleanWorkspace.RData')

analysisSet <- as_tibble(analysisSet)
priorDataSet <-  as_tibble(priorDataSet)

analysisSet <- analysisSet %>% 
  dplyr::select(Agency_Type, Primary_Agency_State, Ending_HH_Count2010, Ending_HH_Count2012)



priorDataSet <- priorDataSet %>% 
  dplyr::select(Agency_Type, Primary_Agency_State, Ending_HH_Count2008, Ending_HH_Count2010)


#Change State Level Names -----
set.seed(124)
perm <- sample(1:(length(levels(analysisSet$Primary_Agency_State))-1))

analysisSet <- analysisSet %>% 
  mutate(State = factor(Primary_Agency_State, labels = perm ))

priorDataSet <- priorDataSet %>% 
  mutate(State = factor(Primary_Agency_State, labels = perm ))

ind_check1 <- with(analysisSet, table(Primary_Agency_State, State)) %>% 
  apply(., 2, function(x) which(x != 0))
ind_check2 <- with(priorDataSet, table(Primary_Agency_State, State))%>% 
  apply(., 2, function(x) which(x != 0))
all.equal(ind_check1, ind_check2)
all.equal(names(ind_check1), names(ind_check2))
  
# View(analysisSet %>% select(Primary_Agency_State, State))
# View(priorDataSet %>% select(Primary_Agency_State, State))
(unique(priorDataSet$State))

# Change Agency Type Level Names ----
types <- unique(c(levels(analysisSet$Agency_Type),levels(priorDataSet$Agency_Type)))

type1 <- which(types == 'CAREER')
type12 <- which(types == 'PASPORT2')
set.seed(124)
perm <- sample(seq_along(types))
perm <- sapply(perm, function(x){
  if(x == perm[which(perm == perm[type1])]){
    1
  } else {
    if(x == perm[which(perm == 1)]){
      perm[type1]
    } else {
      if(x == perm[which(perm == perm[type12])]){
        12
      } else {
        if(x == perm[which(perm == 12)]){
          perm[type12]
        } else {
      x
    }
}}}})
perm
types


change_types <- function(cur_name){
  ind <- which(types == cur_name)
  perm[ind]
}

analysisSet <- analysisSet %>% rowwise() %>% 
  mutate(Type = change_types(Agency_Type)) 
priorDataSet <- priorDataSet %>% rowwise() %>% 
  mutate(Type = change_types(Agency_Type)) 


ind_check1 <- with(analysisSet, table(Agency_Type, Type)) %>% 
  apply(., 1, function(x) colnames(.)[which(x != 0)])
ind_check2 <- with(priorDataSet, table(Agency_Type, Type)) %>% 
  apply(., 1, function(x) colnames(.)[which(x != 0)])

ind_check1
ind_check2
analysisSet %>% ungroup() %>% 
       dplyr::select(Agency_Type, Type) %>% 
  group_by(Agency_Type, Type) %>% 
  summarise(n = n())

priorDataSet %>% ungroup() %>% 
  dplyr::select(Agency_Type, Type) %>% 
  group_by(Agency_Type, Type) %>% 
  summarise(n = n())


prior_data <- priorDataSet %>% 
  select(Type, State, Ending_HH_Count2008, Ending_HH_Count2010)
prior_data$Type <- as_factor(as.character(prior_data$Type))

analysis_data <- analysisSet %>% 
  select(Type, State, Ending_HH_Count2010, Ending_HH_Count2012)
analysis_data$Type <- as_factor(as.character(analysis_data$Type))

levels(prior_data$Type)
levels(analysis_data$Type)

# reorder levels of factors -----
order_numeric_levels <- function(x) {
  factor(x, levels = levels(x)[order(as.numeric(levels(x)))])
}

prior_data <- prior_data %>% ungroup() %>% 
  mutate(Type  = order_numeric_levels(Type), State = order_numeric_levels(State))

analysis_data <- analysis_data %>% ungroup() %>% 
  mutate(Type  = order_numeric_levels(Type), State = order_numeric_levels(State))

# rename some vars ----
prior_data <- prior_data %>% dplyr::rename(Count_2008 = Ending_HH_Count2008, Count_2010 = Ending_HH_Count2010)

analysis_data <- analysis_data %>% dplyr::rename(Count_2010 = Ending_HH_Count2010, Count_2012 = Ending_HH_Count2012)

# plots ----
theme_set(theme_bw())

#by Agency Type
ggplot(prior_data, aes(x = sqrt(Count_2008), y = sqrt(Count_2010), col = Type)) + geom_point(size = .3)  + xlim(c(0,150)) + ylim(c(0, 130)) + stat_smooth(method = "rlm", formula = y ~ x -1, col = 1, method.args = list(psi = psi.bisquare, scale.est = 'Huber'), size = .5, lty = 2) + scale_color_discrete(guide_legend(title = 'Agency Type')) 
  

ggplot(analysis_data, aes(x = sqrt(Count_2010), y = sqrt(Count_2012), col = Type)) + geom_point(size = .3) + xlim(c(0,150)) +  ylim(c(0, 130)) + stat_smooth(method = "rlm", formula = y ~ x -1, col = 1, method.args = list(psi = psi.bisquare, scale.est = 'Huber'), size = .5, lty = 2) + scale_color_discrete(guide_legend(title = 'Agency Type')) 



# by state ---
ggplot(prior_data, aes(x = sqrt(Count_2008), y = sqrt(Count_2010), col = State)) + geom_point(size = .3) + xlim(c(0,150)) +  ylim(c(0, 130)) + stat_smooth(method = "rlm", formula = y ~ x -1, method.args = list(psi = psi.bisquare, scale.est = 'Huber'), size = .5, lty = 2, se = FALSE) 


ggplot(analysis_data, aes(x = sqrt(Count_2010), y = sqrt(Count_2012), col = State)) + geom_point(size = .3) + xlim(c(0,150)) +  ylim(c(0, 130)) + stat_smooth(method = "rlm", formula = y ~ x -1, method.args = list(psi = psi.bisquare, scale.est = 'Huber'), size = .5, lty = 2, se = FALSE)



# write data to rds files in approp directory ----

write_rds(prior_data, file.path('/Users/john/Dropbox/school/osu/dissertationResearch/snm/journal_papers/bayes_rest_like_methods/revision_1/code/data_analysis/data/prior_data.rds'))

write_rds(analysis_data, file.path('/Users/john/Dropbox/school/osu/dissertationResearch/snm/journal_papers/bayes_rest_like_methods/revision_1/code/data_analysis/data/analysis_data.rds'))






