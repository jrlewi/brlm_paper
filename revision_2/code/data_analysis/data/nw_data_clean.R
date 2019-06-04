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
  
# View(analysisSet %>% dplyr::select(Primary_Agency_State, State))
# View(priorDataSet %>% dplyr::select(Primary_Agency_State, State))
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
  dplyr::select(Type, State, Ending_HH_Count2008, Ending_HH_Count2010)
prior_data$Type <- as_factor(as.character(prior_data$Type))

analysis_data <- analysisSet %>% 
  dplyr::select(Type, State, Ending_HH_Count2010, Ending_HH_Count2012)
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

# Center and Scale -------- 


cs <- function(x){
  cs <- (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
  cs + abs(min(cs))
  } # = (x -min(x))/sd(s)

Count_2008 <- cs(prior_data$Count_2008)
Count_2010 <- cs(c(prior_data$Count_2010, analysis_data$Count_2010))
Count_2012 <- cs(analysis_data$Count_2012)

prior_data2 <- prior_data
nr_prior <- nrow(prior_data2)
prior_data2$Count_2008 <- Count_2008
prior_data2$Count_2010 <- Count_2010[1:nr_prior]

analysis_data2 <- analysis_data
nr_analyses <- nrow(analysis_data2)
analysis_data2$Count_2010 <- Count_2010[(nr_prior + 1):(nr_prior + nr_analyses)]
analysis_data2$Count_2012 <- Count_2012



plot(prior_data$Count_2008, prior_data2$Count_2008)
plot(prior_data$Count_2010, prior_data2$Count_2010)
plot(analysis_data$Count_2010, analysis_data2$Count_2010)
plot(analysis_data$Count_2012, analysis_data2$Count_2012)


summary(fit1 <- lm(Count_2010 ~ Count_2008, data = prior_data))
summary(fit2 <- lm(Count_2010 ~ Count_2008, data = prior_data2))
plot(residuals(fit1))
plot(residuals(fit2))

summary(fit1 <- lm(Count_2012 ~ Count_2010, data = analysis_data))
summary(fit2 <- lm(Count_2012 ~ Count_2010, data = analysis_data2))
plot(residuals(fit1))
plot(residuals(fit2))



# plots ----
theme_set(theme_bw())

#by Agency Type
#+ xlim(c(0,150)) + ylim(c(0, 130))
ggplot(prior_data2, aes(x = sqrt(Count_2008), y = sqrt(Count_2010), col = Type)) + geom_point(size = .3) + stat_smooth(method = "rlm", formula = y ~ x -1, col = 1, method.args = list(psi = psi.bisquare, scale.est = 'Huber'), size = .5, lty = 2) + scale_color_discrete(guide_legend(title = 'Agency Type')) 
  
#+ xlim(c(0,150)) +  ylim(c(0, 130)) 
ggplot(analysis_data2, aes(x = sqrt(Count_2010), y = sqrt(Count_2012), col = Type)) + geom_point(size = .3)+ stat_smooth(method = "rlm", formula = y ~ x -1, col = 1, method.args = list(psi = psi.bisquare, scale.est = 'Huber'), size = .5, lty = 2) + scale_color_discrete(guide_legend(title = 'Agency Type')) 



# by state ---
ggplot(prior_data, aes(x = sqrt(Count_2008), y = sqrt(Count_2010), col = State)) + geom_point(size = .3) + xlim(c(0,150)) +  ylim(c(0, 130)) + stat_smooth(method = "rlm", formula = y ~ x -1, method.args = list(psi = psi.bisquare, scale.est = 'Huber', maxit = 100), size = .5, lty = 2, se = FALSE) 


ggplot(analysis_data, aes(x = sqrt(Count_2010), y = sqrt(Count_2012), col = State)) + geom_point(size = .3) + xlim(c(0,150)) +  ylim(c(0, 130)) + stat_smooth(method = "rlm", formula = y ~ x -1, method.args = list(psi = psi.bisquare, scale.est = 'Huber', maxit = 100), size = .5, lty = 2, se = FALSE)



# write data to rds files in approp directory ----

write_rds(prior_data2, file.path(getwd(), 'data','prior_data.rds'))

write_rds(analysis_data2, file.path(getwd(), 'data','analysis_data.rds'))


# write_rds(prior_data, file.path('/Users/john/Dropbox/school/osu/dissertationResearch/snm/journal_papers/bayes_rest_like_methods/revision_1/code/data_analysis/data/prior_data.rds'))
# 
# write_rds(analysis_data, file.path('/Users/john/Dropbox/school/osu/dissertationResearch/snm/journal_papers/bayes_rest_like_methods/revision_1/code/data_analysis/data/analysis_data.rds'))






