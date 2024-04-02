####This script runs all final survival models####
####First the script fits univariate models and then extracts BLUP values used in survival analysis####

##Packages required to run script
require(parallel)
require(remotes)
require(ggplot2)
require(knitr)
require(readr)
require(ggeffects)
require(dplyr)
require(tidybayes)
require(coda)
require(lubridate)
require(posterior)
require(brms)
require(MCMCglmm)

## theme for plots
theme_turkey<- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size=12),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(size = 15),
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          #legend.position = c(0.9, 0.9),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.background = element_rect(color = "black",
                                           fill = "transparent",
                                           size = 4, linetype = "blank"))
}


# run this is a fresh R session or restarting your current session
library(rstan)


# Bayes stuff/options to make model run more quickly
options(mc.cores = 32,  # Use 4 cores
        brms.backend = "rstan")
options(mc.cores = parallel::detectCores(),
        brms.backend = "rstan")
rstan_options(auto_write=T)
bayes_seed <- 1234


# import data and add tag ID
library(readr)
d1<- read_csv("Data/Male_gps_Webb.csv")#load data
d2<- read_csv("Data/Male_gps_GA_sites.csv") #load data

require(dplyr)
data<-full_join(d1,d2)


# filter out data points in which turkeys died that particular day
# there is wonky values on these days since hunters/predators moved carcasses around
# or the bird died on the roost and didn't move at all 
# so better to filter these out so they don't influence the results
data<-data%>%filter(Status=="Live")
data<-data%>%filter(Age=="A")
# filter out values where the total distance traveled in a given day is less than 100m 
# gps error rate is typically ~30m 
# to account for error in every cardinal direction we picked 100m as a reasonable distance
data<-data%>%filter(total_dist>100)

# manipulate stage so the intercept is 'Pre-hunt' 
data<-data%>% group_by(ID)%>%
  mutate(Stage= replace(Stage, Stage== "Pre_Hunt", "1Pre_Hunt"))%>%
  mutate(Stage= replace(Stage, Stage == "Hunt", "2Hunt"))%>%
  mutate(Stage= replace(Stage, Stage == "Post_Hunt", "3Post_Hunt"))
#Change site names
data<-data%>% group_by(ID)%>%
  mutate(Site= replace(Site, Site== "BFG", "GA"))%>%
  mutate(Site= replace(Site, Site == "CC", "GA"))
# make ID a factor
data$ID<-as.factor(data$ID)
data$Site<-as.factor(data$Site)
data$Stage<-as.factor(data$Stage)


# standardize trait so model converges more quickly 
# and are easier to compare distances to speed for subsequent analysis
data$Dist_access_c<-scale(data$Dist_Access, scale = T, center = T)#risk-taking
data$mean_speed_c<-scale(data$mean_speed, scale = T, center = T)#exploration
data$Open_c<-scale(data$Open, scale = T, center = T)#risk-taking
data$Edge_c<-scale(data$dist_edge, scale = T, center = T)#risk-taking



##Univariate model used in survival analysis----
require(brms)
data$Dist_access_bs<-data$Dist_access_c*-1#multiple by -1 to ease interpretation
                                                                  #now a postive log-likelihood = more risky individuals have greater survival
#1. Public Access-Risk-taking univariate model to extract BLUPs----
#run both populations together to obtain unbiased BLUP values
model.traits1<-brm(Dist_access_bs~ Stage*Site + Site + (1|ID),
                   data= data, family="gaussian",
                   warmup=300, iter=15300,
                   chains=2,thin= 15, init="random",
                   cores = 32,
                   seed=bayes_seed,
                   backend = "rstan",
                   control = list(adapt_delta = 0.99, max_treedepth=15)) 
#model summary
summary(model.traits1) 
#check model convergence via trace/denisty plots
# plot(model.traits1)#suppress so script will run when hitting source
require(bayestestR)
#compute median, 95% CrI, and Rhat values
describe_posterior(model.traits1)

# variables(model.traits1)

#variance for ID
var.ID<- posterior_samples(model.traits1)$"sd_ID__Intercept"^2
mean(var.ID);HPDinterval(as.mcmc(var.ID),0.95)#this computes CrI 
#residual variance
var.res<-posterior_samples(model.traits1)$"sigma"^2
mean(var.res);HPDinterval(as.mcmc(var.res),0.95)#this computes CrI

#repeatability (ID variance/ (ID variance + residual variance))
repeatability<-var.ID/(var.ID + var.res)
describe_posterior(repeatability)

#extact posterior samples
o1 <- as.data.frame(model.traits1) #extracts all posterior samples
o1<-o1%>%select(9:117) #extract posterior samples of BLUP values for each individual
o1<-as.mcmc(o1) #make mcmc object

#load survival dataframe
survival_final<-read_csv("Data/survival_dataframe.csv") 
survival_final<-survival_final%>%group_by(Age)%>% filter(Age=="A") #subset to only include adults
require(bayestestR)

#Georgia----
# 1a. Access + survival (harvest)----
# Filter survival data frame for birds in Georgia
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="GA")
# Filter for individuals that survived or were harvested
harv_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Harvested")

#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, harv_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

#extract coefficients from 2000 models 
harv_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)


#### GA-Log-odd ratio for the effect of human harvest on risk-taking (distance to hunter access)----
# plot(harv_ga_da) #plot posterior distribution of log-odd ratios
                 #suppress plot so script will run when hitting source
harv_ga_da_sur<-data.frame(harv_ga_da) #make coefficients data frame so they can viewed using function below
describe_posterior(harv_ga_da_sur)



# 1b. Access + survival (predation)----

# Filter for Georgia
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="GA")
# Filter for individuals that survived or were predated
pred_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Predated")

#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, pred_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

pred_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)


#### GA-Log-odd ratio for the effect of predation on risk-taking (distance to hunter access)----
# plot(pred_ga_da) #plot posterior distribution of log-odd ratios
                 #suppress plot so script will run when hitting source
pred_ga_da_sur<-data.frame(pred_ga_da) #make posterior a data frame to compute median and 95% CrI
describe_posterior(pred_ga_da_sur)



#South Carolina----
# 1aa. Access + survival (harvest)----

# Filter survival data frame for birds in South Carolina
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="Webb")
# Filter for individuals that survived or were harvested
harv_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Harvested")

#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, harv_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

#extract coefficients from 2000 models 
harv_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)



#### SC-Log-odd ratio for the effect of human harvest on risk-taking (distance to hunter access)----
# plot(harv_ga_da) #plot posterior distribution of log-odd ratios
                 #suppress plot so script will run when hitting source
harv_ga_da_sur<-data.frame(harv_ga_da) #make coefficients data frame so they can viewed using function below
describe_posterior(harv_ga_da_sur)



# 1bb. Access + survival (predation)----

# Filter for South Carolina
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="Webb")
# Filter for individuals that survived or were predated
pred_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Predated")

#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, pred_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

pred_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)


#### SC-Log-odd ratio for the effect of predation on risk-taking (distance to hunter access)----
# plot(pred_ga_da) #plot posterior distribution of log-odd ratios
                 #suppress plot so script will run when hitting source
pred_ga_da_sur<-data.frame(pred_ga_da) #make posterior a data frame to compute median and 95% CrI
describe_posterior(pred_ga_da_sur)




#2. Exploration (speed) univariate model to extract BLUPs----
#run both populations together to obtain unbiased BLUP values
model.traits2<-brm(mean_speed_c~ Stage*Site + Site +  (1|ID) ,
                   data= data, family="gaussian",
                   warmup=300, iter=10300,
                   chains=2,thin= 10, init="random",
                   cores = 12,
                   seed=bayes_seed,
                   backend = "rstan",
                   control = list(adapt_delta = 0.99, max_treedepth=15)) # within chain parallelization, number of threads + number of cores should match total cores you intend to use) 

#model summary
summary(model.traits2)
# #plot markov chains and check model convergence
# plot(model.traits2)#suppress plot so script will run when hitting source
require(bayestestR)
#compute median, 95% CrI and Rhat values
describe_posterior(model.trait2)
variables(model.traits2)

#variance for ID
var.ID<- posterior_samples(model.traits2)$"sd_ID__Intercept"^2
mean(var.ID);HPDinterval(as.mcmc(var.ID),0.95)#this computes CrI 
#residual variance
var.res<-posterior_samples(model.traits2)$"sigma"^2
mean(var.res);HPDinterval(as.mcmc(var.res),0.95)#this computes CrI

#repeatability estimate
repeatability<-var.ID/(var.ID + var.res)
describe_posterior(repeatability)
library(emmeans)


#extract posterior samples
o1 <- as.data.frame(model.traits2)
o1<-o1%>%select(9:117)
o1<-as.mcmc(o1)
#load survival dataframe
survival_final<-read_csv("Data/survival_dataframe.csv")
survival_final<-survival_final%>%group_by(Age)%>% filter(Age=="A")
require(bayestestR)

#Georgia----
# 2a. Exploration + survival (harvest)----
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="GA")
# Filter for individuals that survived or were harvested
harv_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Harvested")

#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, harv_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

harv_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)


#### GA-Log-odd ratio for the effect of harvest on exploration----
# plot(harv_ga_da)#plot posterior distribution of log-odd ratios of survival
               #suppress plot so script will run when hitting source 
harv_ga_da_sur<-data.frame(harv_ga_da) #make data frame to compute survival estimates
describe_posterior(harv_ga_da_sur)


# 2b. Exploration + speed (predation)----
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="GA")
# Filter for individuals that survived or were predated
pred_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Predated")


#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, pred_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

pred_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)



#### GA-Log-odd ratio for the effect of predation on exploration (speed)----
# plot(pred_ga_da) #posterior distribution of log odd ratios for survival
                   #suppress plot so script will run when hitting source
pred_ga_da_sur<-data.frame(pred_ga_da) #make data frame to compute median, 95% CrI, and Rhat values
describe_posterior(pred_ga_da_sur)



#South Carolina----
# 2aa. Exploration + survival----
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="Webb")
# Filter for individuals that survived or were harvested
harv_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Harvested")


#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, harv_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

harv_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)



#### SC-Log-odd ratio for the effect of harvest on exploration----
# plot(harv_ga_da)#plot posterior distribution of log-odd ratios of survival
                #suppress plot so script will run when hitting source
harv_ga_da_sur<-data.frame(harv_ga_da) #make data frame to compute survival estimates
describe_posterior(harv_ga_da_sur)


# 2bb. Speed + survival (predation))----
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="Webb")
# Filter for individuals that survived or were predated
pred_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Predated")


#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, pred_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

pred_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)



#### Log-odd ratio for the effect of predation on exploration (speed)----
# plot(pred_ga_da) #posterior distribution of log odd ratios for survival
                 #suppress plot so script will run when hitting source
pred_ga_da_sur<-data.frame(pred_ga_da) #make data frame to compute median, 95% CrI, and Rhat values
describe_posterior(pred_ga_da_sur)



#3. Open landcover-Risk-taking univariate model to extract BLUPs----
#run both populations together to obtain unbiased BLUP values
data$Open_bs<-data$Open_c*-1#mutltiply bt -1 to ease interpretation
                            #postive log-likelihood = greater survival for riskier individuals
model.traits3<-brm(Open_bs~ Stage*Site + Site +  (1|ID) ,
                   data= data, family="gaussian",
                   warmup=300, iter=10300,
                   chains=2,thin= 10, init="random",
                   cores = 12,
                   seed=bayes_seed,
                   backend = "rstan",
                   control = list(adapt_delta = 0.99, max_treedepth=15)) # within chain parallelization, number of threads + number of cores should match total cores you intend to use) 

#summary of model
summary(model.traits3)
# plot(model.traits)#plot markov chains and density plots
                  #suppress plot so script will run when hitting source
require(bayestestR)
describe_posterior(model.traits3) #compute median, 95% CrI, and Rhat values
variables(model.traits3)

#variance for ID
var.ID<- posterior_samples(model.traits3)$"sd_ID__Intercept"^2
mean(var.ID);HPDinterval(as.mcmc(var.ID),0.95)#this computes CrI 
#residual variance
var.res<-posterior_samples(model.traits3)$"sigma"^2
mean(var.res);HPDinterval(as.mcmc(var.res),0.95)#this computes CrI

#repeatability estimate
repeatability<-var.ID/(var.ID + var.res)
describe_posterior(repeatability)


#extract posterior samples
o1 <- as.data.frame(model.traits3)
o1<-o1%>%select(9:117)
o1<-as.mcmc(o1)
#load survival data frame
survival_final<-read_csv("Data/survival_dataframe.csv")
survival_final<-survival_final%>%group_by(Age)%>% filter(Age=="A")
require(bayestestR)

#Georgia----
# 3a. open + Harvest----
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="GA")
# Filter for individuals that survived or were harvested
harv_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Harvested")


#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, harv_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}


harv_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)


#### GA-Log-odd ratio for the effect of harvest on risk-taking (distance to open)----
# plot(harv_ga_da)#suppress plot so script will run when hitting source
harv_ga_da_sur<-data.frame(harv_ga_da) 
describe_posterior(harv_ga_da_sur)


# 3b. open + predation----
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="GA")
# Filter for individuals that survived or were predated
pred_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Predated")


#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, pred_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

pred_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)


#### GA-Log-odd ratio for the effect of predation on risk-taking (distance to open)----
# plot(pred_ga_da)#suppress plot so script will run when hitting source
pred_ga_da_sur<-data.frame(pred_ga_da) 
describe_posterior(pred_ga_da_sur)




#South Carolina----
# 3aa. open + Harvest----
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="Webb")
# Filter for individuals that survived or were harvested
harv_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Harvested")

#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, harv_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}


harv_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)


#### SC-Log-odd ratio for the effect of harvest on risk-taking (distance to open)----
# plot(harv_ga_da)#suppress plot so script will run when hitting source
harv_ga_da_sur<-data.frame(harv_ga_da) 
describe_posterior(harv_ga_da_sur)


# 3bb. open + predation----
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="Webb")
# Filter for individuals that survived or were predated
pred_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Predated")


#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, pred_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

pred_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)

#### SC-Log-odd ratio for the effect of predation on risk-taking (distance to open)----
# plot(pred_ga_da)#suppress plot so script will run when hitting source
pred_ga_da_sur<-data.frame(pred_ga_da) 
describe_posterior(pred_ga_da_sur)




#4. Edge landcover- risk-taking univariate model to extract BLUPs----
#run both populations together to obtain unbiased BLUP values
model.traits4<-brm(Edge_c~ Stage*Site + Site +  (1|ID) ,
                   data= data, family="gaussian",
                   warmup=300, iter=10300,
                   chains=2,thin= 10, init="random",
                   cores = 12,
                   seed=bayes_seed,
                   backend = "rstan",
                   control = list(adapt_delta = 0.99, max_treedepth=15)) # within chain parallelization, number of threads + number of cores should match total cores you intend to use) 
summary(model.traits4)
# plot(model.traits4)#suppress plot so script will run when hitting source
require(bayestestR)
describe_posterior(model.traits4)
variables(model.traits4)
library(emmeans)



#variance for ID
var.ID<- posterior_samples(model.traits4)$"sd_ID__Intercept"^2
mean(var.ID);HPDinterval(as.mcmc(var.ID),0.95)#this computes CrI 
#residual variance
var.res<-posterior_samples(model.traits4)$"sigma"^2
mean(var.res);HPDinterval(as.mcmc(var.res),0.95)#this computes CrI

#repeatability estimate
repeatability<-var.ID/(var.ID + var.res)
describe_posterior(repeatability)

#extact posterior samples
o1 <- as.data.frame(model.traits4)
o1<-o1%>%select(9:117)
o1<-as.mcmc(o1)
#load survival dataframe
survival_final<-read_csv("Data/survival_dataframe.csv")
survival_final<-survival_final%>%group_by(Age)%>% filter(Age=="A")
require(bayestestR)


#Georgia----
require(tidybayes)
# 4a. edge + Harvest----
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="GA")
# Filter for individuals that survived or were harvested
harv_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Harvested")


#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, harv_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

harv_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)


#### GA-Log-odd ratio for the effect of harvest on risk-taking (distance to edge)----
# plot(harv_ga_da)#suppress plot so script will run when hitting source
harv_ga_da_sur<-data.frame(harv_ga_da)
describe_posterior(harv_ga_da_sur)



# 4b. edge + predation----
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="GA")
# Filter for individuals that survived or were predated
pred_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Predated")


#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, pred_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

pred_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)


#### GA-Log-odd ratio for the effect of predation on risk-taking (distance to edge)----
# plot(pred_ga_da)#suppress plot so script will run when hitting source
pred_ga_da_sur<-data.frame(pred_ga_da) 
describe_posterior(pred_ga_da_sur)


#South Carolina----
require(tidybayes)
# 4a. edge + Harvest----
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="Webb")
# Filter for individuals that survived or were harvested
harv_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Harvested")


#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, harv_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

harv_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)


#### SC-Log-odd ratio for the effect of harvest on risk-taking (distance to edge)----
# plot(harv_ga_da)#suppress plot so script will run when hitting source
harv_ga_da_sur<-data.frame(harv_ga_da)
describe_posterior(harv_ga_da_sur)



# 4b. edge + predation----
s_ga<-survival_final%>%group_by(Site)%>% filter(Site=="Webb")
# Filter for individuals that survived or were predated
pred_ga<-s_ga%>% filter(Year_Status=="Live" | Year_Status=="Predated")


#for loop that iterates through entire distribution of blup values
#results in 2000 unique glms that correspond to each BLUP value extracted from univariate models
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(
    as_tibble(o1),
    contains("r_ID")
  )[iter, ] 
  blups_FR <- tibble(
    ID = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  
  data_FR2_lm <- merge(blups_FR, pred_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

pred_ga_da <- as.mcmc(
  sapply(
    FR2_lm_list,
    function(x) {
      summary(x)$coefficients[2, 1]
    }
  )
)


#### SC-Log-odd ratio for the effect of predation on risk-taking (distance to edge)----
# plot(pred_ga_da)#suppress plot so script will run when hitting source
pred_ga_da_sur<-data.frame(pred_ga_da) 
describe_posterior(pred_ga_da_sur)



