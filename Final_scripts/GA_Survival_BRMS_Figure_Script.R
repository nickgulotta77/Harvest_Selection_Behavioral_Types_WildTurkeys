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
# require(posterior)

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
d1<- read_csv("data/Male_gps_Webb.csv")#load data
d2<- read_csv("data/Male_gps_GA_sites.csv") #load data

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
# to account for error in every cardinal direction we picked 75m as a reasonable distance
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
#1. Risk-taking (distance to hunter access) univariate model to extract BLUPs----
#run both populations together to obtain unbiased BLUP values
model.traits1<-brm(Dist_access_c~ Stage*Site + Site + (1|ID),
                  data= data, family="gaussian",
                  warmup=300, iter=15300,
                  chains=2,thin= 15, init="random",
                  cores = 32,
                  seed=bayes_seed,
                  backend = "rstan",
                  control = list(adapt_delta = 0.99, max_treedepth=15)) # within chain parallelization, number of threads + number of cores should match total cores you intend to use) 

#model summary
summary(model.traits1)
#check model covergence via trace/denisty plots
plot(model.traits1)
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
library(emmeans)

#extact posterior samples
o1 <- as.data.frame(model.traits1) #extracts all posterior samples
o1<-o1%>%select(9:117) #extract posterior samples of BLUP values for each individual
o1<-as.mcmc(o1) #make mcmc object

#load survival dataframe
survival_final<-read_csv("Data/survival_dataframe.csv") 
survival_final<-survival_final%>%group_by(Age)%>% filter(Age=="A") #subset to only include adults
require(bayestestR)


# 1a. Human influence on survival----

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



#### Log-odd ratio for the effect of human harvest on risk-taking (distance to hunter access)----
plot(harv_ga_da) #plot posterior distribution of log-odd ratios
harv_ga_da_sur<-data.frame(harv_ga_da) #make coefficients data frame so they can viewed using function below
describe_posterior(harv_ga_da_sur) # results presented here are sign reversed to produce plots.
                                   #  Multiple by -1 to ease interpretation, so positive log-odds = riskier individuals have higher survival

###extracts all blup values for plot
###each column is an ID, 1000 rows per ID
library(data.table)
FR <- select(
  as_tibble(o1),
  contains("r_ID")
)



#converts from wide to long format
require(tidyr)
FR_long<-
  FR %>%
  pivot_longer(
    everything(),
    names_to=c("ID")
  )

FR_long<-tibble(FR_long)
#removes ID. from column
FR_long$ID<-str_remove_all(FR_long$ID, "r_ID")
FR_long<-tibble(FR_long)

#mean blup value for each individual
FR<-FR_long%>%group_by(ID)%>%
  summarise(BLUP=mean(value))

#join survival data frame and BLUPs
FR<-left_join(harv_ga,FR)


## Create a new data frame to store predictions
predictions_df <- data.frame()

# Loop over the list of models to generate predictions
for (i in seq_along(FR2_lm_list)) {
  model <- FR2_lm_list[[i]]
  
  # Generate predictions for the current model
  data_FR2_lm$prediction <- predict(model, newdata = data_FR2_lm, type = "response")
  
  # Add the predictions to the predictions data frame
  predictions_df <- rbind(predictions_df, data.frame(model = rep(i, nrow(data_FR2_lm)), data_FR2_lm))
}

#back transform scaled and centered values to get raw distance values 
predictions_df$BLUP_raw <-
  predictions_df$blups_FR*744+825  #back transform

#make data frame with mean and 95% CrI for plot
df_summary <- predictions_df %>%
  group_by(ID) %>%
  summarise(mean = mean(prediction, na.rm = TRUE),
            ci_lower = quantile(prediction,0.025),
            ci_upper = quantile(prediction,0.975),
            BLUP_raw= mean(BLUP_raw))

#select data for plot
dr<-FR %>% select(ID, BLUP,Survival)
df_summary1<-left_join(df_summary, dr)
df_summary1$BLUP_r <-
  df_summary1$BLUP*744+825


#plot for impact of human harvest on distance to access in Georgia 
require(ggplot2)
p1<-ggplot(df_summary1, aes(x = BLUP_raw, y = mean)) +
  geom_line(color="cornsilk4", linewidth=5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),fill="cornsilk4",linetype=0, alpha = 0.2, linewidth=0) +
  geom_point(aes(x=BLUP_r, y= Survival), color="cornsilk4",alpha=0.6,size=8)+ theme_turkey()+
  xlab("Distance to Public Access (m)")+
  ylab("Survival Probability")+
  theme_turkey()+
  theme(axis.text=element_text(size=14),  
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+ 
  scale_x_continuous(breaks = c(500,1500,2500,3500))+
  ggtitle(label ="A.",
          subtitle = "") +  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(), 
        plot.margin = margin(r=45, l=20, b=35, t=35),
        panel.background = element_blank(),
        
        panel.grid.minor = element_blank())
p1





# 1b. Predation influence on survival----


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


#### Log-odd ratio for the effect of predation on risk-taking (distance to hunter access)----
plot(pred_ga_da) #plot posterior distribution of log-odd ratios
pred_ga_da_sur<-data.frame(pred_ga_da) #make posterior a data frame to compute median and 95% CrI
describe_posterior(pred_ga_da_sur)# results presented here are sign reversed to produce plots.
                                  #  Multiple by -1 to ease interpretation, so positive log-odds = riskier individuals have higher survival

###extracts all blup values
###each column is an ID, 1000 rows per ID
library(data.table)
FR <- select(
  as_tibble(o1),
  contains("r_ID")
)



#converts from wide to long format
require(tidyr)
FR_long<-
  FR %>%
  pivot_longer(
    everything(),
    names_to=c("ID")
  )

FR_long<-tibble(FR_long)

#removes ID. from column
FR_long$ID<-str_remove_all(FR_long$ID, "r_ID")
FR_long<-tibble(FR_long)

#mean blup value for each individual
FR<-FR_long%>%group_by(ID)%>%
  summarise(BLUP=mean(value))

#join survival data frame and BLUPs
FR<-left_join(pred_ga,FR)


## Create a new data frame to store predictions
predictions_df <- data.frame()

# Loop over the list of models to generate predictions
for (i in seq_along(FR2_lm_list)) {
  model <- FR2_lm_list[[i]]
  
  # Generate predictions for the current model
  data_FR2_lm$prediction <- predict(model, newdata = data_FR2_lm, type = "response")
  
  # Add the predictions to the predictions data frame
  predictions_df <- rbind(predictions_df, data.frame(model = rep(i, nrow(data_FR2_lm)), data_FR2_lm))
}

#back transform
predictions_df$BLUP_raw <-
  predictions_df$blups_FR*744+825 


df_summary <- predictions_df %>%
  group_by(ID) %>%
  summarise(mean = mean(prediction, na.rm = TRUE),
            ci_lower = quantile(prediction,0.025),
            ci_upper = quantile(prediction,0.975),
            BLUP_raw= mean(BLUP_raw))

dr<-FR %>% select(ID, BLUP,Survival)
df_summary1<-left_join(df_summary, dr)
df_summary1$BLUP_r <-
  df_summary1$BLUP*744+825

#plot for the effect of distance to public access (risk-taking) on survival from predation in Georgia
require(ggplot2)
p4<-ggplot(df_summary1, aes(x = BLUP_raw, y = mean)) +
  geom_line(color="cornsilk4", linewidth=5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),fill="cornsilk4",linetype=0, alpha = 0.2, linewidth=0) +
  geom_point(aes(x=BLUP_raw, y= Survival), color="cornsilk4",alpha=0.6,size=8)+ theme_turkey()+
  xlab("Distance to Public Access (m)")+
  ylab("Survival Probability")+
  theme_turkey()+
  theme(axis.text=element_text(size=14),  
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+ 
  scale_x_continuous(breaks = c(500,1500,2500,3500))+
  ggtitle(label ="E.",
          subtitle = "") +  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(), 
        plot.margin = margin(r=45, l=20, b=35, t=35),
        panel.background = element_blank(),
        
        panel.grid.minor = element_blank())
p4



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
#plot markov chains and check model convergence
plot(model.traits2)
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


# 2a. Effect of survival from harvest on exploration (speed)----
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



####GA human distance to access + survival----
plot(harv_ga_da)#plot posterior distribution of log-odd ratios of survival
harv_ga_da_sur<-data.frame(harv_ga_da) #make data frame to compute survival estimates
describe_posterior(harv_ga_da_sur)

###extracts all blup values
###each column is an ID, 1000 rows per ID
library(data.table)
FR <- select(
  as_tibble(o1),
  contains("r_ID")
)



#converts from wide to long format
require(tidyr)
FR_long<-
  FR %>%
  pivot_longer(
    everything(),
    names_to=c("ID")
  )


FR_long<-tibble(FR_long)
#removes ID. from column
FR_long$ID<-str_remove_all(FR_long$ID, "r_ID")
FR_long<-tibble(FR_long)

#mean blup value for each individual
FR<-FR_long%>%group_by(ID)%>%
  summarise(BLUP=mean(value))

#join survival data frame and BLUPs
FR<-left_join(harv_ga,FR)


## Create a new data frame to store predictions
predictions_df <- data.frame()

# Loop over the list of models to generate predictions
for (i in seq_along(FR2_lm_list)) {
  model <- FR2_lm_list[[i]]
  
  # Generate predictions for the current model
  data_FR2_lm$prediction <- predict(model, newdata = data_FR2_lm, type = "response")
  
  # Add the predictions to the predictions data frame
  predictions_df <- rbind(predictions_df, data.frame(model = rep(i, nrow(data_FR2_lm)), data_FR2_lm))
}

predictions_df$BLUP_raw <-
  predictions_df$blups_FR*98.1+198  #bac

df_summary <- predictions_df %>%
  group_by(ID) %>%
  summarise(mean = mean(prediction, na.rm = TRUE),
            ci_lower = quantile(prediction,0.025),
            ci_upper = quantile(prediction,0.975),
            BLUP_raw= mean(BLUP_raw))

dr<-FR %>% select(ID, BLUP,Survival)
df_summary1<-left_join(df_summary, dr)
df_summary1$BLUP_r <-
  df_summary1$BLUP*98.1+198

#plot for effect of harvest on exploration 
require(ggplot2)
p2<-ggplot(df_summary1, aes(x = BLUP_raw, y = mean)) +
  geom_line(color="firebrick4", linewidth=5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),fill="firebrick4",linetype=0, alpha = 0.25, linewidth=0) +
  geom_point(aes(x=BLUP_r, y= Survival), color="firebrick4",alpha=0.6,size=8)+ theme_turkey()+
  xlab("Average speed (m/hr)")+
  ylab("Survival Probability")+
  theme_turkey()+
  theme(axis.text=element_text(size=14),  
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+ 
  scale_x_continuous(breaks = c(100,150,200,250,300))+
  ggtitle(label ="B.",
          subtitle = "") +  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(), 
        plot.margin = margin(r=45, l=20, b=35, t=35),
        panel.background = element_blank(),
        
        panel.grid.minor = element_blank())
p2







# 2b. Effect of survival from predation on exploration (speed)----
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



#### Log-odd ratio for the effect of predation on exploration (speed)----
plot(pred_ga_da) #posterior distribution of log odd ratios for survival
pred_ga_da_sur<-data.frame(pred_ga_da) #make data frame to compute median, 95% CrI, and Rhat values
describe_posterior(pred_ga_da_sur)



###extracts all blup values
###each column is an ID, 1000 rows per ID
library(data.table)
FR <- select(
  as_tibble(o1),
  contains("r_ID")
)



#converts from wide to long format
require(tidyr)
FR_long<-
  FR %>%
  pivot_longer(
    everything(),
    names_to=c("ID")
  )


FR_long<-tibble(FR_long)
#removes ID. from column
FR_long$ID<-str_remove_all(FR_long$ID, "r_ID")
FR_long<-tibble(FR_long)

#mean blup value for each individual
FR<-FR_long%>%group_by(ID)%>%
  summarise(BLUP=mean(value))

#join survival data frame and BLUPs
FR<-left_join(pred_ga,FR)


## Create a new data frame to store predictions
predictions_df <- data.frame()

# Loop over the list of models to generate predictions
for (i in seq_along(FR2_lm_list)) {
  model <- FR2_lm_list[[i]]
  
  # Generate predictions for the current model
  data_FR2_lm$prediction <- predict(model, newdata = data_FR2_lm, type = "response")
  
  # Add the predictions to the predictions data frame
  predictions_df <- rbind(predictions_df, data.frame(model = rep(i, nrow(data_FR2_lm)), data_FR2_lm))
}

predictions_df$BLUP_raw <-
  predictions_df$blups_FR*98.1+198  #back transform

df_summary <- predictions_df %>%
  group_by(ID) %>%
  summarise(mean = mean(prediction, na.rm = TRUE),
            ci_lower = quantile(prediction,0.025),
            ci_upper = quantile(prediction,0.975),
            BLUP_raw= mean(BLUP_raw))

dr<-FR %>% select(ID, BLUP,Survival)
df_summary1<-left_join(df_summary, dr)
df_summary1$BLUP_r <-
  df_summary1$BLUP*98.1+198


#plot for the effect of exploratory behavior on survival from predation
require(ggplot2)
p5<-ggplot(df_summary1, aes(x = BLUP_raw, y = mean)) +
  geom_line(color="firebrick4", linewidth=5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),fill="firebrick4",linetype=0, alpha = 0.25, linewidth=0) +
  geom_point(aes(x=BLUP_raw, y= Survival), color="firebrick4",alpha=0.6,size=8)+ theme_turkey()+
  xlab("Average speed (m/hr)")+
  ylab("Survival Probability")+
  theme_turkey()+
  theme(axis.text=element_text(size=14),  
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+ 
  scale_x_continuous(breaks = c(100,150,200,250,300))+
  ggtitle(label ="F.",
          subtitle = "") +  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(), 
        plot.margin = margin(r=45, l=20, b=35, t=35),
        panel.background = element_blank(),
        
        panel.grid.minor = element_blank())
p5




#3. Risk-taking (distnce to open landcover) univariate model to extract BLUPs----
#run both populations together to obtain unbiased BLUP values
model.traits3<-brm(Open_c~ Stage*Site + Site +  (1|ID) ,
                  data= data, family="gaussian",
                  warmup=300, iter=10300,
                  chains=2,thin= 10, init="random",
                  cores = 12,
                  seed=bayes_seed,
                  backend = "rstan",
                  control = list(adapt_delta = 0.99, max_treedepth=15)) # within chain parallelization, number of threads + number of cores should match total cores you intend to use) 

#summary of model
summary(model.traits3)
plot(model.traits)#plot markov chains and density plots
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


#### Log-odd ratio for the effect of harvest on risk-taking (distance to open)----
plot(harv_ga_da)
harv_ga_da_sur<-data.frame(harv_ga_da) 
describe_posterior(harv_ga_da_sur)# results presented here are sign reversed to produce plots.
                                  #  Multiple by -1 to ease interpretation, so positive log-odds = riskier individuals have higher survival



###extracts all blup values
###each column is an ID, 1000 rows per ID
library(data.table)
FR <- select(
  as_tibble(o1),
  contains("r_ID")
)



#converts from wide to long format
require(tidyr)
FR_long<-
  FR %>%
  pivot_longer(
    everything(),
    names_to=c("ID")
  )

FR_long<-tibble(FR_long)

#removes ID. from column
FR_long$ID<-str_remove_all(FR_long$ID, "r_ID")
FR_long<-tibble(FR_long)

#mean blup value for each individual
FR<-FR_long%>%group_by(ID)%>%
  summarise(BLUP=mean(value))

#join survival data frame and BLUPs
FR<-left_join(harv_ga,FR)


## Create a new data frame to store predictions
predictions_df <- data.frame()

# Loop over the list of models to generate predictions
for (i in seq_along(FR2_lm_list)) {
  model <- FR2_lm_list[[i]]
  
  # Generate predictions for the current model
  data_FR2_lm$prediction <- predict(model, newdata = data_FR2_lm, type = "response")
  
  # Add the predictions to the predictions data frame
  predictions_df <- rbind(predictions_df, data.frame(model = rep(i, nrow(data_FR2_lm)), data_FR2_lm))
}

predictions_df$BLUP_raw <-
  predictions_df$blups_FR*212+276  #back transform

df_summary <- predictions_df %>%
  group_by(ID) %>%
  summarise(mean = mean(prediction, na.rm = TRUE),
            ci_lower = quantile(prediction,0.025),
            ci_upper = quantile(prediction,0.975),
            BLUP_raw= mean(BLUP_raw))

dr<-FR %>% select(ID, BLUP,Survival)
df_summary1<-left_join(df_summary, dr)
df_summary1$BLUP_r <-
  df_summary1$BLUP*212+276


#plot the effect of risk-taking on survival from harvest
require(ggplot2)
p3<-ggplot(df_summary1, aes(x = BLUP_raw, y = mean)) +
  geom_line(color="gold4", linewidth=5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),fill="gold4",linetype=0, alpha = 0.25, linewidth=0) +
  geom_point(aes(x=BLUP_r, y= Survival), color="gold4",alpha=0.6,size=8)+ theme_turkey()+
  xlab("Distance to open landcover (m)")+
  ylab("Survival Probability")+
  theme_turkey()+
  theme(axis.text=element_text(size=14),  
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+ 
  scale_x_continuous(breaks = c(200,400,600,800))+
  ggtitle(label ="C.",
          subtitle = "") +  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(), 
        plot.margin = margin(r=45, l=20, b=35, t=35),
        panel.background = element_blank(),
        
        panel.grid.minor = element_blank())
p3





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


####distance to open + survival----
plot(pred_ga_da)
pred_ga_da_sur<-data.frame(pred_ga_da) 
describe_posterior(pred_ga_da_sur)# results presented here are sign reversed to produce plots.
                                  #  Multiple by -1 to ease interpretation, so positive log-odds = riskier individuals have higher survival




###extracts all blup values
###each column is an ID, 1000 rows per ID
library(data.table)
FR <- select(
  as_tibble(o1),
  contains("r_ID")
)



#converts from wide to long format
require(tidyr)
FR_long<-
  FR %>%
  pivot_longer(
    everything(),
    names_to=c("ID")
  )


FR_long<-tibble(FR_long)
#removes ID. from column
FR_long$ID<-str_remove_all(FR_long$ID, "r_ID")
FR_long<-tibble(FR_long)

#mean blup value for each individual
FR<-FR_long%>%group_by(ID)%>%
  summarise(BLUP=mean(value))

#join survival data frame and BLUPs
FR<-left_join(pred_ga,FR)


## Create a new data frame to store predictions
predictions_df <- data.frame()

# Loop over the list of models to generate predictions
for (i in seq_along(FR2_lm_list)) {
  model <- FR2_lm_list[[i]]
  
  # Generate predictions for the current model
  data_FR2_lm$prediction <- predict(model, newdata = data_FR2_lm, type = "response")
  
  # Add the predictions to the predictions data frame
  predictions_df <- rbind(predictions_df, data.frame(model = rep(i, nrow(data_FR2_lm)), data_FR2_lm))
}

predictions_df$BLUP_raw <-
  predictions_df$blups_FR*212+276  #bac

df_summary <- predictions_df %>%
  group_by(ID) %>%
  summarise(mean = mean(prediction, na.rm = TRUE),
            ci_lower = quantile(prediction,0.025),
            ci_upper = quantile(prediction,0.975),
            BLUP_raw= mean(BLUP_raw))

dr<-FR %>% select(ID, BLUP,Survival)
df_summary1<-left_join(df_summary, dr)
df_summary1$BLUP_r <-
  df_summary1$BLUP*212+276

#plot for the effect of risk-taking on survival from predation
require(ggplot2)
p6<-ggplot(df_summary1, aes(x = BLUP_raw, y = mean)) +
  geom_line(color="gold4", linewidth=5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),fill="gold4",linetype=0, alpha = 0.25, linewidth=0) +
  geom_point(aes(x=BLUP_r, y= Survival), color="gold4",alpha=0.6,size=8)+ theme_turkey()+
  xlab("Distance to open landcover (m)")+
  ylab("Survival Probability")+
  theme_turkey()+
  theme(axis.text=element_text(size=14),  
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+ 
  scale_x_continuous(breaks = c(200,400,600,800))+
  ggtitle(label ="G.",
          subtitle = "") +  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(), 
        plot.margin = margin(r=45, l=20, b=35, t=35),
        panel.background = element_blank(),
        
        panel.grid.minor = element_blank())
p6



#4. Risk-taking (distance to edge landcover) univariate model to extract BLUPs----
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
plot(model.traits4)
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



####GA human distance to edge + survival----
plot(harv_ga_da)
harv_ga_da_sur<-data.frame(harv_ga_da) # results presented here are sign reversed to produce plots.
                                       #  Multiple by -1 to ease interpretation, so positive log-odds = riskier individuals have higher survival

describe_posterior(harv_ga_da_sur)

###extracts all blup values
###each column is an ID, 1000 rows per ID
library(data.table)
FR <- select(
  as_tibble(o1),
  contains("r_ID")
)



#converts from wide to long format
require(tidyr)
FR_long<-
  FR %>%
  pivot_longer(
    everything(),
    names_to=c("ID")
  )


FR_long<-tibble(FR_long)
#removes ID. from column
FR_long$ID<-str_remove_all(FR_long$ID, "r_ID")
FR_long<-tibble(FR_long)

#mean blup value for each individual
FR<-FR_long%>%group_by(ID)%>%
  summarise(BLUP=mean(value))

#join survival data frame and BLUPs
FR<-left_join(harv_ga,FR)


## Create a new data frame to store predictions
predictions_df <- data.frame()

# Loop over the list of models to generate predictions
for (i in seq_along(FR2_lm_list)) {
  model <- FR2_lm_list[[i]]
  
  # Generate predictions for the current model
  data_FR2_lm$prediction <- predict(model, newdata = data_FR2_lm, type = "response")
  
  # Add the predictions to the predictions data frame
  predictions_df <- rbind(predictions_df, data.frame(model = rep(i, nrow(data_FR2_lm)), data_FR2_lm))
}

predictions_df$BLUP_raw <-
  predictions_df$blups_FR*122+133  #bac

df_summary <- predictions_df %>%
  group_by(ID) %>%
  summarise(mean = mean(prediction, na.rm = TRUE),
            ci_lower = quantile(prediction,0.025),
            ci_upper = quantile(prediction,0.975),
            BLUP_raw= mean(BLUP_raw))

dr<-FR %>% select(ID, BLUP,Survival)
df_summary1<-left_join(df_summary, dr)
df_summary1$BLUP_r <-
  df_summary1$BLUP*122+133 

#plot for the effect of risk-taking (distance to edge) on surviavl from harvest
require(ggplot2)
p33<-ggplot(df_summary1, aes(x = BLUP_raw, y = mean)) +
  geom_line(color="#333300", linewidth=5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),fill="#333300",linetype=0, alpha = 0.25, linewidth=0) +
  geom_point(aes(x=BLUP_r, y= Survival), color="#333300",alpha=0.6,size=8)+ theme_turkey()+
  xlab("Distance to edge landcover (m)")+
  ylab("Survival Probability")+
  theme_turkey()+
  theme(axis.text=element_text(size=14),  
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+ 
  scale_x_continuous(breaks = c(100,200,300,400))+
  ggtitle(label ="D.",
          subtitle = "") +  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(), 
        plot.margin = margin(r=45, l=20, b=35, t=35),
        panel.background = element_blank(),
        
        panel.grid.minor = element_blank())
p33




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


####GA predators distance to access + survival----
plot(pred_ga_da)
pred_ga_da_sur<-data.frame(pred_ga_da) # results presented here are sign reversed to produce plots.
                                       #  Multiple by -1 to ease interpretation, so positive log-odds = riskier individuals have higher survival

describe_posterior(pred_ga_da_sur)



###extracts all blup values
###each column is an ID, 1000 rows per ID
library(data.table)
FR <- select(
  as_tibble(o1),
  contains("r_ID")
)



#converts from wide to long format
require(tidyr)
FR_long<-
  FR %>%
  pivot_longer(
    everything(),
    names_to=c("ID")
  )


FR_long<-tibble(FR_long)
#removes ID. from column
FR_long$ID<-str_remove_all(FR_long$ID, "r_ID")
FR_long<-tibble(FR_long)

#mean blup value for each individual
FR<-FR_long%>%group_by(ID)%>%
  summarise(BLUP=mean(value))

#join survival data frame and BLUPs
FR<-left_join(pred_ga,FR)


## Create a new data frame to store predictions
predictions_df <- data.frame()

# Loop over the list of models to generate predictions
for (i in seq_along(FR2_lm_list)) {
  model <- FR2_lm_list[[i]]
  
  # Generate predictions for the current model
  data_FR2_lm$prediction <- predict(model, newdata = data_FR2_lm, type = "response")
  
  # Add the predictions to the predictions data frame
  predictions_df <- rbind(predictions_df, data.frame(model = rep(i, nrow(data_FR2_lm)), data_FR2_lm))
}

predictions_df$BLUP_raw <-
  predictions_df$blups_FR*122+133  #bac

df_summary <- predictions_df %>%
  group_by(ID) %>%
  summarise(mean = mean(prediction, na.rm = TRUE),
            ci_lower = quantile(prediction,0.025),
            ci_upper = quantile(prediction,0.975),
            BLUP_raw= mean(BLUP_raw))

dr<-FR %>% select(ID, BLUP,Survival)
df_summary1<-left_join(df_summary, dr)
df_summary1$BLUP_r <-
  df_summary1$BLUP*122+133 

#plot for the effect of risk-taking (distance to edge) on survival from predation
require(ggplot2)
p66<-ggplot(df_summary1, aes(x = BLUP_raw, y = mean)) +
  geom_line(color="#333300", linewidth=5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),fill="#333300",linetype=0, alpha = 0.25, linewidth=0) +
  geom_point(aes(x=BLUP_r, y= Survival), color="#333300",alpha=0.6,size=8)+ theme_turkey()+
  xlab("Distance to edge landcover (m)")+
  ylab("Survival Probability")+
  theme_turkey()+
  theme(axis.text=element_text(size=14),  
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+ 
  scale_x_continuous(breaks = c(100,200,300,400))+
  ggtitle(label ="H.",
          subtitle = "") +  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(), 
        plot.margin = margin(r=45, l=20, b=35, t=35),
        panel.background = element_blank(),
        
        panel.grid.minor = element_blank())
p66


# Combine plots for final plot----
pharv<-cowplot::plot_grid(p1,p2,p3,p33, label_size = 73, ncol=4)
title1<-"Human induced selection on behavioral types"
require(ggpubr)
require(gridExtra)
require(cowplot)
title_grob<-text_grob(title1, face="bold", size=72)
pharv<-gridExtra::grid.arrange(pharv,top=title_grob)


ppred<-cowplot::plot_grid(p4,p5,p6,p66, label_size = 73, ncol=4)
title1<-"Predation induced selection on behavioral types"
require(ggpubr)
title_grob<-text_grob(title1, face="bold", size=72)
ppred<-gridExtra::grid.arrange(ppred, top=title_grob)


#final survival plot
cowplot::plot_grid(pharv, ppred, nrow = 2)
