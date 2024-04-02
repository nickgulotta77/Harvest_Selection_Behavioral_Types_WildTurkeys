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


# Bayes stuff
# Use the cmdstanr backend for Stan because it's faster and more modern than
# the default rstan. You need to install the cmdstanr package first
# (https://mc-stan.org/cmdstanr/) and then run cmdstanr::install_cmdstan() to
# install cmdstan on your computer.
options(mc.cores = 12,  # Use 4 cores
        brms.backend = "rstan")
bayes_seed <- 1234


# import data and add tag ID
library(readr)
d1<- read_csv("data/Male_gps_Webb.csv")
d2<- read_csv("data/Male_gps_GA_sites.csv")

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
data$Age<-as.factor(data$Age)
data$Year<-as.factor(data$Year)

# standardize traits so models converge more quickly 
# and are easier to compare distances to speed for subsequent analysis
data$mean_speed_c<-scale(data$mean_speed,scale = T, center = T) #proxy for exploration
data$Dist_Access_c<-scale(data$Dist_Access,scale = T, center = T)#risk-taking
data$Open_c<-scale(data$Open, scale = T, center = T)#risk-taking
data$edge_c<-scale(data$dist_edge, scale = T, center = T)#risk-taking
# histogram of each behavioral trait 
hist(data$mean_speed_c)
hist(data$Open_c)
hist(data$Dist_Access_c)



# how many tests per hunting stage?
table(data$Stage, data$ID)
po<-data%>%
  group_by(Stage)%>%
  count()

barplot(po$n~ po$Stage)

# how many tests per individual?
po<-data%>%
  group_by(ID,Year, Stage)%>%
  count()

counts<-po%>%group_by(ID)%>%
  summarise(sum(n))

max(counts$`sum(n)`)
median(counts$`sum(n)`)
min(counts$`sum(n)`)

# how many individuals in each stage?


counts<-data%>%group_by(Site,Stage,ID)%>%
  count()
counts<-counts%>%group_by(Site,ID,Stage)%>%
  summarise(sum(n))
barplot(table(counts$Site, counts$Stage))
table(counts$Site, counts$Stage)

###**Univarite models for each trait**----


#1. Distance to access (risk-taking)----
# Individual differences in BT
require(brms)

m222_brm<-brm(Dist_Access_c~ Stage*Site + Site + (1|ID) ,
              data= data, family="gaussian",
              warmup=1000, iter=10000,
              chains=2, init="random",
              cores = 32,
              seed=bayes_seed,
              backend = "rstan",
              control = list(adapt_delta = 0.99, max_treedepth=15)) # within chain parallelization, number of threads + number of cores should match total cores you intend to use) 
summary(m222_brm)
plot(m222_brm)
require(bayestestR)
describe_posterior(m222_brm)
variables(m222_brm)
library(emmeans)
#pair wise differences for hunting stage
(group_diff2 <- emmeans(m222_brm, pairwise~Site ))
(group_diff2 <- emmeans(m222_brm, pairwise~Stage*Site ))

describe_posterior(group_diff2)
#variance for ID
var.ID<- posterior_samples(m222_brm)$"sd_ID__Intercept"^2
mean(var.ID);HPDinterval(as.mcmc(var.ID),0.95)#this computes CrI 
#residual variance
var.res<-posterior_samples(m222_brm)$"sigma"^2
mean(var.res);HPDinterval(as.mcmc(var.res),0.95)#this computes CrI

#repeatability estimate
repeatability<-var.ID/(var.ID + var.res)
describe_posterior(repeatability)


#conditional effect plots
me<-conditional_effects(m222_brm, "Stage:Site")
plot(me,theme = theme_classic()) 
me<-conditional_effects(m222_brm, "Site")
plot(me,theme = theme_classic()) 


##Site plot

nn<-m222_brm %>% 
  spread_draws(b_Intercept, b_SiteWebb) %>%
  head(20000)%>%
  mutate(Webb_mean = (b_Intercept + b_SiteWebb))

data9<- data.frame(nn)
data9<-data9%>%select(b_Intercept, Webb_mean)%>%
  mutate(Georgia = b_Intercept,
         Webb = Webb_mean)%>%
  select(Georgia, Webb)

#convert data to long format for plotting
require(tidyr)
d1<-data9%>% 
  pivot_longer(c("Georgia", "Webb") , 
               names_to= "Site", values_to= "Estimate")

require(dplyr)
d11<- d1 %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(mean = mean(Estimate))%>%
  dplyr::ungroup()

d11$Estimate <-
  d11$Estimate*744 + 825 #back transfrom


#plot for figure 1
plot_1a<-ggplot(data=d11,aes(y = Site, x = Estimate, fill=Site)) +
  labs(y = "Study Site",
       x = "Distance to hunter access (m)",
       element_text(face="bold"))+ 
  scale_fill_manual(values = c("#993333",
                               "#999933"))+
  stat_halfeye(point_interval = median_qi,scale=0.9, 
               slab_alpha= 1, point_size=24,stroke=1, interval_alpha=1,
               .width = c(0.67, 0.95),
               interval_size_range = c(5.5, 9.5))+
  scale_x_continuous(breaks = c(400,600,800,1000,1200))+
  theme(plot.title = element_text(face = "bold", hjust=1.5),
        plot.subtitle=element_text(hjust=0.9))+
  theme_turkey()+
  ggtitle(label ="A.",
          subtitle = "")+theme(legend.position="none")+  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 74, face = "bold"),
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
        
        panel.grid.minor = element_blank()) + scale_y_discrete(labels = c('Georgia','South Carolina'))
plot_1a



#plot site vs hunting stage
me<-conditional_effects(m222_brm, "Stage:Site")

dp<-data.frame(me[["Stage:Site"]])
dp$Site<-as.factor(dp$Site)
dp$estimate__ <-
  dp$estimate__*744 + 825  #backtransform
dp$lower__ <-
  dp$lower__*744 + 825 #backtransform
dp$upper__ <-
  dp$upper__*744 + 825  #backtransform

#plot for figure 1
plot_1<-ggplot(dp, aes(x = Stage, y = estimate__, group=Site, color= Site)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__), width = .5,
                position=position_dodge(width=0.65),size=4.85) +
  geom_point(position=position_dodge(width=0.65),size=24) +
  scale_y_continuous("Distance to hunter access (m)", breaks= c(400,600,800,1000,1200)) + theme_classic() + theme_turkey()+
  scale_x_discrete(labels = c('Pre-Hunt','Hunt','Post-Hunt'))+
  theme(legend.position="none")+
  scale_color_manual(values = c("#993333",
                                "#999933"))+
  ggtitle(label ="B.",
          subtitle = "")+theme(legend.position="none")+  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 74, face = "bold"),
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

plot_1


#2. Exploration (mean-daily speed)----
# Individual differences in BT
require(brms)

m223_brm<-brm(mean_speed_c~ Stage*Site  +Site + (1|ID) ,
              data= data, family="gaussian",
              warmup=1000, iter=10000, thin=4,
              chains=4, init="random",
              cores = 12,
              seed=bayes_seed,
              backend = "rstan",
              control = list(adapt_delta = 0.99, max_treedepth=15)) # within chain parallelization, number of threads + number of cores should match total cores you intend to use) 
summary(m223_brm)
plot(m223_brm)
require(bayestestR)
describe_posterior(m223_brm)
variables(m223_brm)
library(emmeans)

#pair wise differences for nest stage
(group_diff2 <- emmeans(m223_brm, pairwise~Site ))
(group_diff2 <- emmeans(m223_brm, pairwise~Stage*Site ))
describe_posterior(group_diff2)

#variance for ID
var.ID<- posterior_samples(m223_brm)$"sd_ID__Intercept"^2
mean(var.ID);HPDinterval(as.mcmc(var.ID),0.95)#this computes CrI 
#residual variance
var.res<-posterior_samples(m223_brm)$"sigma"^2
mean(var.res);HPDinterval(as.mcmc(var.res),0.95)#this computes CrI

#repeatability estimate
repeatability<-var.ID/(var.ID+var.res)
describe_posterior(repeatability)
#pair wise differences for nest stage
(group_diff2 <- emmeans(m223_brm, pairwise~Age ))
(group_diff2 <- emmeans(m223_brm, pairwise~Stage ))

#conditional effect plots
me2<-conditional_effects(m223_brm, "Stage:Site")
plot(me2,theme = theme_classic()) 
me2<-conditional_effects(m223_brm, "Site")
plot(me2,theme = theme_classic()) 



##Plot for figure 1

nn<-m223_brm %>%
  spread_draws(b_Intercept, b_SiteWebb) %>%
  head(20000)%>%
  mutate(Webb_mean = (b_Intercept + b_SiteWebb))

data9<- data.frame(nn)
data9<-data9%>%select(b_Intercept, Webb_mean)%>%
  mutate(Georgia = b_Intercept,
         Webb = Webb_mean)%>%
  select(Georgia, Webb)

#convert data to long format for plotting
require(tidyr)
d1<-data9%>% 
  pivot_longer(c("Georgia", "Webb") , 
               names_to= "Site", values_to= "Estimate")

require(dplyr)
d11<- d1 %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(mean = mean(Estimate))%>%
  dplyr::ungroup()


require(dplyr)
d11<- d1 %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(mean = mean(Estimate))%>%
  dplyr::ungroup()

d11$Estimate <-
  d11$Estimate*98.1 + 198 #back transform

#plot for figure 1
plot_2a<-ggplot(data=d11,aes(y = Site, x = Estimate, fill=Site)) +
stat_halfeye(point_interval = median_qi,scale=0.9, 
            slab_alpha= 1, point_size=24,stroke=1, interval_alpha=1,
            .width = c(0.67, 0.95),
            interval_size_range = c(5.5, 9.5))+
        labs(y = "Study Site",
        x = "Average speed (m/hr)",
       element_text(face="bold"))+ 
  scale_fill_manual(values = c("#993333",
                               "#999933"))+
  theme(plot.title = element_text(face = "bold", hjust=0.9), plot.subtitle=element_text(hjust=0.9))+
  scale_x_continuous(n.breaks = 6)+theme_turkey()+
  ggtitle(label ="C.",
          subtitle = "")+theme(legend.position="none")+theme(legend.position="none")+  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 74, face = "bold"),
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
        
        panel.grid.minor = element_blank()) + scale_y_discrete(labels = c('Georgia','South Carolina'))
plot_2a

#extract estimates for plot in figure 1
me2<-conditional_effects(m223_brm, "Stage:Site")
plot(me2,theme = theme_classic())
dp<-data.frame(me2[["Stage:Site"]])
dp$Site<-as.factor(dp$Site)
dp$estimate__ <-
  dp$estimate__*98.1 + 198  #back transform
dp$lower__ <-
  dp$lower__*98.1 + 198   #back transform
dp$upper__ <-
  dp$upper__*98.1 + 198   #back transform

#plot for figure 1
plot_2<-ggplot(dp, aes(x = Stage, y = estimate__, colour=Site, group=Site)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__), width = .5,
                position=position_dodge(width=0.65),size=4.85) +
  geom_point(position=position_dodge(width=0.65),size=24) +
  scale_y_continuous("Average speed (m/hr)", n.breaks = 6) + theme_classic() + theme_turkey()+
  scale_x_discrete(labels = c('Pre-Hunt','Hunt','Post-Hunt'))+
  theme(legend.position="none")+
  scale_color_manual(values = c("#993333",
                                "#999933")) +
  ggtitle(label ="D.",
          subtitle = "")+theme(legend.position="none")+  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 74, face = "bold"),
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
plot_2


#3. Risk-taking (Avg daily distance to open landcover)----
# Individual differences in BT
require(brms)

m224_brm<-brm(Open_c~ Stage*Site +Site + (1|ID) ,
              data= data, family="gaussian",
              warmup=1000, iter=10000,
              chains=2, init="random",
              cores = 12,
              seed=bayes_seed,
              backend = "rstan",
              control = list(adapt_delta = 0.99, max_treedepth=15)) # within chain parallelization, number of threads + number of cores should match total cores you intend to use) 
summary(m224_brm)
plot(m224_brm)
require(bayestestR)
describe_posterior(m224_brm)
variables(m224_brm)
library(emmeans)


#pair wise differences for nest stage
(group_diff2 <- emmeans(m224_brm, pairwise~Site ))
(group_diff2 <- emmeans(m224_brm, pairwise~Stage*Site ))
describe_posterior(group_diff2)
#variance for ID
var.ID<- posterior_samples(m224_brm)$"sd_ID__Intercept"^2
mean(var.ID);HPDinterval(as.mcmc(var.ID),0.95)#this computes CrI 
#residual variance
var.res<-posterior_samples(m224_brm)$"sigma"^2
mean(var.res);HPDinterval(as.mcmc(var.res),0.95)#this computes CrI

#repeatability estimate
repeatability<-var.ID/ (var.ID + var.res)
describe_posterior(repeatability)
#pair wise differences for nest stage
(group_diff2 <- emmeans(m224_brm, pairwise~Age ))
(group_diff2 <- emmeans(m224_brm, pairwise~Stage ))

#conditional effect plots
me3<-conditional_effects(m224_brm, "Stage:Site")
plot(me3,theme = theme_classic()) 
me3<-conditional_effects(m224_brm, "Site")
plot(me3,theme = theme_classic()) 


#Posterior required for plotting figure 1
nn<-m224_brm %>% spread_draws(b_Intercept, b_SiteWebb) %>%
  head(20000)%>%
  mutate(Webb_mean = (b_Intercept + b_SiteWebb))

data9<- data.frame(nn)
data9<-data9%>%select(b_Intercept, Webb_mean)%>%
  mutate(Georgia = b_Intercept,
         Webb = Webb_mean)%>%
  select(Georgia, Webb)

#convert data to long format for plotting
require(tidyr)
d1<-data9%>% 
  pivot_longer(c("Georgia", "Webb") , 
               names_to= "Site", values_to= "Estimate")

require(dplyr)
d11<- d1 %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(mean = mean(Estimate))%>%
  dplyr::ungroup()


require(dplyr)
d11<- d1 %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(mean = mean(Estimate))%>%
  dplyr::ungroup()

d11$Estimate <-
  d11$Estimate*212 + 276

#plot for figure 1
plot_3a<-ggplot(data=d11,aes(y = Site, x = Estimate, fill=Site)) +
  stat_halfeye(point_interval = median_qi,scale=0.9, 
               slab_alpha= 1, point_size=24,stroke=1, interval_alpha=1,
               .width = c(0.67, 0.95),
               interval_size_range = c(5.5, 9.5))+
  labs(y = "Study Site",
       x = "Distance to open landcover (m)",
       element_text(face="bold"))+ 
  scale_fill_manual(values = c("#993333",
                               "#999933"))+
  theme(plot.title = element_text(face = "bold", hjust=0.9), plot.subtitle=element_text(hjust=0.9))+
  scale_x_continuous(n.breaks = 6)+theme_turkey()+
  ggtitle(label ="E.",
          subtitle = "")+theme(legend.position="none")+theme(legend.position="none")+  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 74, face = "bold"),
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
        
        panel.grid.minor = element_blank()) + scale_y_discrete(labels = c('Georgia','South Carolina'))
plot_3a



#Hunting stage plot for figure 1
me3<-conditional_effects(m224_brm, "Stage:Site")

dp<-data.frame(me3[["Stage:Site"]])
dp$Site<-as.factor(dp$Site)
dp$estimate__ <-
  dp$estimate__*212 + 276  #back transform
dp$lower__ <-
  dp$lower__*212 + 276  #back transform
dp$upper__ <-
  dp$upper__*212 + 276  #back transform

#plot for figure 1
plot_3<-ggplot(dp, aes(x = Stage, y = estimate__, colour=Site, group=Site)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__), width = .5,
                position=position_dodge(width=0.65),size=4.85) +
  geom_point(position=position_dodge(width=0.65),size=24) +
  scale_y_continuous("Distance to open landcover (m)",n.breaks = 6)+
  theme_classic() + theme_turkey()+
  scale_x_discrete(labels = c('Pre-Hunt','Hunt','Post-Hunt'))+
  theme(legend.position="none")+
  scale_color_manual(values = c("#993333",
                                "#999933"))+
  ggtitle(label ="F.",
          subtitle = "")+theme(legend.position="none")+  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 74, face = "bold"),
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


plot_3



#4. Risk-taking (Avg daily distance to edge landcover)----
# Individual differences in BT
require(brms)

m225_brm<-brm(edge_c~ Stage*Site +Site + (1|ID) ,
              data= data, family="gaussian",
              warmup=1000, iter=10000,
              chains=2, init="random",
              cores = 12,
              seed=bayes_seed,
              backend = "rstan",
              control = list(adapt_delta = 0.99, max_treedepth=15)) # within chain parallelization, number of threads + number of cores should match total cores you intend to use) 
summary(m225_brm)
plot(m225_brm)
require(bayestestR)
describe_posterior(m225_brm)
variables(m225_brm)
library(emmeans)

#pair wise differences for nest stage
(group_diff2 <- emmeans(m225_brm, pairwise~Site ))
(group_diff2 <- emmeans(m225_brm, pairwise~Stage*Site ))
describe_posterior(group_diff2)

#variance for ID
var.ID<- posterior_samples(m225_brm)$"sd_ID__Intercept"^2
mean(var.ID);HPDinterval(as.mcmc(var.ID),0.95)#this computes CrI 
describe_posterior(var.ID)
#residual variance
var.res<-posterior_samples(m225_brm)$"sigma"^2
mean(var.res);HPDinterval(as.mcmc(var.res),0.95)#this computes CrI
describe_posterior(var.res)

#repeatability
repeatability<-var.ID/(var.ID+var.res)
describe_posterior(repeatability)
#pair wise differences for nest stage
(group_diff2 <- emmeans(m225_brm, pairwise~Stage*Site ))
(group_diff2 <- emmeans(m225_brm, pairwise~Stage ))

#conditional effect plots
me3<-conditional_effects(m225_brm, "Stage:Site")
plot(me3,theme = theme_classic()) 
me3<-conditional_effects(m225_brm, "Age")
plot(me3,theme = theme_classic()) 



##Extract posterior for site  plot

nn<-m225_brm %>%  spread_draws(b_Intercept, b_SiteWebb) %>%
  head(20000)%>%
  mutate(Webb_mean = (b_Intercept + b_SiteWebb))

data9<- data.frame(nn)
data9<-data9%>%select(b_Intercept, Webb_mean)%>%
  mutate(Georgia = b_Intercept,
         Webb = Webb_mean)%>%
  select(Georgia, Webb)

#convert data to long format for plotting
require(tidyr)
d1<-data9%>% 
  pivot_longer(c("Georgia", "Webb") , 
               names_to= "Site", values_to= "Estimate")

require(dplyr)
d11<- d1 %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(mean = mean(Estimate))%>%
  dplyr::ungroup()


require(dplyr)
d11<- d1 %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(mean = mean(Estimate))%>%
  dplyr::ungroup()

d11$Estimate <-
  d11$Estimate*122+ 133 #back transform

#plot for figure 1
plot_4a<-ggplot(data=d11,aes(y = Site, x = Estimate, fill=Site)) +
  stat_halfeye(point_interval = median_qi,scale=0.9, 
               slab_alpha= 1, point_size=24,stroke=1, interval_alpha=1,
               .width = c(0.67, 0.95),
               interval_size_range = c(5.5, 9.5)) +
  labs(y = "Study Site",
       x = "Distance to edge landcover (m)",
       element_text(face="bold"))+ 
  scale_fill_manual(values = c("#993333",
                               "#999933"))+
  theme(plot.title = element_text(face = "bold", hjust=0.9), plot.subtitle=element_text(hjust=0.9))+
  scale_x_continuous(breaks = c(75,125,175,225))+theme_turkey()+
  ggtitle(label ="G.",
          subtitle = "")+theme(legend.position="none")+theme(legend.position="none")+theme(legend.position="none")+  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 74, face = "bold"),
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
        
        panel.grid.minor = element_blank()) + scale_y_discrete(labels = c('Georgia','South Carolina'))
plot_4a



#stage plot for figure 1
me3<-conditional_effects(m225_brm, "Stage:Site")

dp<-data.frame(me3[["Stage:Site"]])
dp$Site<-as.factor(dp$Site)
dp$estimate__ <-
  dp$estimate__*122 + 133  #back transform
dp$lower__ <-
  dp$lower__*122+ 133  #back transform
dp$upper__ <-
  dp$upper__*122 + 133 #back transform

#plot for figure 1
plot_4<-ggplot(dp, aes(x = Stage, y = estimate__, colour=Site, group=Site)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__), width = .5,
                position=position_dodge(width=0.65),size=4.85) +
  geom_point(position=position_dodge(width=0.65),size=24) +
  scale_y_continuous("Distance to edge landcover (m)", breaks = c(75,125,175,225))+
  theme_classic() + theme_turkey()+
  scale_x_discrete(labels = c('Pre-Hunt','Hunt','Post-Hunt'))+
  theme(legend.position="none")+
  scale_color_manual(values = c("#993333",
                                "#999933"))+
  ggtitle(label ="H.",
          subtitle = "")+theme(legend.position="none")+  theme(legend.title = element_text( size=74, face="bold"), legend.position = "bottom") + 
  theme(axis.title = element_text(size = 74, face = "bold"),
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


plot_4


#this plots a common legend that is used across all panels in figure 1
plot_legend<-ggplot(dp, aes(x = Stage, y = estimate__, colour=Site)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__), width = .5,
                position=position_dodge(width=0.65),size=4.85) +
  geom_point(position=position_dodge(width=0.65),size=24) +
  scale_y_continuous("Distance to open landcover (m)") + theme_classic() + theme_turkey()+
  scale_x_discrete(labels = c('Pre-Hunt','Hunt','Post-Hunt'))+
  theme(legend.position = "bottom", legend.title = element_text(size=84),
        legend.text = element_text(size=84),
        legend.key.size = unit(14.7,"cm"),
        legend.key.width = unit(14, "cm"))+
  scale_color_manual(name= "Population",labels = c("Georgia", "South Carolina") ,values = c("#993333",
                                "#999933")) #change legend text font size


plot_legend
#extract legend
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}
# Apply user-defined function to extract legend
shared_legend <- extract_legend(plot_legend)


# Final plot for interaction
require(gridExtra)
g6<-grid.arrange(arrangeGrob(plot_1, plot_2,plot_3,plot_4, ncol = 4),
             shared_legend, nrow = 2, heights = c(10, 1))



#plots
g5<-grid.arrange(arrangeGrob(plot_1a, plot_2a,plot_3a,plot_4a, ncol = 4))

#final plot for figure 1
grid.arrange(g5, g6)

