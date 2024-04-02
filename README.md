
## The role of human hunters and natural predators in shaping the selection of behavioral types in male wild turkeys






## Project Description
The expression of behavior can vary both among (i.e., animal personality) and within individuals (i.e., plasticity), and investigating causes and consequences of variation has garnered significant attention. Conversely, studies quantifying harvest-induced selection on behavioral traits have received significantly less attention, and work investigating harvest-induced selection and natural selection simultaneously is rare. We studied sources of variation in three movement traits that represented risk-taking and one trait that represented exploration in male Eastern Wild Turkeys (Meleagris gallopavo silvestris). We used data from 109 males in two hunted populations located in Georgia and South Carolina, USA. We assessed how both hunters and natural predators simultaneously influenced the selection of male turkey personalities. We found significant among-individual variation in all movement traits and considerable plasticity in risk-taking and exploration relative to whether hunting was occurring. We observed that predators selected against similar personalities across both populations, whereas hunters selected for different personalities across populations. We also demonstrated that significant harvest-induced selection acts on risk-taking behaviors in both populations, which could render turkeys more difficult to harvest if these traits are indeed heritable. 


## General Project File Structure

```
├── README.md                                                      <- The top-level README including general project descriptions
|
├── Data
│   ├── Male_gps_Webb.csv                                              <- Data from South Carolina ready for modeling
│   ├── Male_gps_GA_sites.csv                                          <- Data from Georgia ready for modeling
├── Final_Figures
│   ├── GA_Survival_Final                                              <- Selection on BTs from hunters/predators in Georgia
│   ├── Webb_Survival_Final                                            <- Selection on BTs from hunters/predators in South Carolina
│   ├── Main_effect_plot_final                                         <- Univariate results for Georgia and South Carolina
└── Final_scripts
    ├── All_Sites_Combined_Univariate_Models_Figure_1_Script.R         <- Univariate models with main effects and repeatabilty estimates reported in table 1 and figure 1
    ├── GA_SC_Survival_BRMS_Final_Script.R                             <- Script contains survival models and results reported in table 2
    ├── GA_Survival_BRMS_Figure_Script.R                               <- Script produces all plots reported in figure 2
    └── SouthCarolina_Survival_BRMS_Figure_Script.R                    <- Script produces all plots reported in figure 3

```

## Data 
  * **Male_gps_Webb**-      data used for all analysis containing South Carolina males
  * **Male_gps_GA_sites**-  data used for all analysis containing Georgia males


## Markdown file
* Click the link below to access markdown file with full workflow and output from each model. 
* Coming soon

  
![turkey](https://github.com/nickgulotta77/GA_Webb_Male_Survival/assets/56907107/f4e2fab7-0fb8-4c02-9c11-a0fe33888a1b)
