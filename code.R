# ~~~~~~~~~~ loading required libraries ~~~~~~~~~~

library(tidyverse)
library(sjlabelled)
library(R2jags)
library(ggridges)
library(ggplot2)
library(eurostat)
library(sf)

# ~~~~~~~~~~ reading SHARE raw data & data preparation for modelling ~~~~~~~~~~

# Loading the dataset with country names and their corresponding codes as defined by SHARE, formatted into a table.
country_ids <- readRDS("country_ids.rds")

# before covid (before 2020)
sharew8_ph <- haven::read_sav("sharew8_rel8-0-0_ALL_datasets_spss/sharew8_rel8-0-0_ph.sav") %>% 
  select(mergeid, country, ph006d1, ph006d5, ph005_, ph006d2, ph006d3, ph006d6) %>% 
  rename(diabetes = ph006d5,
         heart = ph006d1,
         blood_pressure = ph006d2,
         chronic_lung_disease = ph006d6) 

sharew8_cv_r <- haven::read_sav("sharew8_rel8-0-0_ALL_datasets_spss/sharew8_rel8-0-0_cv_r.sav") %>% 
  mutate(wave = 8) %>% 
  select(mergeid, gender, age_int, wave) %>% 
  rename(age = age_int)

sharew8_gv_health <- haven::read_sav("sharew8_rel8-0-0_ALL_datasets_spss/sharew8_rel8-0-0_gv_health.sav") %>% 
  select(mergeid, bmi)

df_before_covid <- left_join(sharew8_ph, sharew8_cv_r, by = "mergeid")
df_before_covid <- left_join(df_before_covid, sharew8_gv_health, by = "mergeid")
bmi <- df_before_covid %>% select(mergeid, bmi) %>% filter(bmi > 0)

# during covid (since 2020)
sharew8ca <- haven::read_sav("sharew8ca_rel8-0-0_ALL_datasets_spss/sharew8ca_rel8-0-0_ca.sav") %>% 
  mutate(wave = 9) %>% 
  mutate(age = 2020 - cadn003_) %>% 
  select(mergeid, 
         country, 
         gender = cadn042_, 
         age, 
         diabetes = cah004_2, 
         heart = cah004_4, 
         covid = cac005_1, 
         blood_pressure = cah004_3, 
         chronic_lung_disease = cah004_5,
         wave) %>% 
  mutate(heart = ifelse(heart == 5, 0 , heart),
         diabetes = ifelse(diabetes == 5, 0, diabetes),
         chronic_lung_disease = ifelse(chronic_lung_disease == 5, 0, chronic_lung_disease),
         blood_pressure = ifelse(blood_pressure == 5, 0, blood_pressure)) %>% 
  mutate(across(country:wave, ~as.numeric(.)))

sharew9ca <- haven::read_sav("sharew9ca_rel8-0-0_ALL_datasets_spss/sharew9ca_rel8-0-0_ca.sav") %>% 
  mutate(wave = 10,
         age = 2021 - cadn003_) %>% 
  select(mergeid, 
         country, 
         gender = cadn042_, 
         age, 
         diabetes = cah004_2, 
         heart = cah004_4, 
         covid = cac105_1, 
         blood_pressure = cah004_3, 
         chronic_lung_disease = cah004_5, 
         wave) %>% 
  mutate(heart = ifelse(heart == 5, 0 , heart),
         diabetes = ifelse(diabetes == 5, 0, diabetes),
         chronic_lung_disease = ifelse(chronic_lung_disease == 5, 0, chronic_lung_disease),
         blood_pressure = ifelse(blood_pressure == 5, 0, blood_pressure)) %>% 
  mutate(across(country:wave, ~as.numeric(.)))


# merging during covid datasets 
df_during_covid <- rbind(sharew8ca, sharew9ca) %>% 
  group_by(mergeid) %>% 
  summarise(across(everything(), ~max(., na.rm = T))) %>% 
  filter_all(all_vars(!is.infinite(.)))


# filtering individuals who haven't had heart attack before covid pandemic
healthy_people <- df_before_covid %>% filter(heart != 1) %>% pull(mergeid)


df_jags <- df_during_covid %>% 
  filter(mergeid %in% healthy_people, age > 49, country != 25) %>%     # Filtering out individuals younger than 50 years old and those from the country 'Israel'
  mutate(covid_nan = ifelse(covid < 0, 1, 0),                           # Set any negative values for covid to 0
         covid = ifelse(covid < 0, 0, covid)) %>% 
  filter(across(country:wave, ~ . >= 0)) %>%
  mutate(country_id = as.numeric(as.factor(country))) %>% 
  left_join(bmi) %>% 
  na.omit()


#  ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ 
# ~~~~~~~~~~ model fitting  ~~~~~~~~~~
#  ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ 

model_code <- "
model
{
 # Likelihood

 for (i in 1:N){
     y[i] ~ dbin(p[i], 1)

          logit(p[i]) <- theta[country[i]] +
                         alpha[country[i]] * gender[i] + 
                         beta_diab[country[i]] * diab[i] + 
                         beta_lung[country[i]] * lung[i] + 
                         beta_pressure[country[i]] * pressure[i] +
                         beta_BMI[country[i]] * BMI[i] + 
                         beta_covid[country[i]] * covid[i]  + 
                         beta_age[country[i]] * age[i]           
                         
 }

 # Priors
 
 mu_covid  ~ dnorm(0.0, 0.01)
 mu_theta  ~ dnorm(0.0, 0.01)
 mu_alpha  ~ dnorm(0.0, 0.01)
 mu_diab  ~ dnorm(0.0, 0.01)
 mu_lung  ~ dnorm(0.0, 0.01)
 mu_pressure  ~ dnorm(0.0, 0.01)
 mu_BMI  ~ dnorm(0.0, 0.01)
 mu_age  ~ dnorm(0.0, 0.01)
 sigma_covid ~ dt(0, 1, 1)T(0,)
 sigma_theta ~ dt(0, 1, 1)T(0,)
 sigma_alpha ~ dt(0, 1, 1)T(0,)
 sigma_diab ~ dt(0, 1, 1)T(0,)
 sigma_lung ~ dt(0, 1, 1)T(0,)
 sigma_pressure ~ dt(0, 1, 1)T(0,)
 sigma_BMI ~ dt(0, 1, 1)T(0,)
 sigma_age ~ dt(0, 1, 1)T(0,)
 
 
 for (i in 1:K){
    beta_covid[i] ~ dnorm(mu_covid, sigma_covid^-2)
    theta[i] ~ dnorm(mu_theta, sigma_theta^-2)
    alpha[i] ~ dnorm(mu_alpha, sigma_alpha^-2)
    beta_diab[i] ~ dnorm(mu_diab, sigma_diab^-2)
    beta_lung[i] ~ dnorm(mu_lung, sigma_lung^-2)
    beta_pressure[i] ~ dnorm(mu_pressure, sigma_pressure^-2)
    beta_BMI[i] ~ dnorm(mu_BMI, sigma_BMI^-2)
    beta_age[i] ~ dnorm(mu_age, sigma_age^-2)
 }
      
         
}
"
# Choose the parameters to watch
model_parameters <- c('alpha', 
                      'theta',
                      'beta_diab', 
                      'beta_lung', 
                      'beta_pressure', 
                      'beta_BMI', 
                      'beta_covid', 
                      'beta_age',
                      'mu_covid',
                      'sigma_covid',
                      'mu_theta',
                      'sigma_theta',
                      'mu_alpha',
                      'sigma_alpha',
                      'mu_diab',
                      'sigma_diab',
                      'mu_lung',
                      'sigma_lung',
                      'mu_pressure',
                      'sigma_pressure',
                      'mu_BMI',
                      'sigma_BMI',
                      'mu_age',
                      'sigma_age')


# Set up the data
model_data <- list(
  y = df_jags$heart,     
  N = nrow(df_jags),
  age = df_jags$age,
  covid = df_jags$covid,
  gender = df_jags$gender-1,
  gender_id = df_jags$gender,
  diab = df_jags$diabetes,
  lung = df_jags$chronic_lung_disease,
  pressure = df_jags$blood_pressure,
  BMI = df_jags$bmi,
  K = max(df_jags$country_id),
  country = df_jags$country_id)



# Run the model
set.seed(123)
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code),
  n.iter = 5000,
  n.burnin = 500)


plot(model_run)
check = model_run$BUGSoutput$summary
#  ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ 
#  ~~~~~~~~~~ results   ~~~~~~~~~~ 
#  ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ 

# Extracting the country codes assigned by our model for each country and merging them with the corresponding codes assigned by SHARE
country_code_converter <- df_jags %>% 
  select(country, country_id) %>% 
  unique() %>% 
  left_join(country_ids)

fit = model_run

# extracting the simulations matrix from the output of the model fit
post_samples = fit$BUGSoutput$sims.matrix

# plotting the effects
# extracting the effect of the variables 
europe_effect_covid = as.data.frame(post_samples) %>% select(contains('mu_covid'))
europe_effect_diab = as.data.frame(post_samples) %>% select(contains('mu_diab'))
europe_effect_lung = as.data.frame(post_samples) %>% select(contains('mu_lung'))
europe_effect_pressure = as.data.frame(post_samples) %>% select(contains('mu_pressure'))
europe_effect_bmi = as.data.frame(post_samples) %>% select(contains('mu_BMI'))

# binding all the effects together
europe_effect <- cbind(europe_effect_covid, 
                       europe_effect_diab,
                       europe_effect_lung,
                       europe_effect_pressure,
                       europe_effect_bmi) %>% 
                mutate(across(everything(.), ~(exp(.)-1)))

# calculating 95 and 80 percentile values 
europe_summary = apply(europe_effect, 2, function(x){quantile(x, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))})
europe_summary_df <- as.data.frame(t(europe_summary)) %>% 
  mutate(var = c("Covid", "Diabetes",  "Chronic lung disease", "Hypertension", "BMI")) #

# plotting (Figure 1 in paper)
europe_summary_df %>% 
  mutate(across(`2.5%`:`97.5%`, ~. * 100)) %>% 
  ggplot(aes(x = `50%`, y = reorder(var, `50%`))) +
  geom_pointrange(aes(xmin = `2.5%`, xmax = `97.5%`), position = position_dodge(width = 0.75), size = 0.5, fatten = 1.25, color = "#85C1E9") + 
  geom_pointrange(aes(xmin = `10%`, xmax = `90%`), position = position_dodge(width = 0.75), size = 2, fatten = 1.25, color = "#2980B9", linewidth = 1.25) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  # xlim(-0.7, 97) +
  geom_point(color = '#A9CCE3', size = 1.5) +
  scale_y_discrete(expand = expansion(mult = 0.0000000000001)) +
  theme_bw() + 
  theme(axis.text.y = element_text(face = "bold", size = 13), 
        axis.text.x = element_text(face = "bold", size = 13),
        axis.title.x = element_text(face = "bold", size = 13)) + # Adjust this size as needed
  labs(x = "Change in odds of CVD (%)", y = NULL)

# to plot the Fig. 2 from the paper
# extracting beta_covid
covid_effects = as.data.frame(post_samples) %>% select(contains('beta_covid'))

# calculating 95 and 80 percentile values 
covid_effects_summary = apply(covid_effects, 2, function(x){quantile(x, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))})

# a function to remove the unwanted characters from the name of columns
extract_first_number <- function(str) {
  # Remove everything before and including the first '['
  str <- gsub(".*\\[", "", str)
  # Remove everything after and including the first comma or closing bracket
  str <- gsub(",.*|\\].*", "", str)
  # Convert to numeric
  as.numeric(str)
}

# transposing the data (changing rows into columns) and extracting the identifiers 
covid_effects_summary_df = as.data.frame(t(covid_effects_summary)) %>% 
  mutate(var_name = rownames(.),
         country_id = extract_first_number(var_name)) %>% 
  select(-var_name)

# remove the row names 
rownames(covid_effects_summary_df) <- NULL

# Join the covid_effects_summary_df dataframe to the dataset containing country names
df_final <- covid_effects_summary_df %>% left_join(country_code_converter)

# Calculating the probability of the COVID effect for each country
countries_probs = as.data.frame(sort(round(apply(covid_effects, 2, function(x){mean(x >= europe_effect_covid$mu_covid)}),2))) %>% 
  rename(probs = `sort(round(apply(covid_effects, 2, function(x) {     mean(x >= europe_effect_covid$mu_covid) }), 2))`) %>% t()


df_final1 <- left_join(df_final, country_prob, by = "country_id") %>% 
  mutate(probs = probs * 100)  %>% 
  select(-country) %>% 
  rename(name = country_name) %>% 
  mutate(Odds = round(`50%` * 100, 2))

df_final$name[12] <- "Czechia"

# to check if all the countries we have match the countries in the shapefile
df_merge <- left_join(df_final1, SHP_28, by = "name")


df_merge %>% 
  st_as_sf() %>% 
  ggplot(aes(fill = probs)) +
  geom_sf(lwd = 0.5) +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "red", midpoint = 50, 
                       limits = c(0, 100)) +
  scale_x_continuous(limits = c(-9, 35)) +
  scale_y_continuous(limits = c(35, 70)) +
  theme_bw() + 
  theme(legend.position = c(0.11,0.8),
        panel.background = element_rect(fill = "grey90", colour="blue"),
        legend.text = element_text(color = "black", size = 10),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12, margin = margin(t = 6)),
        axis.text.y = element_text(size = 12, margin = margin(t = 6))) +
  labs(fill = "Probability (%)")


