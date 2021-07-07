
# Power Analysis


# 1. Housekeeping ---------------------------------------------------------

# Load packages

library(tidyverse)
library(jbmisc)
library(furrr)
library(broom)
library(haven)
library(brms)
library(here)


# Set random seed

set.seed(666)



# 2. Use pilot data to develop priors -------------------------------------

# Load pilot data

dta <- 
  read_sav(
    here(
      "_data",
      "raw-data",
      "BES_Results_210420_CLIENT.sav"
    )
  )


# Select variables

dta <- 
  dta %>% 
  select(
    retro = econGenRetroExp,
    retro_treat = retroText
  )


# Convert outcomes and treatments to factors and remove "Don't know" responses

dta <- 
  dta %>% 
  mutate(
    retro =
      retro %>% 
      as_factor(ordered = T) %>% 
      mark_na(c("Don't know", "Skipped", "Not Asked")),
    retro_treat =
      retro_treat %>% 
      as_factor() %>% 
      fct_relevel(
        "3 months",
        "6 months",
        "12 months",
        "2 years",
        "5 years"
      )
  )


# Fit model to the retrospective data

retro_pilot_model <- 
  brm(
    formula =
      bf(retro ~ 1 + retro_treat) +
      lf(disc ~ 0 + retro_treat, cmc = FALSE),
    family = 
      cumulative(
        link = "probit",
        link_disc = "log"
      ),
    prior = 
      prior(normal(0, 1), class = "Intercept") +
      prior(normal(0, 0.25), class = "b") +
      prior(normal(0, 0.25), class = "b", dpar = "disc"),
    data = dta,
    backend = "cmdstanr",
    seed = 666,
    inits = 0,
    cores = 4,
    chains = 4,
    threads = threading(3)
  )


# Get point estimates from the model

priors <- data.frame(fixef(retro_pilot_model))


# Save power priors to disk

saveRDS(priors, here("_output", "power_priors.rds"))



# 3. Create functions and objects for simulation --------------------------

# Set small sample size to fit model that we can refit to the simulated data

N <- 100


# Simulate treatment

sim_dta <- 
  tibble(
    treatment = 
      sample(
        0:1,
        size = N,
        replace = T
      )
  )


# Simulate outcomes with upper-bound effect size

sim_dta <- 
  sim_dta %>% 
  mutate(
    outcome = 
      rlikert(
        n = N,
        t1 = priors$Estimate[1],
        t2 = priors$Estimate[2],
        t3 = priors$Estimate[3],
        t4 = priors$Estimate[4],
        mean = 0.3*treatment
      )
  )


# We'll then fit a model to these data to avoid recompiling later

fit <- 
  brm(
    formula = bf(outcome ~ treatment),
    family = cumulative(link = "probit"),
    prior = 
      prior(normal(0, 1), class = "Intercept") +
      prior(normal(0, 0.25), class = "b"),
    data = sim_dta,
    backend = "cmdstanr",
    seed = 666,
    inits = 0,
    cores = 2,
    chains = 2
  )


# Create function to simulate data and fit the model to it

sim_and_fit <- 
  function(
    seed, 
    n,
    fx = 0.3
  ){
    
    # Simulate data
    
    dta <- 
      tibble(
        treatment = 
          sample(
            0:1,
            size = n,
            replace = T
          )
      ) %>% 
      mutate(
        outcome = 
          rlikert(
            n = n,
            t1 = priors$Estimate[1],
            t2 = priors$Estimate[2],
            t3 = priors$Estimate[3],
            t4 = priors$Estimate[4],
            mean = {{fx}}*treatment
          )
      )
    
    
    # Print simulation number
    
    print(paste0("Simulation ", seed, "/1000, Effect Size = ", fx))
    
    
    # Fit model
    
    update(fit,
           newdata = dta,
           seed = seed,
           backend = "cmdstanr",
           warmup = 250,
           iter = 750,
           inits = 0,
           cores = 1,
           chains = 1,
           silent = T,
           refresh = 0) %>% 
      fixef() %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "term") %>% 
      rename(
        est = Estimate,
        error = Est.Error,
        lower = `Q2.5`,
        upper = `Q97.5`
      ) %>% 
      filter(term %in% c("treatment")) %>% 
      mutate(effect_size = {{fx}})
    
  }



# 4. Power analysis -------------------------------------------------------

# Specify future plan to run simulations in parallel

plan(multisession)


# Run simulation assuming 7,000 respondents split into 5 random groups with
# an effect size one half of that as in the BES 2019 data

power_sim <-
  tibble(
    effect = rep(seq(0.025, 0.3, by = 0.025), each = 1000),
    seed = rep(1:1000, 12)
  ) %>% 
  mutate(
    tidy = 
      future_map2(
        seed,
        effect,
        sim_and_fit,
        n = 7000/5,
        .options = furrr_options(seed = 666),
        .progress = TRUE
      )
  ) %>% 
  unnest(tidy) %>% 
  mutate(n = 7000/5)


# Save simulated outcomes to disk

saveRDS(power_sim, here("_output", "power_sim.rds"))


# Calculate power

power_sim %>% 
  group_by(effect_size) %>% 
  mutate(check = ifelse(lower > 0, 1, 0)) %>% 
  summarise(power = mean(check))
