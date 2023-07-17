
library(tidyverse)
library(doParallel)
library(rstan)
library(ggpubr)

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)


# function for Kunz' design
STE_ana <- function(n1,n2,s1,r,p1_H0,p1_H1,p12_H0,p12_H1,p2_H0,p2_H1){
  n <- n1+n2
  
  fun.first.sum <- function(x1){choose(n1, x1) * p2_H0^x1 * (1-p2_H0)^(n1-x1)}
  dfgrid <- expand.grid(x1=0:s1)
  first.sum <- sum(mapply(FUN = fun.first.sum, x1=dfgrid$x1))
  
  expr1 <- function(y1,x1,z1){choose(n1,x1) * choose(x1,z1) * choose((n1-x1), y1)}
  expr2 <- function(y1,x1,z1){(p2_H0 - p12_H0)^(x1 - z1) * (p1_H0 - p12_H0)^y1 * p12_H0^z1 * (1-p1_H0-p2_H0+p12_H0)^(n1-x1-y1)}
  expr3 <- function(y1,x1,z1,y2,x2,z2){choose((n-n1),x2) * choose(x2,z2) * choose(((n-n1)-x2), y2)}
  expr4 <- function(y1,x1,z1,y2,x2,z2){(p2_H0 - p12_H0)^(x2 - z2) * (p1_H0 - p12_H0)^y2 * p12_H0^z2 * (1-p1_H0-p2_H0+p12_H0)^((n-n1)-x2-y2)}
  
  
  sum_x1 <- c()
  for(x1 in (s1+1):n1){
    sum_z1 <- c()
    for(z1 in 0:min(x1,r)){
      sum_y1 <- c()
      for(y1 in 0:min(r-z1, n1-x1)){
        
        
        sum_x2 <- c()
        for(x2 in 0:(n-n1)){
          sum_z2 <- c()
          for(z2 in 0:min(x2,r-z1-y1)){ # added -y1
            sum_y2 <- c()
            for(y2 in 0:min((r-y1-z1-z2),((n-n1)-x2))){
              sum_y2[y2+1] <- expr3(y1=y1,x1=x1,z1=z1,x2=x2,z2=z2,y2=y2) * expr4(y1=y1,x1=x1,z1=z1,x2=x2,z2=z2,y2=y2)
            }
            sum_z2[z2+1] <- sum(sum_y2)
          }
          sum_x2[x2+1] <- sum(sum_z2)
        }
        expr34 <- sum(sum_x2)
        
        
        sum_y1[y1+1] <- expr1(y1=y1,x1=x1,z1=z1) * expr2(y1=y1,x1=x1,z1=z1) * expr34
      }
      sum_z1[z1+1] <- sum(sum_y1)
    }
    sum_x1[x1-s1] <- sum(sum_z1)
  }
  second.sum <- sum(sum_x1)
  
  typeI <- 1 - first.sum - second.sum
  
  PET <- first.sum
  ESS <- n1 + (1-PET)*n2
  
  
  # alternative hypothesis
  fun.first.sum <- function(x1){choose(n1, x1) * p2_H1^x1 * (1-p2_H1)^(n1-x1)}
  dfgrid <- expand.grid(x1=0:s1)
  first.sum <- sum(mapply(FUN = fun.first.sum, x1=dfgrid$x1))
  expr1 <- function(y1,x1,z1){choose(n1,x1) * choose(x1,z1) * choose((n1-x1), y1)}
  expr2 <- function(y1,x1,z1){(p2_H1 - p12_H1)^(x1 - z1) * (p1_H1 - p12_H1)^y1 * p12_H1^z1 * (1-p1_H1-p2_H1+p12_H1)^(n1-x1-y1)}
  expr3 <- function(y1,x1,z1,y2,x2,z2){choose((n-n1),x2) * choose(x2,z2) * choose(((n-n1)-x2), y2)}
  expr4 <- function(y1,x1,z1,y2,x2,z2){(p2_H1 - p12_H1)^(x2 - z2) * (p1_H1 - p12_H1)^y2 * p12_H1^z2 * (1-p1_H1-p2_H1+p12_H1)^((n-n1)-x2-y2)}
  sum_x1 <- c()
  for(x1 in (s1+1):n1){
    sum_z1 <- c()
    for(z1 in 0:min(x1,r)){
      sum_y1 <- c()
      for(y1 in 0:min(r-z1, n1-x1)){
        sum_x2 <- c()
        for(x2 in 0:(n-n1)){
          sum_z2 <- c()
          for(z2 in 0:min(x2,r-z1-y1)){ # added -y1
            sum_y2 <- c()
            for(y2 in 0:min((r-y1-z1-z2),((n-n1)-x2))){
              sum_y2[y2+1] <- expr3(y1=y1,x1=x1,z1=z1,x2=x2,z2=z2,y2=y2) * expr4(y1=y1,x1=x1,z1=z1,x2=x2,z2=z2,y2=y2)
            }
            sum_z2[z2+1] <- sum(sum_y2)
          }
          sum_x2[x2+1] <- sum(sum_z2)
        }
        expr34 <- sum(sum_x2)
        sum_y1[y1+1] <- expr1(y1=y1,x1=x1,z1=z1) * expr2(y1=y1,x1=x1,z1=z1) * expr34
      }
      sum_z1[z1+1] <- sum(sum_y1)
    }
    sum_x1[x1-s1] <- sum(sum_z1)
  }
  second.sum <- sum(sum_x1)
  power <- 1 - first.sum - second.sum
  PET_H1 <- first.sum
  
  return(list("Type_I_Error"=typeI, "Power"=power, "PET"=PET, "PET_H1"=PET_H1, "ESS"=ESS))
}


# Stan model for STE-PoS
smodel = "data {
      int<lower=1> K; 
      int<lower=0> N[K];
      int<lower=0> X[K];
      real<lower=0> square_root;
    }
      
      
      parameters {
      real<lower=0> neg_ln_p[K];
      }
      
      transformed parameters {
        real<lower=0,upper=1> p_primary;
        real<lower=0,upper=1> p[K];
        
        for (k in 1:K){
          p[k] = exp(neg_ln_p[k] * (-1));
        }
        p_primary = prod(p);
      }
      
      model {
        for (k in 1:K){
          X[k] ~ binomial(N[k],p[k]);
          neg_ln_p[k] ~ gamma(square_root, 1);
        }
      }
      "

stanmodel <- stan_model(model_code = smodel)

### example run of the model
# N=c(0,0,0,0)
# X=c(0,0,0,0)
# K=length(N)
# square_root=1/K
# 
# data <- list("N"=N, "X"=X, "K"=K, "square_root"=square_root)
# fit <- sampling(stanmodel, data = data, warmup = 2000, iter = 4000, chains = 4, cores = 1, thin = 1, refresh = 0)
# fit



### function for STE-PoS and STE-CP (requires interim data as input)
twostageSTE <- function(pat_data, pat_data_total, timepoints_of_measurement, primary_endpoint, n1, N, r, priors_a, priors_b, STAN=stanmodel, nsamples=2500){
  n2 <- N-n1
 
  # count cases that were survivors at the previous time point AND have an outcome on the current time point
  N_current <- c()
  for(timepoint_index in 1:length(timepoints_of_measurement)){
    N_current[timepoint_index] <- pat_data %>% add_column("PFS_0"=TRUE, .after=0) %>% select(starts_with("PFS_")) %>% filter(.[[timepoint_index]] == TRUE) %>% drop_na(timepoint_index+1) %>% nrow(.)
  }
  
  # count successes for each timepoint
  X_current <- as.numeric(colSums(pat_data, na.rm=TRUE)[grepl("PFS_", names(colSums(pat_data, na.rm=TRUE)))]) 
  
  # pass data into STAN
  data <- list("N"=N_current, "X"=X_current, "K"=length(N_current), "square_root"=1/length(N_current))
  
  fit <- sampling(STAN, data = data, warmup = 1000, iter = 1000+nsamples, chains = 4, cores = 1, thin = 1, refresh = 0, control=list(adapt_delta=0.95, max_treedepth=15))
  rbeta_samples <- data.frame(extract(fit)$p) 

  # there are m (n2 plus patients with short-term endpoint but unknown long-term endpoint) still to evaluate
  # and we need at least r-x11+1 successes among them
  # but the probability differs, depending on which timepoints are known, so determine for each patient first_missing_timepoint
  
  first_missing_timepoints <- pat_data_total %>% rownames_to_column(var = "ID") %>% #filter_at(vars(starts_with("PFS_")), any_vars(!(.==FALSE)))
    pivot_longer(cols=starts_with("PFS_"), names_prefix="PFS_", names_to="timepoint", values_to="PFS") %>%
    group_by(ID) %>% filter(!any(PFS %in% FALSE)) %>%
    filter(is.na(PFS)) %>%
    group_by(ID) %>%
    summarize(first_missing_timepoint = as.numeric(first(timepoint))) %>% 
    filter(first_missing_timepoint <= primary_endpoint)
  
  # open timepoints: patients of second stage all need the first time point
  open_timepoints <- sort(unique(first_missing_timepoints$first_missing_timepoint))
  
  
  ##### PoS
  posterior_prob <- list()
  pred_prim_endpoint <- list()
  for(i in 1:length(open_timepoints)){
    # calculate posterior probability for primary endpoint
    posterior_prob[[i]] <- as.data.frame(rbeta_samples[,match(x=open_timepoints[i], table=timepoints_of_measurement):match(x=primary_endpoint, table=timepoints_of_measurement)]) %>% 
      mutate(product = Reduce(`*`, .)) %>% select(product) # multiplies all columns in this data frame
    
    #posterior_prob[[i]] <- rbeta_samples[,match(x=open_timepoints[i], table=timepoints_of_measurement)]
    
    # sample prediction of number of successes conditional on the number of open patients
    number_pat_timepoint <- first_missing_timepoints %>% filter(first_missing_timepoint==open_timepoints[i]) %>% nrow(.)
    pred_prim_endpoint[[i]] <- rbinom(nsamples, number_pat_timepoint, posterior_prob[[i]] %>% unlist() %>% as.numeric())
  }
  pred_prim_endpoint <- bind_cols(pred_prim_endpoint) %>% rowSums()
  PoS <- sum((pred_prim_endpoint + tail(X_current, 1)) >r)/nsamples
  
  ##### conditional power
  pred_prim_endpoint_cond_power <- list()
  est <- rep(NA,length(open_timepoints))
  for(i in 1:length(open_timepoints)){
    # calculate conditional probability for primary endpoint
    # but this is only possible if we have data observed for this specific time point 
    # otherwise, plug-in H1 (assuming that the probability to survive after each time step until primary endpoint is equal)
    est[i] <- apply(
      data.frame("X_current"=X_current[match(x=open_timepoints[i], table=timepoints_of_measurement):match(x=primary_endpoint, table=timepoints_of_measurement)], 
                 "N_current"=N_current[match(x=open_timepoints[i], table=timepoints_of_measurement):match(x=primary_endpoint, table=timepoints_of_measurement)]) %>% 
        mutate(prob = ifelse(N_current==0, (priors_a/(priors_a+priors_b))^(1/length(timepoints_of_measurement)), X_current/N_current)), 
      2, prod)["prob"]
    
    # sample prediction of number of successes conditional on the number of open patients
    number_pat_timepoint <- first_missing_timepoints %>% filter(first_missing_timepoint==open_timepoints[i]) %>% nrow(.)
    pred_prim_endpoint_cond_power[[i]] <- rbinom(nsamples, number_pat_timepoint, est[i])
  }
  pred_prim_endpoint_cond_power <- bind_cols(pred_prim_endpoint_cond_power) %>% rowSums()
  cond_power <- sum((pred_prim_endpoint_cond_power + tail(X_current, 1)) >r)/nsamples
  
  return(list("PoS"=PoS, "cond_power"=cond_power))
}


### function to build simulation scenarios 
find_weibull <- function(PFS_1_time, PFS_1_probability, PFS_2_time, PFS_2_probability){
  x1=PFS_1_time; x2=PFS_2_time; p1=1-PFS_1_probability; p2=1-PFS_2_probability
  shape <- (log(-log(1-p2)) - log(-log(1-p1))) / (log(x2)-log(x1))
  scale <- x1 / (-log(1-p1))^(1/shape)
  median_survival <- scale*log(2)^(1/shape)
  return(list("shape"=shape, "scale"=scale, "median"=median_survival, "PFS_1"=1-pweibull(q=PFS_1_time,shape=shape,scale=scale), "PFS_2"=1-pweibull(q=PFS_2_time,shape=shape,scale=scale),
              "PFS_1_time"=PFS_1_time, "PFS_2_time"=PFS_2_time))
}



### function to simulate a single two-stage trial
trial_simulator <- function(N, n1, r, short_term_endpoint, primary_endpoint, timepoints_of_measurement, patients_per_month, weibull_params,
                            priors_a=rep(1,length(timepoints_of_measurement)), priors_b=rep(1,length(timepoints_of_measurement))){
  
  time_stamp <- 0; N_current <- 0; pat_data_total <- NULL
  
  priors_a <<- priors_a 
  priors_b <<- priors_b
  
  while(N_current < N){
    patients <- rpois(n=1, lambda=patients_per_month)
    event_time <- rweibull(n=patients, shape=weibull_params$shape, scale=weibull_params$scale)
    if(patients>0){
      pat_data_total <- rbind(pat_data_total, data.frame(event_time, time_stamp))
    }
    N_current <- ifelse(is.numeric(nrow(pat_data_total)), nrow(pat_data_total), 0)
    time_stamp <- time_stamp+1
  }
  
  # remove over-recruited patients
  pat_data_total <- pat_data_total[1:N,]
  
  # add stopping time during interim due to Kunz' design to patients of second stage
  pat_data_total[(n1+1):N,]$time_stamp <- pat_data_total[(n1+1):N,]$time_stamp + short_term_endpoint
  
  # define interime data
  pat_data <- pat_data_total[1:n1,]
  
  
  current_time <- max(pat_data$time_stamp) + short_term_endpoint
  
  pat_data_list <- list()
  for(i in 1:length(timepoints_of_measurement)){
    pat_data_list[[i]] <- pat_data %>% mutate(time = timepoints_of_measurement[i], 
                                              time_in_study = current_time - time_stamp,
                                              PFS = case_when(
                                                (event_time <= timepoints_of_measurement[i]) & (time_in_study >= timepoints_of_measurement[i]) ~ FALSE,
                                                (event_time > timepoints_of_measurement[i]) & (time_in_study >= timepoints_of_measurement[i]) ~ TRUE,
                                                time_in_study < timepoints_of_measurement[i] ~ NA
                                              ))
  }
  pat_data <- bind_rows(pat_data_list) %>% pivot_wider(., names_from=time, names_prefix="PFS_", values_from=PFS)
  
  pat_data_total <- bind_rows(pat_data, pat_data_total[(n1+1):N,]) #%>% print(n=nrow(.))
  
  PoS_cond_power <- twostageSTE(pat_data=pat_data, pat_data_total=pat_data_total, timepoints_of_measurement=timepoints_of_measurement, primary_endpoint=primary_endpoint,
                                r=r, n1=n1, N=N, priors_a=priors_a, priors_b=priors_b)
  PoS <- PoS_cond_power$PoS
  cond_power <- PoS_cond_power$cond_power
  
  # simulate missing patients from first stage and the second stage to be able to calculate power
  successes <- sum(pat_data_total$event_time > primary_endpoint)
  
  # number of short-term survivors at interim 
  successes_STE_interim <- sum(pat_data$event_time > short_term_endpoint)
  
  # number of long-term survivors at interim 
  successes_LTE_interim <- sum(pat_data$event_time > primary_endpoint)
  
  # number of available short-term data at interim
  available_STE_interim <- sum(pat_data$time_in_study > short_term_endpoint)

  # number of available long-term data at interim
  available_LTE_interim <- sum(pat_data$time_in_study > primary_endpoint)
  
  return(list("pat_data"=pat_data, "pat_data_total"=pat_data_total, "PoS"=PoS, "cond_power"=cond_power, "successes"=successes, 
              "successes_STE_interim"=successes_STE_interim, "successes_LTE_interim"=successes_LTE_interim,
              "available_STE_interim"=available_STE_interim, "available_LTE_interim"=available_LTE_interim))
}



######### define simulation scenarios
### Scenario 1
#  n1    n2    s1     r  typeI power PET_H0 PET_H1   ESS
#  18    20    13     6 0.0861 0.951  0.667 0.0282  24.7

### Scenario 2
# n1    n2    s1     r  typeI power PET_H0 PET_H1   ESS
# 17    15    12    10 0.0875 0.955  0.611 0.0221  22.8

### Scenario 3
# n1    n2    s1     r  typeI power PET_H0 PET_H1   ESS
# 16    15    11    20 0.0881 0.955  0.550 0.0170  22.7


patients_per_month <- c(0.5,1,2,4)
#which_timepoints <- list(c(3,12), c(3,6,9,12)) # could simulate different scenarios for the number of time points, but results were too similar to be relevant for the manuscript 
which_timepoints <- list(c(3,6,9,12))
timepoints_of_measurement_index <- 1:length(which_timepoints)


# design parameters

# paste() is required since seq(0.4,0.9,by=0.1) == 0.7 returns FALSE

PFS_1_probability_input <- list(paste(seq(0.3,0.9,by=0.1)), paste(seq(0.4,0.9,by=0.1)), paste(seq(0.7,0.9,by=0.05)))
PFS_2_probability_input <- list(paste(c(0.127, 0.317)),paste(c(0.25,0.5)),paste(c(0.55,0.8)))

# function to find suitable alpha parameter for the beta distribution given a target mean and beta=1
mean_beta <- function(a,b,target_mean) {abs(target_mean-a/(a+b))}
find_alpha_parm <- function(target_mean) {optimize(f=mean_beta,b=1,target_mean=target_mean,interval=c(0,10))$minimum}

scenario1 <- expand.grid("DGP"=1,N=38,n1=18,r=6,s1=13,"PFS_1_probability"=PFS_1_probability_input[[1]], "PFS_2_probability"=PFS_2_probability_input[[1]],"patients_per_month"=patients_per_month, "timepoints_of_measurement_index"=timepoints_of_measurement_index, priors_a=find_alpha_parm(0.317), priors_b=1, stringsAsFactors=FALSE)
scenario2 <- expand.grid("DGP"=2,N=32,n1=17,r=10,s1=12,"PFS_1_probability"=PFS_1_probability_input[[2]], "PFS_2_probability"=PFS_2_probability_input[[2]],"patients_per_month"=patients_per_month, "timepoints_of_measurement_index"=timepoints_of_measurement_index, priors_a=find_alpha_parm(0.5), priors_b=1, stringsAsFactors=FALSE)
scenario3 <- expand.grid("DGP"=3,N=30,n1=18,r=19,s1=13,"PFS_1_probability"=PFS_1_probability_input[[3]], "PFS_2_probability"=PFS_2_probability_input[[3]],"patients_per_month"=patients_per_month, "timepoints_of_measurement_index"=timepoints_of_measurement_index, priors_a=find_alpha_parm(0.8), priors_b=1, stringsAsFactors=FALSE)

scenarios <- rbind(scenario1, scenario2, scenario3) %>% as_tibble() %>% 
  mutate(PFS_2_probability=as.numeric(PFS_2_probability), PFS_1_probability=as.numeric(PFS_1_probability)) %>% 
  filter(PFS_2_probability<PFS_1_probability)

### perform simulations 

# total number of simulations per scenario
nsim <- 1e5

# to save results in between, split up simulations
nsim_per_run <- 1e4

for(run in 1:(nsim/nsim_per_run)){
  ncores <- if((detectCores()-1)>nrow(scenarios)){nrow(scenarios)} else{detectCores()-1}; 
  cl <- makeCluster(ncores); registerDoParallel(cl);
  start.time <- Sys.time()
  PoS <- list()
  PoS <- foreach(scenario=1:nrow(scenarios), .combine=rbind, .packages=c("tidyverse", "rstan")) %dopar% {
    
    current.scenario <- scenarios[scenario,]
    fit_weibull <- find_weibull(PFS_1_time=3, PFS_1_probability=current.scenario$PFS_1_probability, PFS_2_time=12, PFS_2_probability=current.scenario$PFS_2_probability)
    timepoints_of_measurement <- which_timepoints[[current.scenario$timepoints_of_measurement_index]]
    cond_power <- PoS <- rep(NA, nsim_per_run)
    successes <- successes_STE_interim <- successes_LTE_interim <- available_STE_interim <- available_LTE_interim <- rep(NA, nsim_per_run)
    for(sim in 1:nsim_per_run){
      STE_PoS_sim <- trial_simulator(N=current.scenario$N, n1=current.scenario$n1, r=current.scenario$r, 
                                     short_term_endpoint=fit_weibull$PFS_1_time, primary_endpoint=fit_weibull$PFS_2_time, 
                                     timepoints_of_measurement=timepoints_of_measurement,
                                     patients_per_month=current.scenario$patients_per_month, 
                                     weibull_params=list("shape"=fit_weibull$shape, "scale"=fit_weibull$scale), 
                                     priors_a=current.scenario$priors_a, priors_b=current.scenario$priors_b)
      PoS[sim] <- STE_PoS_sim$PoS
      cond_power[sim] <- STE_PoS_sim$cond_power
      successes[sim] <- STE_PoS_sim$successes
      successes_STE_interim[sim] <- STE_PoS_sim$successes_STE_interim
      successes_LTE_interim[sim] <- STE_PoS_sim$successes_LTE_interim
      available_STE_interim[sim] <- STE_PoS_sim$available_STE_interim
      available_LTE_interim[sim] <- STE_PoS_sim$available_LTE_interim
    }  
    
    return(data.frame("scenario"=rep(scenario,nsim_per_run), "PoS"=PoS, "cond_power"=cond_power, "successes"=successes, "successes_STE_interim"=successes_STE_interim, "successes_LTE_interim"=successes_LTE_interim, 
                      "available_STE_interim"=available_STE_interim, "available_LTE_interim"=available_LTE_interim))
  }
  
  Sys.time() - start.time
  stopCluster(cl) # end parallel computing 
  saveRDS(PoS, paste0("sim_results_PoS",run,".rds"))
}


### analyze simulation results

PoS <- rbind(readRDS("sim_results_PoS1.rds"),readRDS("sim_results_PoS2.rds"),readRDS("sim_results_PoS3.rds"),readRDS("sim_results_PoS4.rds"),readRDS("sim_results_PoS5.rds"),
             readRDS("sim_results_PoS6.rds"),readRDS("sim_results_PoS7.rds"),readRDS("sim_results_PoS8.rds"),readRDS("sim_results_PoS9.rds"),readRDS("sim_results_PoS10.rds"),
             readRDS("sim_results_PoS11.rds"),readRDS("sim_results_PoS12.rds"),readRDS("sim_results_PoS13.rds"),readRDS("sim_results_PoS14.rds"),readRDS("sim_results_PoS15.rds"))

scenarios <- rowid_to_column(scenarios, "scenario")

results <- left_join(PoS, scenarios, by="scenario") %>% as_tibble()

results$timepoints_of_measurement_index <- as.factor(results$timepoints_of_measurement_index)
results$PFS_1_probability <- as.numeric(results$PFS_1_probability)

# if goal is to calibrate stopping rule such that PET is the same as with Kunz' design
# PET_kunz <- STE_ana(n1=19, n2=24, s1=13, r=8, p1_H0=0.127, p1_H1=0.317, p12_H0=0.127, p12_H1=0.317, p2_H0=0.7, p2_H1=0.9)$PET
# critical_PoS <- results %>% filter(PFS_1_probability == "0.7" & PFS_2_probability == 0.127 & patients_per_month==1) %>% 
#   summarise(critical_PoS = quantile(PoS, probs=PET_kunz)) %>% as.numeric()
# critical_cond_power <- results %>% filter(PFS_1_probability == "0.7" & PFS_2_probability == 0.127 & patients_per_month==1) %>% 
#   summarise(critical_cond_power = quantile(cond_power, probs=PET_kunz)) %>% as.numeric()


###################################################################################################
# under H0, type I error should be 0.1
cutoffs <- results %>% 
  filter(PFS_1_probability == "0.7" & patients_per_month==0.5) %>% dplyr::select(DGP, PFS_2_probability, successes, r, PoS, cond_power) %>% 
  filter(PFS_2_probability == c(0.127, 0.25, 0.55)[DGP]) %>% 
  mutate(PoS = case_when(successes<=r ~ 0, TRUE~PoS), cond_power = case_when(successes<=r ~ 0, TRUE~cond_power)) %>% 
  group_by(DGP) %>% summarise(PoS_cutoff = quantile(PoS, probs=0.901), CP_cutoff = quantile(cond_power, probs=0.901))

results <- results %>% mutate(critical_PoS = cutoffs$PoS_cutoff[DGP], critical_cond_power = cutoffs$CP_cutoff[DGP]) %>% 
  mutate(rejectH0_STE_PoS = (successes>r & PoS>critical_PoS), ET_STE_PoS = PoS<critical_PoS,
         rejectH0_cond_power = (successes>r & cond_power>critical_cond_power), ET_cond_power = cond_power<critical_cond_power,
         rejectH0_kunz = (successes>r & successes_STE_interim>s1), ET_kunz = successes_STE_interim<=s1)

###################################################################################################




#######################################################################
# plot PoS

plot_PoS <- function(dgp, min_pfs1, method, add_to_title=""){
  if(dgp!=3){
    scenario_labeller <- labeller(
      patients_per_month = c(`0.5` = "Patients/month:\n0.5", `1` = "Patients/month:\n1", `2` = "Patients/month:\n2",
                             `4` = "Patients/month:\n4", `10` = "Patients/month:\n10"),
      PFS_1_probability = c(`0.4` = "3-months survival:\n0.4", `0.5` = "3-months survival:\n0.5", `0.6` = "3-months survival:\n0.6",
                            `0.7` = "3-months survival:\n0.7", `0.8` = "3-months survival:\n0.8", `0.9` = "3-months survival:\n0.9")
    )} else{
      scenario_labeller <- labeller(
        patients_per_month = c(`0.5` = "Patients/month:\n0.5", `1` = "Patients/month:\n1", `2` = "Patients/month:\n2",
                               `4` = "Patients/month:\n4", `10` = "Patients/month:\n10"),
        PFS_1_probability = c(`0.8` = "3-months survival:\n0.8", `0.85` = "3-months survival:\n0.85", `0.9` = "3-months survival:\n0.9")
      )}
  
  if(method=="cond_power"){
    plottitle <- paste0("STE-CP", add_to_title)
    xlabel <- "Conditional Power (CP)"
  }else{
    plottitle <- paste0("STE-PoS", add_to_title)
    xlabel <- "Probability of success (PoS)"
    }
  p <- results %>% filter(DGP==dgp & PFS_1_probability>=min_pfs1) %>% 
    ggplot(.) + geom_histogram(aes(x=.data[[method]], y=100*(..count..)/nsim, fill=as.factor(PFS_2_probability)), alpha=0.5, position="identity", binwidth=0.05) +
    facet_grid(patients_per_month ~ PFS_1_probability, labeller = scenario_labeller) + 
    theme_bw() + scale_fill_discrete(name="", labels=c("H0", "H1")) +
    scale_x_continuous(name=xlabel, breaks=c(0,0.5,1)) + scale_y_continuous("Percentage of simulated trials") +
    theme(legend.position = "right") + ggtitle(plottitle)
  return(p)
}

p1 <- plot_PoS(dgp=1, min_pfs1=0.4, method="cond_power")
p2 <- plot_PoS(dgp=1, min_pfs1=0.4, method="PoS")
ggarrange(p1,p2,ncol=1,common.legend=TRUE,legend="right",labels=c("A","B"))
ggsave("CP_and_PoS.svg", device=svg, width=10, height=10)

p1 <- plot_PoS(dgp=2, min_pfs1=0.6, method="cond_power", add_to_title=" (Scenario 2)")
p2 <- plot_PoS(dgp=2, min_pfs1=0.6, method="PoS", add_to_title=" (Scenario 2)")
p3 <- plot_PoS(dgp=3, min_pfs1=0.85, method="cond_power", add_to_title=" (Scenario 3)")
p4 <- plot_PoS(dgp=3, min_pfs1=0.85, method="PoS", add_to_title=" (Scenario 3)")
ggarrange(p1,p2,p3,p4,nrow=2,ncol=2,common.legend=TRUE,legend="right",labels=c("A","B","C","D"))
ggsave("CP_and_PoS_scenario23.svg", device=svg, width=12, height=10)


# single scenario
results %>% 
  filter(PFS_1_probability == 0.4, patients_per_month==2) %>%
  ggplot(.) + geom_histogram(aes(x=cond_power, y=100*(..count..)/nsim, fill=as.factor(PFS_2_probability)), alpha=0.5, position="identity", binwidth=0.05) +
  theme_bw() + scale_fill_discrete(name="", labels=c("H0", "H1")) +
  scale_x_continuous("Posterior predictive probability of success (PoS)", breaks=c(0,0.5,1)) + scale_y_continuous("Percentage of simulated trials") +
  theme(legend.position = "right") + ggtitle("STE-CP: the conditional power (CP)\nunder the null and alternative hypothesis\nfor pSTE = 0.4 and 2 patients per month.")
ggsave(filename="single_scenario_CP.svg", device=svg, width=4, height=3.5)

###############################################################################################################################################################
# Operating characteristics

results <- results %>% mutate(hypothesis = case_when(DGP==1 & PFS_2_probability==0.127 ~ "Null hypothesis",
                                                     DGP==1 & PFS_2_probability==0.317 ~ "Alternative hypothesis",
                                                     DGP==2 & PFS_2_probability==0.25 ~ "Null hypothesis",
                                                     DGP==2 & PFS_2_probability==0.5 ~ "Alternative hypothesis",
                                                     DGP==3 & PFS_2_probability==0.55 ~ "Null hypothesis",
                                                     DGP==3 & PFS_2_probability==0.8 ~ "Alternative hypothesis"))

scenario_labeller2 <- labeller(
  `patients_per_month` = c(`0.5` = "Patients per month:\n0.5", `1` = "Patients per month:\n1", `2` = "Patients per month:\n2",
                           `4` = "Patients per month:\n4", `10` = "Patients per month:\n10"),
  `hypothesis` = c(`Null hypothesis` = "H0", `Alternative hypothesis` = "H1"))
scenario_labeller3 <- labeller(
  `patients_per_month` = c(`0.5` = "Patients per month:\n0.5", `1` = "Patients per month:\n1", `2` = "Patients per month:\n2",
                           `4` = "Patients per month:\n4", `10` = "Patients per month:\n10"),
  `hypothesis` = c(`Null hypothesis` = "H0", `Alternative hypothesis` = "H1"))

plot_PETs <- function(dgp){
  results %>% filter(DGP==dgp) %>% 
    group_by(scenario, patients_per_month, PFS_1_probability, hypothesis) %>%
    summarise(PET_kunz = sum(ET_kunz) / n(),
              PET_STE_PoS = sum(ET_STE_PoS)/n(),
              PET_cond_power = sum(ET_cond_power)/n()) %>%
    pivot_longer(., cols=starts_with("PET_"), names_to="design", values_to="PET", names_prefix="PET") %>%
    mutate(design = factor(design, levels=c("_kunz","_cond_power","_STE_PoS")),
           hypothesis = factor(hypothesis, levels=c("Null hypothesis","Alternative hypothesis"))) %>%
    ggplot(.) + 
    geom_line(aes(x=PFS_1_probability, y=PET, color=as.factor(design))) + 
    scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.25)) + scale_x_continuous(limits=c(0.3,0.9), breaks = seq(0.3,0.9,by=0.1)) +
    facet_grid(hypothesis~patients_per_month, labeller = scenario_labeller3) +
    scale_color_manual(name="Design", labels=c("Kunz' Design", "STE-CP", "STE-PoS"), values=c("darkgreen","blue","red")) +
    theme_bw() + ggtitle(paste0("Scenario ", dgp)) +
    xlab("3-months survival rate") + ylab("Probability of early termination (PET)")
}
plot_PETs1 <- plot_PETs(1)
plot_PETs2 <- plot_PETs(2)
plot_PETs3 <- plot_PETs(3)

ggarrange(plot_PETs1, plot_PETs2, plot_PETs3, ncol=1, common.legend=FALSE, legend="right")
ggsave("PETs.eps", device="eps", width=14, height=14)


plot_OCs <- function(dgp){
  results %>% filter(DGP==dgp) %>%
    pivot_longer(., cols=starts_with("reject"), values_to="rejectH0", names_to="design", names_prefix="rejectH0") %>%
    group_by(scenario, patients_per_month, PFS_1_probability, hypothesis, design) %>% 
    summarise(rejectH0 = sum(rejectH0)/n()) %>% 
    mutate(y_min = 0, y_max = case_when(hypothesis==0.127 ~ 0.2, hypothesis==0.317 ~ 1),
           design = factor(design, levels=c("_kunz","_cond_power","_STE_PoS")),
           hypothesis = factor(hypothesis, levels=c("Null hypothesis","Alternative hypothesis"))) %>%
    ggplot(.) + 
    geom_line(aes(x=PFS_1_probability, y=rejectH0, color=as.factor(design))) + 
    geom_hline(data=data.frame(hypothesis=factor(c("Null hypothesis","Alternative hypothesis"), levels=c("Null hypothesis","Alternative hypothesis")),hline=c(0.1,0.95)), aes(yintercept=hline), linetype="dashed") + 
    facet_grid(hypothesis~patients_per_month, labeller = scenario_labeller2, scales="free_y") +
    geom_blank(data=data.frame(hypothesis=factor(c("Null hypothesis","Alternative hypothesis"), levels=c("Null hypothesis","Alternative hypothesis")), max.height=c(0.2,1)), aes(y=max.height)) + # ensures common scales for same facet levels
    geom_blank(data=data.frame(hypothesis=factor(c("Null hypothesis","Alternative hypothesis"), levels=c("Null hypothesis","Alternative hypothesis")), min.height=c(0,0)), aes(y=min.height)) + # ensures common scales for same facet levels
    scale_x_continuous(limits=c(0.3,0.9), breaks = seq(0.3,0.9,by=0.1)) + 
    scale_color_manual(name="Design", labels=c("Kunz' Design", "STE-CP", "STE-PoS"), values=c("darkgreen","blue","red")) +
    theme_bw() + ggtitle(paste0("Scenario ", dgp)) +
    xlab("3-months survival rate") + ylab("Probability of rejecting H0")
}
plot_OCs1 <- plot_OCs(1)
plot_OCs2 <- plot_OCs(2)
plot_OCs3 <- plot_OCs(3)

ggarrange(plot_OCs1, plot_OCs2, plot_OCs3, ncol=1, common.legend=FALSE, legend="right")
ggsave("OCs.eps", device="eps", width=14, height=14)
