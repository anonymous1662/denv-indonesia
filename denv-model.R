#### Definition of discrete time stochastic compartmental model
### Structure of script
## Model parameters
## Core equations
## Compute probabilities of transition
## Draw numbers moving between compartments
## Ageing & mortality
## Burden tracking
## Initial states and dimensions

### Compartments S, I and R are stratified by 
## i for the age group
## j for region
## Compartments I and R are also stratified by 
## k for number of past infections
# - k = 1 primary infection
# - k = 2 secondary infection
# - k = 3 tertiary infection
# - k = 4 quaternary infection
# Note that we don't keep track of which serotype is infecting

#### User-input model parameters
## Time-keeping
initial(sim_year) <- 1
# time/365 converts days to year, +1 to start index at year 1 (not 0), floor() rounds down to integer
update(sim_year) <- floor((time)/365 + 1) # track year of simulation, assuming all years are 365 days long. Note this updates in next time step so sim_year == 2 happens when time = 366.
scenario_years <- parameter()

## Demog & serotype parameters
n_age <- parameter()
n_regions <- parameter()
n_serotypes <- parameter()
demog <- parameter()
amp_seas <- parameter(0.2) # Amplitude of seasonality
phase_seas <- parameter(1.56) # Parameter to phase seasonality
beta <- parameter()
gamma <- parameter(0.2) # so mean infectious period is 5 days
seeding_times <- parameter()
seeding_ages <- parameter()
seeding_value <- parameter()

## Burden user inputs
symp_prop <- parameter()
hosp_prop <- parameter()
death_prop <- parameter()

## Initial states
N_init <- parameter()

## Vaccination parameters
# Campaign
vacc_on <- parameter()
vacc_start <- parameter()
catch_up_end <- parameter()
vacc_age <- parameter()
vacc_regions <- parameter()
catch_up_age <- parameter()
vacc_coverage <- parameter()

# Efficacy
ve_inf <- parameter() # efficacy against infection
ve_symp <- parameter() # efficacy against symptomatic dengue
ve_hosp <- parameter() # efficacy against hospitalisation
ve_death <- parameter() # efficacy against death
ve_wane <- parameter() # duration of vaccine protection

# Calculate conditional vaccine efficacies
# ve_dc is the vaccine efficacy against symptomatic disease given infection has occurred
ve_dc[] <- (ve_symp[i] - ve_inf[i]) / (1 - ve_inf[i])
# ve_sd is the vaccine efficacy against hospitalisation (severe dengue) given infection has occurred
ve_sd[] <- (ve_hosp[i] - ve_inf[i]) / (1 - ve_inf[i]) 

## Wolbachia parameters
wol_on <- parameter()
wol_start <- parameter()
wol_reduction <- parameter()
wol_regions <- parameter()

## Track population size by age group
print("sim_year: {sim_year}", when = time %% 365 == 0)

#### Core equations ####

# %% is remainder
# all individuals age on the same day each year
# vaccination occurs the day after aging
ageing_day <- (time %% 365 == 0)
vaccination_day <- (time %% 365 == 1)
day_of_year <- time %%365 + 1

# Updates occur daily
# Update susceptible compartment
# If it is an ageing day, S[1,]= births that day
# Else = births + old S - (s to i) - vaccinated + waned vaccinated, shifted to adjust demographic
update(S[1,]) <- floor((if(ageing_day) births[day_of_year, j] else 
  births[day_of_year,j] + S[i,j] - n_SI[i,j] - vacc_S[i,j] + n_Vwn[i,j,1]) * shift[i,j])

# Update susceptible for all other ages
update(S[2:n_age,]) <- floor((if(ageing_day) S[i-1,j] - n_SI[i-1,j] - vacc_S[i-1,j] + n_Vwn[i-1,j,1] else
      S[i,j] - n_SI[i,j] - vacc_S[i,j] + n_Vwn[i,j,1]) * shift[i,j])

## Primary infection - individuals enter I[i,1] from S and V[1]
# Changed from if (i==1) to if (ageing_day)
update(I[1,,1]) <- floor((if(ageing_day) 0 else # cannot age into the youngest group with an existing infection
    I[i,j,k] + n_SI[i,j] + n_VI[i,j,k] - n_IC[i,j,k]) * shift[i,j]) 

update(I[2:n_age,,1]) <- floor((
  if(ageing_day) I[i-1,j,k] + n_SI[i-1,j] + n_VI[i-1,j,k] - n_IC[i-1,j,k] else
    I[i,j,k] + n_SI[i,j] + n_VI[i,j,k] - n_IC[i,j,k]) * shift[i,j])

# Secondary / Tertiary / Quaternary infection 
## individuals enter[i, 2:4] from R[1:3] and from V[2:4]
update(I[1,,2:n_serotypes]) <- floor((if(ageing_day) 0 else
    I[i,j,k] + n_RI[i,j,k-1] + n_VI[i,j,k] - n_IC[i,j,k]) * shift[i,j])

update(I[2:n_age,,2:n_serotypes]) <- floor((if(ageing_day) I[i-1,j,k] + n_RI[i-1,j,k-1] + n_VI[i-1,j,k]- n_IC[i-1,j,k] else
    I[i,j,k] + n_RI[i,j,k-1] + n_VI[i,j,k] - n_IC[i,j,k]) * shift[i,j])

## Individuals enter C[i,1] from I[i,1]
update(C[1,,]) <- floor((if(ageing_day) 0 else
    C[i,j,k] + n_IC[i,j,k] - n_CR[i,j,k]) * shift[i,j])

update(C[2:n_age,,]) <- floor((if(ageing_day && i > 1) C[i-1,j,k] + n_IC[i-1,j,k] - n_CR[i-1,j,k] else
    C[i,j,k] + n_IC[i,j,k] - n_CR[i,j,k]) * shift[i,j])

## Individuals enter R[i,1] from C[i,1]
update(R[1,,]) <- floor((if(ageing_day) 0 else
    R[i,j,k] + n_CR[i,j,k] - n_RI[i,j,k] - vacc_R[i,j,k] + n_Vwn[i,j,k+1]) * shift[i,j])

update(R[2:n_age,,]) <- floor((if(ageing_day && i > 1) R[i-1,j,k] + n_CR[i-1,j,k] - n_RI[i-1,j,k] - vacc_R[i-1,j,k] + n_Vwn[i-1,j,k+1]  else
    R[i,j,k] + n_CR[i,j,k] - n_RI[i,j,k] - vacc_R[i,j,k] + n_Vwn[i,j,k+1]) * shift[i,j])

# Vaccine compartments
update(V[1,,1]) <- floor((if(ageing_day) 0 else
    V[i,j,1] + vacc_S[i,j] - n_VI[i,j,k] - n_Vwn[i,j,k]) * shift[i,j]) # vaccinated while seronegative

update(V[2:n_age,,1]) <- floor((if(ageing_day && i > 1) V[i-1,j,1] + vacc_S[i-1,j] - n_VI[i-1,j,k] - n_Vwn[i-1,j,k] else
    V[i,j,1] + vacc_S[i,j] - n_VI[i,j,k] - n_Vwn[i,j,k]) * shift[i,j]) # vaccinated while seronegative

update(V[1,,2:5]) <- floor((if(ageing_day && i == 1) 0 else
    V[i,j,k] + vacc_R[i,j,k-1] - n_VI[i,j,k] - n_Vwn[i,j,k]) * shift[i,j]) # vaccinated post primary infection

update(V[2:n_age,,2:5]) <- floor((if(ageing_day && i > 1) V[i-1,j,k] + vacc_R[i-1,j,k-1] - n_VI[i-1,j,k] - n_Vwn[i-1,j,k] else
    V[i,j,k] + vacc_R[i,j,k-1] - n_VI[i,j,k] - n_Vwn[i,j,k]) * shift[i,j]) # vaccinated post primary infection

update(N[,]) <- if(ageing_day) demog[i,j, sim_year + 1] else N[i,j]

#### Calculate individual probabilities of transition ####

## Beta seasonality + Wolbachia (acts on beta) 
# wol_reduction is a fractional reduction in transmission
beta_adj[] <- if(wol_on == 1 && wol_regions[i] == 1 && time >= wol_start) beta[i] * (1-wol_reduction) else beta[i]
beta_t[] <- beta_adj[i] * (1 + amp_seas * cos(2 * 3.14159 * time/365 - phase_seas)) # add seasonality to beta

## Force of infection
# calculate the FOI for each region daily
lambda[] <- if(sum(N[,i]) == 0) 0 else beta_t[i] * sum(I[,i,])/sum(N[,i]) # overall force of infection (all serotypes)
update(local_foi_region[]) <- lambda[i]
update(annual_foi_region[]) <- if(ageing_day) lambda[i] else lambda[i] + annual_foi_region[i]

## Probabilities of transition
# Note dt is set as 1 in run-model.R
p_inf[,1] <- 1 - exp(-lambda[i] * dt) # Probability of transition S to I
p_inf[,2] <- 1 - exp(-lambda[i] * 3/4 * dt) # Probability of transition R1 to I2 
p_inf[,3] <- 1 - exp(-lambda[i] * 2/4 * dt) # Probability of transition R2 to I3 
p_inf[,4] <- 1 - exp(-lambda[i] * 1/4 * dt) # Probability of transition R3 to I4 
p_inf[,5] <- 0 # need this for overall probability of leaving V5 below

p_inf[,] <- if(p_inf[i,j] <1) p_inf[i,j] else 1 # no infection probability exceeds 1

p_IC <- 1 - exp(-gamma * dt) # I to C, gamma is rate of leaving I
p_CR <- 1 - exp(-nu * dt) # C to R
nu <- dt/365 # duration of cross protection of a year (can make this user defined)
p_vw <- 1 - exp(-1 * dt/(ve_wane*365)) # probability of waning from vaccine compartment

## Rt

update(rt_region[]) <- if(sum(N[,i]) > 0) (beta_t[i]/gamma)*((sum(S[,i]) + sum(R[,i,1])*3/4 + sum(R[,i,2])*2/4 + sum(R[,i,3])*1/4 + sum(V[,i,1])*(1-ve_inf[1]) + sum(V[,i,2])*(1-ve_inf[2])*3/4 + sum(V[,i,3])*(1-ve_inf[2])*2/4))/sum(N[,i]) else 0
update(r0_region[]) <- (beta_t[i]/gamma)

#### Draw number moving between compartments ####
# problem encountered where p exceeds 1 for some binomial draws - need to identify where

n_SI[, ] <- Binomial(S[i, j], p_inf[j, 1]) # draw number moving from S to I
# seed additional infections only in chosen age group on chosen day each year
# and as long as the number of people infected is not more than the number in S
n_SI[,] <- if(i == seeding_ages[j, sim_year] && day_of_year == seeding_times[j, sim_year] && seeding_value[j] <= (S[i,j] - n_SI[i,j])) n_SI[i,j] + seeding_value[j] else n_SI[i,j] # add seeding to S compartment
n_IC[, , ] <- Binomial(I[i, j, k], p_IC)
n_CR[, , ] <- Binomial(C[i, j, k], p_CR)
n_RI[, , 1:(n_serotypes - 1)] <- Binomial(R[i, j, k], p_inf[j, k + 1]) # reinfections
n_RI[,,n_serotypes] <- 0 # maximum number of infections is 4, only leave R4 compartment through death

### Rewritten vaccine waning code - sequential vaccine waning then breakthrough infection
# Draw number waning first
n_Vwn[ , , ] <- Binomial(V[i, j, k], p_vw)

# Calculate remaining after waned cases leave
V_remaining[ , , ] <- V[i, j, k] - n_Vwn[i, j, k]

# Draw breakthrough infections from those in V
n_VI[ , ,1]   <- Binomial(V_remaining[i, j, 1],   p_inf[j, 1] * (1 - ve_inf[1]))
n_VI[ , ,2:3] <- Binomial(V_remaining[i, j, k],   p_inf[j, k] * (1 - ve_inf[2]))
n_VI[ , ,4:5] <- 0  # no reinfection from vaccine compartment level 4 and 5


#### Ageing & demography ####

births_per_day[] <- floor(demog[1,i,sim_year] / 365) ## spread births out throughout the year
remainder[] <- demog[1,i,sim_year] %% 365
births[1:365, ] <- if(sim_year != 1 && i <= remainder[j]) births_per_day[j] + 1 else if(sim_year!= 1 && i > remainder[j]) births_per_day[j] else 0

## each year adjust each compartment to match demography data
sum_comp[1:n_age,] <- S[i,j] + sum(I[i,j,]) + sum(C[i,j,]) + sum(R[i,j,]) + sum(V[i,j,])
shift[1,] <- 1 # for youngest age group all enter as susceptible
shift[2:n_age,] <- if(ageing_day && sum_comp[i-1,j] > 0) demog[i, j, sim_year+1]/sum_comp[i-1,j] else 1.0 # adjust population sizes so they match demographic data from upcoming year

#### Interventions ####

## Vaccination (occurs as part of ageing process)
vacc_prop[,] <- if(vacc_on == 1 && vacc_age == i && vacc_regions[j] == 1 && vaccination_day && time >= vacc_start) vacc_coverage else 0 # only vaccinate once a year in line with ageing
vacc_prop[,] <- if(vacc_on == 1 && catch_up_age[i] == 1 && vacc_regions[j] == 1 && vaccination_day && time >= vacc_start && time <= catch_up_end) vacc_coverage else vacc_prop[i,j] # catch-up campaign targets wider age group

vacc_S[,] <-  (S[i,j] - n_SI[i,j])*vacc_prop[i,j] # of those not already leaving compartment, how many are vaccinated
vacc_R[,,] <- (R[i,j,k] - n_RI[i,j,k])*vacc_prop[i,j]

#### Track burden ####

### Infections
inf_i[,] <- n_SI[i,j] + sum(n_RI[i,j,]) + sum(n_VI[i,j,])
update(inf_year[,]) <- if(ageing_day) inf_i[i,j] else inf_year[i,j] + inf_i[i,j]
update(inf_cumulative[,]) <- inf_i[i,j] + inf_cumulative[i,j]

### Cumulative vaccinations
update(vacc_i[,]) <- vacc_S[i,j] + sum(vacc_R[i,j,])
update(vacc_cumulative[,]) <- vacc_S[i,j] + sum(vacc_R[i,j,]) + vacc_cumulative[i,j]
update(vacc_cumulative_seronegative[,]) <- vacc_S[i,j] + vacc_cumulative_seronegative[i,j]
update(vacc_cumulative_seropositive[,]) <-  sum(vacc_R[i,j,]) + vacc_cumulative_seropositive[i,j]

### Cases, hospitalisations and deaths
## Symptomatic cases, hospitalisations and deaths in whole population by age - daily incidence, yearly total and cumulative total 

non_vacc_cases_i[,] <- (n_SI[i,j]*symp_prop[1]) + (n_RI[i,j,1]*symp_prop[2]) + (n_RI[i,j,2]*symp_prop[3]) + (n_RI[i,j,3]*symp_prop[4])
non_vacc_hosp_i[,] <- (n_SI[i,j]*symp_prop[1]*hosp_prop[1]) + (n_RI[i,j,1]*symp_prop[2]*hosp_prop[2]) + (n_RI[i,j,2]*symp_prop[3]*hosp_prop[3]) + (n_RI[i,j,3]*symp_prop[4]*hosp_prop[4])
non_vacc_deaths_i[,] <- (n_SI[i,j]*symp_prop[1]*hosp_prop[1]*death_prop[i]) + (n_RI[i,j,1]*symp_prop[2]*hosp_prop[2]*death_prop[i]) + (n_RI[i,j,2]*symp_prop[3]*hosp_prop[3]*death_prop[i]) + (n_RI[i,j,3]*symp_prop[4]*hosp_prop[4]*death_prop[i])

vacc_cases_i[,] <- n_VI[i,j,1]*symp_prop[1]*(1-ve_dc[1]) + n_VI[i,j,2]*symp_prop[2]*(1-ve_dc[2]) + n_VI[i,j,3]*symp_prop[3]*(1-ve_dc[2])
vacc_hosp_i[,]<- n_VI[i,j,1]*symp_prop[1]*hosp_prop[1]*(1-ve_sd[1]) + n_VI[i,j,2]*symp_prop[2]*hosp_prop[2]*(1-ve_sd[2]) + n_VI[i,j,3]*symp_prop[3]*hosp_prop[3]*(1-ve_sd[2])
vacc_deaths_i[,] <- n_VI[i,j,1]*symp_prop[1]*hosp_prop[1]*death_prop[i]*(1-ve_death[1]) + n_VI[i,j,2]*symp_prop[2]*hosp_prop[2]*death_prop[i]*(1-ve_death[2]) + n_VI[i,j,3]*symp_prop[3]*hosp_prop[3]*death_prop[i]*(1-ve_death[2])

update(cases_i[,]) <- non_vacc_cases_i[i,j] + vacc_cases_i[i,j]
update(hosp_i[,]) <- non_vacc_hosp_i[i,j] + vacc_hosp_i[i,j]
update(deaths_i[,]) <- non_vacc_deaths_i[i,j] + vacc_deaths_i[i,j]

update(cases_year[,]) <- if(ageing_day) non_vacc_cases_i[i,j] + vacc_cases_i[i,j] else cases_year[i,j] + non_vacc_cases_i[i,j] + vacc_cases_i[i,j]
update(hosp_year[,]) <- if(ageing_day) non_vacc_hosp_i[i,j] + vacc_hosp_i[i,j] else hosp_year[i,j] + non_vacc_hosp_i[i,j] + vacc_hosp_i[i,j]
update(deaths_year[,]) <- if(ageing_day) non_vacc_deaths_i[i,j] + vacc_deaths_i[i,j] else deaths_year[i,j] + non_vacc_deaths_i[i,j] + vacc_deaths_i[i,j]

update(cases_cumulative[,]) <- cases_cumulative[i,j] + non_vacc_cases_i[i,j] + vacc_cases_i[i,j]
update(hosp_cumulative[,]) <- hosp_cumulative[i,j] + non_vacc_hosp_i[i,j] + vacc_hosp_i[i,j]
update(deaths_cumulative[,]) <- deaths_cumulative[i,j] + non_vacc_deaths_i[i,j] + vacc_deaths_i[i,j]

## Symptomatic cases, hospitalisations and deaths in vaccinated individuals by age - daily incidence, yearly total and cumulative total

update(vacc_cases_year[,]) <- if(ageing_day) vacc_cases_i[i,j] else vacc_cases_year[i,j] + vacc_cases_i[i,j]
update(vacc_hosp_year[,]) <- if(ageing_day) vacc_hosp_i[i,j] else vacc_hosp_year[i,j] + vacc_hosp_i[i,j]
update(vacc_deaths_year[,]) <- if(ageing_day) vacc_deaths_i[i,j] else vacc_deaths_year[i,j] + vacc_deaths_i[i,j]

update(vacc_cases_cumulative[,]) <- vacc_cases_cumulative[i,j] + vacc_cases_i[i,j]
update(vacc_hosp_cumulative[,]) <- vacc_hosp_cumulative[i,j] + vacc_hosp_i[i,j]
update(vacc_deaths_cumulative[,]) <- vacc_deaths_cumulative[i,j] + vacc_deaths_i[i,j]



##### ADDED: Tracking the burden for seronegative and seropositive vaccine recipients separately #####
# Does not change model, uses values already calculated
### Burden in seronegative vaccine recipients
vacc_seroneg_cases_i[,] <- floor(n_VI[i,j,1]*symp_prop[1]*(1-ve_dc[1]))
vacc_seroneg_hosp_i[,] <- floor(n_VI[i,j,1]*symp_prop[1]*hosp_prop[1]*(1-ve_sd[1]))
vacc_seroneg_deaths_i[,] <- floor(n_VI[i,j,1]*symp_prop[1]*hosp_prop[1]*death_prop[i]*(1-ve_death[1]))

update(vacc_seroneg_cases_year[,]) <- if(ageing_day) vacc_seroneg_cases_i[i,j] else vacc_seroneg_cases_year[i,j] + vacc_seroneg_cases_i[i,j]
update(vacc_seroneg_hosp_year[,]) <- if(ageing_day) vacc_seroneg_hosp_i[i,j] else vacc_seroneg_hosp_year[i,j] + vacc_seroneg_hosp_i[i,j]
update(vacc_seroneg_deaths_year[,]) <- if(ageing_day) vacc_seroneg_deaths_i[i,j] else vacc_seroneg_deaths_year[i,j] + vacc_seroneg_deaths_i[i,j]

update(vacc_seroneg_cases_cumulative[,]) <- vacc_seroneg_cases_cumulative[i,j] + vacc_seroneg_cases_i[i,j]
update(vacc_seroneg_hosp_cumulative[,]) <- vacc_seroneg_hosp_cumulative[i,j] + vacc_seroneg_hosp_i[i,j]
update(vacc_seroneg_deaths_cumulative[,]) <- vacc_seroneg_deaths_cumulative[i,j] + vacc_seroneg_deaths_i[i,j]

### Burden in seropositive vaccine recipients
vacc_seropos_cases_i[,] <- floor(n_VI[i,j,2]*symp_prop[2]*(1-ve_dc[2]) + n_VI[i,j,3]*symp_prop[3]*(1-ve_dc[2]))
vacc_seropos_hosp_i[,] <- floor(n_VI[i,j,2]*symp_prop[2]*hosp_prop[2]*(1-ve_sd[2]) + n_VI[i,j,3]*symp_prop[3]*hosp_prop[3]*(1-ve_sd[2]))
vacc_seropos_deaths_i[,] <- floor(n_VI[i,j,2]*symp_prop[2]*hosp_prop[2]*death_prop[i]*(1-ve_death[2]) + n_VI[i,j,3]*symp_prop[3]*hosp_prop[3]*death_prop[i]*(1-ve_death[2]))

update(vacc_seropos_cases_year[,]) <- if(ageing_day) vacc_seropos_cases_i[i,j] else vacc_seropos_cases_year[i,j] + vacc_seropos_cases_i[i,j]
update(vacc_seropos_hosp_year[,]) <- if(ageing_day) vacc_seropos_hosp_i[i,j] else vacc_seropos_hosp_year[i,j] + vacc_seropos_hosp_i[i,j]
update(vacc_seropos_deaths_year[,]) <- if(ageing_day) vacc_seropos_deaths_i[i,j] else vacc_seropos_deaths_year[i,j] + vacc_seropos_deaths_i[i,j]

update(vacc_seropos_cases_cumulative[,]) <- vacc_seropos_cases_cumulative[i,j] + vacc_seropos_cases_i[i,j]
update(vacc_seropos_hosp_cumulative[,]) <- vacc_seropos_hosp_cumulative[i,j] + vacc_seropos_hosp_i[i,j]
update(vacc_seropos_deaths_cumulative[,]) <- vacc_seropos_deaths_cumulative[i,j] + vacc_seropos_deaths_i[i,j]




#### Initial states & dimensions ####

initial(S[,]) <- if(i == 11) N_init[i,j] - 10 else N_init[i,j]
initial(I[,,]) <- if(i == 11 && k == 1) 10 else 0 
initial(C[,,]) <- 0 
initial(R[,,]) <- 0
initial(N[,]) <- N_init[i,j]
initial(V[,,]) <- 0 
initial(vacc_i[,]) <- 0
initial(vacc_cumulative[,]) <-0
initial(vacc_cumulative_seropositive[,]) <- 0
initial(vacc_cumulative_seronegative[,]) <- 0
initial(inf_year[,]) <- 0
initial(inf_cumulative[,]) <- 0
initial(cases_i[,]) <- 0
initial(hosp_i[,]) <- 0
initial(deaths_i[,]) <- 0
initial(cases_year[,]) <- 0
initial(hosp_year[,]) <- 0
initial(deaths_year[,]) <- 0
initial(cases_cumulative[,]) <- 0
initial(hosp_cumulative[,]) <- 0
initial(deaths_cumulative[,]) <- 0
initial(vacc_cases_year[,]) <- 0
initial(vacc_hosp_year[,]) <- 0
initial(vacc_cases_cumulative[,]) <- 0
initial(vacc_deaths_year[,]) <- 0
initial(vacc_hosp_cumulative[,]) <- 0
initial(vacc_deaths_cumulative[,]) <- 0
initial(local_foi_region[]) <- 0
initial(annual_foi_region[]) <- 0
initial(rt_region[]) <- 0
initial(r0_region[]) <- 0

dim(S) <- c(n_age, n_regions)
dim(I) <- c(n_age, n_regions, n_serotypes)
dim(C) <- c(n_age, n_regions, n_serotypes)
dim(R) <- c(n_age, n_regions, n_serotypes)
dim(V) <- c(n_age, n_regions, 5)
dim(N) <- c(n_age, n_regions)
dim(sum_comp) <- c(n_age, n_regions)
dim(shift) <- c(n_age, n_regions)
dim(n_SI) <- c(n_age, n_regions)
dim(n_IC) <- c(n_age, n_regions, n_serotypes)
dim(n_CR) <- c(n_age, n_regions, n_serotypes)
dim(n_RI) <- c(n_age, n_regions, n_serotypes)
dim(n_VI) <- c(n_age, n_regions, 5)
dim(beta) <- n_regions
dim(beta_adj) <- n_regions
dim(beta_t) <- n_regions
dim(p_inf) <- c(n_regions, 5)
dim(lambda) <- n_regions
dim(seeding_times) <- c(n_regions, scenario_years)
dim(seeding_ages) <- c(n_regions, scenario_years)
dim(N_init) <- c(n_age, n_regions)
dim(demog) <- c(n_age, n_regions, scenario_years)
dim(births) <- c(365, n_regions)
dim(births_per_day) <- n_regions
dim(vacc_prop) <- c(n_age, n_regions)
dim(ve_inf) <- 2
dim(ve_symp) <- 2
dim(ve_hosp) <- 2
dim(ve_death) <- 2
dim(ve_dc) <- 2
dim(ve_sd) <- 2
dim(vacc_S) <- c(n_age, n_regions)
dim(vacc_R) <- c(n_age, n_regions, n_serotypes)
dim(V_remaining) <- c(n_age, n_regions, 5)
dim(n_Vwn) <- c(n_age, n_regions, 5)
dim(remainder) <- n_regions
dim(wol_regions) <- c(n_regions)
dim(vacc_regions) <- c(n_regions)
dim(catch_up_age) <- c(n_age)
dim(symp_prop) <- 4
dim(hosp_prop) <- 4
dim(death_prop) <- n_age
dim(vacc_cumulative) <- c(n_age, n_regions)
dim(vacc_cumulative_seronegative) <- c(n_age, n_regions)
dim(vacc_cumulative_seropositive) <- c(n_age, n_regions)
dim(inf_i) <- c(n_age, n_regions)
dim(inf_year) <- c(n_age, n_regions)
dim(inf_cumulative) <- c(n_age, n_regions)
dim(non_vacc_cases_i) <- c(n_age, n_regions)
dim(non_vacc_hosp_i) <- c(n_age, n_regions)
dim(non_vacc_deaths_i) <- c(n_age, n_regions)
dim(cases_i) <- c(n_age, n_regions)
dim(hosp_i) <- c(n_age, n_regions)
dim(deaths_i) <- c(n_age, n_regions)
dim(cases_year) <- c(n_age, n_regions)
dim(hosp_year) <- c(n_age, n_regions)
dim(deaths_year) <- c(n_age, n_regions)
dim(cases_cumulative) <- c(n_age, n_regions)
dim(hosp_cumulative) <- c(n_age, n_regions)
dim(deaths_cumulative) <- c(n_age, n_regions)
dim(vacc_i) <- c(n_age, n_regions)
dim(vacc_cases_i) <- c(n_age, n_regions)
dim(vacc_hosp_i) <- c(n_age, n_regions)
dim(vacc_deaths_i) <- c(n_age, n_regions)
dim(vacc_cases_year) <- c(n_age, n_regions)
dim(vacc_hosp_year) <-c(n_age, n_regions)
dim(vacc_deaths_year) <- c(n_age, n_regions)
dim(vacc_cases_cumulative) <- c(n_age, n_regions)
dim(vacc_hosp_cumulative) <- c(n_age, n_regions)
dim(vacc_deaths_cumulative) <- c(n_age, n_regions)
dim(local_foi_region) <- n_regions
dim(annual_foi_region) <- n_regions
dim(rt_region) <- n_regions
dim(r0_region) <- n_regions
dim(seeding_value) <- n_regions



###### Initial states and dimensions for added code #####
initial(vacc_seroneg_cases_year[,]) <- 0
initial(vacc_seroneg_hosp_year[,]) <- 0
initial(vacc_seroneg_deaths_year[,]) <- 0
initial(vacc_seroneg_cases_cumulative[,]) <- 0
initial(vacc_seroneg_hosp_cumulative[,]) <- 0
initial(vacc_seroneg_deaths_cumulative[,]) <- 0
initial(vacc_seropos_cases_year[,]) <- 0
initial(vacc_seropos_hosp_year[,]) <- 0
initial(vacc_seropos_deaths_year[,]) <- 0
initial(vacc_seropos_cases_cumulative[,]) <- 0
initial(vacc_seropos_hosp_cumulative[,]) <- 0
initial(vacc_seropos_deaths_cumulative[,]) <- 0

dim(vacc_seroneg_cases_i) <- c(n_age, n_regions)
dim(vacc_seroneg_hosp_i) <- c(n_age, n_regions)
dim(vacc_seroneg_deaths_i) <- c(n_age, n_regions)
dim(vacc_seroneg_cases_year) <- c(n_age, n_regions)
dim(vacc_seroneg_hosp_year) <- c(n_age, n_regions)
dim(vacc_seroneg_deaths_year) <- c(n_age, n_regions)
dim(vacc_seroneg_cases_cumulative) <- c(n_age, n_regions)
dim(vacc_seroneg_hosp_cumulative) <- c(n_age, n_regions)
dim(vacc_seroneg_deaths_cumulative) <- c(n_age, n_regions)
dim(vacc_seropos_cases_i) <- c(n_age, n_regions)
dim(vacc_seropos_hosp_i) <- c(n_age, n_regions)
dim(vacc_seropos_deaths_i) <- c(n_age, n_regions)
dim(vacc_seropos_cases_year) <- c(n_age, n_regions)
dim(vacc_seropos_hosp_year) <- c(n_age, n_regions)
dim(vacc_seropos_deaths_year) <- c(n_age, n_regions)
dim(vacc_seropos_cases_cumulative) <- c(n_age, n_regions)
dim(vacc_seropos_hosp_cumulative) <- c(n_age, n_regions)
dim(vacc_seropos_deaths_cumulative) <- c(n_age, n_regions)

