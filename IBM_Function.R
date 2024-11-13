#' Run the simulation (wrapper)
#' 
#' Runs one simulation replicate for the IBM given a set of parameters.
#'
#' @param init_pop_size number of individuals (total of males and females) at the begining of the simulation
#' @param start year in which the simulation starts
#' @param end year in which the simulation ends 
#' @param fecundity number of surviving offspring per female and per year. In case of selection on juvenile survival it may be more consistent to input fertility here (that is, the number of offpsring, irrespective of survival)
#' @param reprovarf unexplained variance among males in latent reproductive success
#' @param reprovarm unexplained variance among males in latent reproductive success
#' @param VP phenotypic variance
#' @param VA additive genetic variance
#' @param maturation_rate maturation, a.k.a. recruitment, rate
#' @param surv_j baseline survival probability for the first stage class ("juveniles") with a phenotype of 0
#' @param surv_a baseline survival probability for the second stage class ("adults") with a phenotype of 0
#' @param beta_s_j slope of selection on juvenile survival
#' @param beta_s_a slope of selection on adult survival
#' @param beta_f slope of selection on fecundity
#' @param beta_m slope of selection on maturation
#' @param init_prop_juv initial proportion of juveniles in the population 
#' @param K "carrying capacity" or rather, population size threshold above which mean survival and reproduction are reduced in proportion of the excess above that threshold
#' @param repro_model reproductive model among "poisson", "bernoulli", "foldednormal"
#'
#' @return
#' @export
#'
#' @examples
main_simul <- function(init_pop_size=200, start=1, end=100, 
                       fecundity=2, reprovarf=0.001, reprovarm=0.001, 
                       VP=1, VA=0.2,
                       maturation_rate=0.3,
                       surv_j=0.4, surv_a=0.83,
                       beta_s_j=0, beta_s_a=0, beta_f=0, beta_m=0,
                       init_prop_juv=0.58, K=K, repro_model="poisson"){
  
  if(VP<VA){stop("Phenotypic variance cannot be less than additive genetic variance. Adjust VP and VA.")}
  
  pop <- InitialisePopulation(nind=init_pop_size, VA = VA, VP = VP, year = start, init_prop_juv=init_prop_juv)
  pop <- RunPopulation(pop = pop, start=start, end=end,
                       fecundity = fecundity, reprovarf = reprovarf, reprovarm = reprovarm,
                       VP=VP, VA=VA,
                       maturation_rate = maturation_rate, surv_j = surv_j, surv_a = surv_a,
                       beta_s_j=beta_s_j, beta_s_a=beta_s_a, beta_f=beta_f, beta_m=beta_m,
                       K=K, repro_model=repro_model)
  pop <- add_repro(pop)
  return(pop)
}#end main_simul()

#' Initialise the population
#' 
#' Initialise the population before the main simulation. Not intended for direct use.
#'  
#' @param nind number of individuals (total of males and females)
#' @param VA additive genetic variance
#' @param VP phenotypic variance
#' @param year starting year (e.g., 1 or 0)
#' @param init_age age of individuals the founding population
#' @param init_prop_juv proportion of immature individuals in the founding population
#'
#' @return
#' @export
#'
#' @examples
#' 
#' InitialisePopulation()
#' 
InitialisePopulation <- function(nind=100, VA=0.2, VP=1, year=1, init_age=2, init_prop_juv=0.5)
{
  VE <- VP-VA
  pop <- data.frame(id=1:nind, sex=c(rep("F", times=round(nind/2)),
                                     rep("M", times=nind - round(nind/2))),
                    birthyear = year-init_age,
                    deathyear = NA,
                    alive = 1,
                    mature = sample(c(0,1), size = nind, replace = TRUE, prob=c(init_prop_juv, 1-init_prop_juv)),
                    maturation = sample(c(NA,year-init_age + 1), size = nind, replace = TRUE, prob=c(init_prop_juv, 1-init_prop_juv)),
                    dam = NA,
                    sire = NA,
                    a = rnorm(nind, 0, sd= sqrt(VA)),
                    e = rnorm(nind, 0, sd= sqrt(VE)))
  pop$z <- pop$a+pop$e
  
  return(pop)
}#end InitialisePopulation()

#' Compute survival probability
#'
#' Compute survival probability on a given time-step given an individual age class, their phenotype and selective pressures. Not intended for direct use.
#' 
#' @param z phenotypic value
#' @param mature boolean, is the individual mature
#' @param beta_s_j slope of selection on juvenile survival
#' @param beta_s_a slope of selection on adult survival
#' @param surv_j baseline survival probability for the first stage class ("juveniles") with a phenotype of 0
#' @param surv_a baseline survival probability for the second stage class ("adults") with a phenotype of 0
#'
#' @return
#' @export
#'
#' @examples
#' survfuction()
survfuction <- function(z=0, mature=0, beta_s_j=0, beta_s_a=0,
                        surv_j=0.4, surv_a=0.83) {
  
  surv <- ifelse(mature==0, plogis(log(surv_j/(1-surv_j)) + beta_s_j*z), plogis(log(surv_a/(1-surv_a)) + beta_s_a*z) )
  
  return(surv)
}#end survfunction()


#' Main simulation function
#' 
#' Main simulation function but requires (i) the input of InitialisePopulation() and (ii) its output to be completed by add_repro(). Use main_simul() to run the IBM instead.
#'
#' @param start 
#' @param end 
#' @param fecundity 
#' @param reprovarf 
#' @param reprovarm 
#' @param VP 
#' @param VA 
#' @param maturation_rate 
#' @param surv_j 
#' @param surv_a 
#' @param beta_s_j 
#' @param beta_s_a 
#' @param beta_f 
#' @param beta_m 
#' @param pop 
#' @param K 
#' @param repro_model 
#'
#' @return
#' @export
#'
#' @examples
RunPopulation <- function(start=1, end=100,
                          fecundity=2, reprovarf=0, reprovarm=0.001, 
                          VP=1, VA=0.2,
                          maturation_rate=0.3,
                          surv_j=0.4, surv_a=0.83, 
                          beta_s_j=0, beta_s_a=0, beta_f=0, beta_m=0,
                          pop, K=5000, repro_model="poisson") 
{
  VE <- VP-VA
  
  for (year in start:end)
  {
    # reproduction
    newborn <- data.frame()
    
    if(sum(pop$alive)>K)
    {
      warning(paste("Population is above carrying capacity at the beginning of year", year, ". I am taking counter-measures."))
      regulation <- K/sum(pop$alive)
    }else{regulation <- 1}
    
    adultfemales <- which(pop$alive==1 & pop$sex=="F" & pop$mature==1)
    adultmales <- which(pop$alive==1 & pop$sex=="M" & pop$mature==1)
    
    if(length(adultfemales)*length(adultmales) > 0) #is any repro possible at all?
    {
      if(repro_model=="poisson")
      {offsp <- rpois(length(adultfemales),
                     lambda = exp(log(fecundity*regulation) + beta_f * pop$z[adultfemales] + rnorm(length(adultfemales), 0, sqrt(reprovarf))))
      }else{
        if(repro_model=="bernoulli")
        {
          if(fecundity*regulation > 1){stop("bernoulli reproductive mode is incompatible with a fecundity of 1 or more")}
          
          offsp <- rbinom(length(adultfemales), size = 1, prob = plogis(log(fecundity*regulation / (1-fecundity*regulation))  + beta_f * pop$z[adultfemales] + rnorm(length(adultfemales), 0, sqrt(reprovarf)) ) )
          
        }else{
          if(repro_model=="foldednormal")
          {
            offsp <- round(abs(rnorm(length(adultfemales), mean = fecundity*regulation + beta_f * pop$z[adultfemales], sd = sqrt(reprovarf))))
          }else{stop("repro_model not known")}
        }
      }
      
      malereprofitness <- cumsum(exp(rnorm(length(adultmales), 0, sqrt(reprovarm))))
      
      if( sum(offsp)>0) #is there realised reproduction at all
      {
        newborn <- data.frame(id=(max(pop$id)+1):(max(pop$id)+sum(offsp)), 
                              sex=sample(c("F", "M"), size = sum(offsp), replace = TRUE),
                              birthyear = year,
                              deathyear = NA,
                              alive = 1,
                              mature = 0,
                              maturation=NA,
                              dam = NA,
                              sire = NA,
                              a = 0,
                              e = rnorm(sum(offsp), 0, sd= sqrt(VE)))
        
        for(o in 1:sum(offsp))
        {
          currentoffsp <- min(which(cumsum(offsp) >= runif(1, min=0, max=sum(offsp))))
          newborn$dam[o] <- adultfemales[currentoffsp]
          offsp[currentoffsp] <- offsp[currentoffsp] - 1
          newborn$sire[o] <-adultmales[min(which(malereprofitness >= runif(1, min=0, max=max(malereprofitness))))]
          
          newborn$a[o] <- rnorm(n=1,
                                mean = (pop$a[pop$id==newborn$dam[o]] + pop$a[pop$id==newborn$sire[o]]) /2, 
                                sd = sqrt(VA/2))
          
        }# for(o in 1:sum(offsp))
        
        newborn$z <- newborn$a + newborn$e
        
        if(beta_s_j != 0) # facultative step to avoid unnecessary computations in case no juv selection
        {
          surviving_newborn <- as.logical(rbinom(n = nrow(newborn), 
                              size = 1,
                              prob = survfuction(newborn$z,
                                                 mature = rep(0,times=nrow(newborn)), # rep is crucial to return the right number of survival probabilities
                                                 surv_j = surv_j*regulation,
                                                 surv_a = 0,
                                                 beta_s_j = beta_s_j, beta_s_a = 0)))
          newborn <- newborn[surviving_newborn,]
          # need to rename ids to avoid shifting frame between ids and positions later on
          newborn$id <- (max(pop$id)+1):(max(pop$id)+sum(surviving_newborn))
        }
      }#end if( sum(offsp)>0) 
    } # end any repro at all
    
    # viability regulation and selection
    
    survivors <- rbinom(n = sum(pop$alive), 
                        size = 1,
                        prob = survfuction(pop$z[pop$alive==1],
                                           pop$mature[pop$alive==1],
                                           surv_j = surv_j*regulation,
                                           surv_a = surv_a*regulation,
                                           beta_s_j = beta_s_j, beta_s_a = beta_s_a ))
    
    pop$deathyear[pop$alive==1][which(survivors==0)] <- year
    pop$alive[pop$alive==1] <- survivors
    
    #maturation
    chance_to_mature <- which(pop$mature==0 & pop$alive==1)
    
    if(length(chance_to_mature)>0){
      pop$mature[chance_to_mature] <- rbinom(n = length(chance_to_mature), size = 1, 
                                             prob = plogis(log(maturation_rate/(1-maturation_rate)) + beta_m * pop$z[chance_to_mature] ) )
      pop$maturation[chance_to_mature][pop$mature[chance_to_mature]==1] <- year
    }
    
    if( nrow(newborn)>0) #if any reproduction at all
    {
      pop <- rbind(pop, newborn)
      rm(newborn) #crucial to avoid bug when pop close to extinction and not reproducing
    }
    
    sexofsurvivors <- pop$sex[pop$alive==1]
    if((sum(sexofsurvivors=="F")*sum(sexofsurvivors=="M"))==0){
      warning(paste("Population went functionally exctinct on year", year))
      return(pop)
    }
    
  }#end  for (year in start:end)
  
  return(pop)
}#end RunPopulation()


#' Add information about individual lifetime reproductive success at the end of the simulation
#'
#' @param pop 
#'
#' @return
#' @export
#'
#' @examples
add_repro <- function(pop)  
{
  pop$repro <- 0
  pop2 <- pop[!is.na(pop$dam) & !is.na(pop$sire),]
  for (i in 1:nrow(pop))
  {
    if(pop$sex[i]=="F")
    {
      pop$repro[i] <- nrow(pop2[pop2$dam==pop$id[i],])
    }else{
      pop$repro[i] <- nrow(pop2[pop2$sire==pop$id[i],])
    }
  }
  return(pop)
}#end add_repro()



#' Compute realized generation time
#'
#' @param pop 
#'
#' @return
#' @export
#'
#' @examples
generation_time <- function(pop){
  agesatbirth <- vector()
  for (i in 1:nrow(pop))
  {
    if(!is.na(pop$dam[i]))
    {
      agesatbirth <- c(agesatbirth, pop$birthyear[i] - pop$birthyear[pop$id==pop$dam[i]])
    }
    if(!is.na(pop$sire[i]))
    {
      agesatbirth <- c(agesatbirth, pop$birthyear[i] - pop$birthyear[pop$id==pop$sire[i]])
    }
  }
  return(mean(agesatbirth))
}#end generation_time()    


#' Compute population size for every year in the simulation
#'
#' @param pop 
#' @param maxyear 
#' @param minyear 
#'
#' @return
#' @export
#'
#' @examples
produce_pop_data <- function(pop, maxyear=NULL, minyear=NULL)
{
  years <- vector()
  
  pop$deathyear[is.na(pop$deathyear)] <- max(pop$birthyear)
  for(i in 1:nrow(pop))
  {
    years <- c(years, pop$birthyear[i]:pop$deathyear[i])
  }
  popsizes <- as.data.frame(table(years))
  popsizes$years <- as.numeric(as.character(popsizes$years))
  
  if(!is.null(minyear)){
    popsizes <- popsizes[popsizes$years>=minyear,]
  }
  
  colnames(popsizes)[2] <- "size"
  if(!is.null(maxyear)){
    if(max(popsizes$years)>maxyear){
      popsizes <- popsizes[popsizes$years<=maxyear]
    }
    if(max(popsizes$years)<maxyear){
      popsizes <- rbind(popsizes, data.frame(years=(max(popsizes$years)+1):maxyear, size=0))
    }
  }
  return(popsizes)
}#end produce_pop_data()

#' Compute geometric average population growth rate
#'
#' @param pop_data 
#' @param exclude 
#'
#' @return
#' @export
#'
#' @examples
geom_pop_growth <- function(pop_data, exclude=5)
{
  y_gr <- pop_data$size[(exclude+1):nrow(pop_data)]/pop_data$size[exclude:(nrow(pop_data)-1)]
  prod(y_gr)^(1/length(y_gr))
}#end geom_pop_growth()


#' Compute breeding value average per cohort
#'
#' @param pop 
#' @param maxyear 
#' @param minyear 
#' @param burnin 
#'
#' @return
#' @export
#'
#' @examples
produce_evolution_data <- function(pop, maxyear=NULL, minyear=NULL, burnin=0)
{
  
  if(!is.null(minyear)){
    pop <- pop[pop$birthyear>=burnin,]}
  if(!is.null(maxyear)){
    pop <- pop[pop$birthyear<=maxyear,]}
  
  require(magrittr)
  require(dplyr)
  
  evol <- pop %>% 
    group_by(birthyear, .drop=FALSE) %>%
    summarise(mean_pheno = mean(z), mean_bv = mean(a))
  
  year_range <- min(pop$birthyear):max(pop$birthyear)
  not_sampled <-  year_range[!year_range %in% evol$birthyear]
  if(length(not_sampled) > 0)
  {
    evol <- rbind(evol, data.frame(birthyear = not_sampled, mean_pheno=NA, mean_bv=NA))
  }
  
  if(!is.null(maxyear)){
    if(max(evol$birthyear)<maxyear){
      evol <- rbind(evol, data.frame(birthyear=(max(evol$birthyear)+1):maxyear, mean_pheno=NA, mean_bv=NA))
    }
  }
  return(evol)
}#end produce_pop_data()



#' Compute breeding value average per year
#'
#' @param pop 
#' @param maxyear 
#' @param minyear 
#' @param burnin 
#'
#' @return
#' @export
#'
#' @examples
produce_evolution_data_byyear <- function(pop, maxyear=NULL, minyear=NULL, burnin=0)
{
  
  if(!is.null(minyear)){
    pop <- pop[pop$birthyear>=burnin,]}
  if(!is.null(maxyear)){
    pop <- pop[pop$birthyear<=maxyear,]}
  
  require(magrittr)
  require(dplyr)
  
  pop$deathyear[is.na(pop$deathyear)] <- max(pop$deathyear, na.rm = TRUE)
  
  pop_expanded <- pop %>%
    group_by(id) %>% 
    uncount(deathyear-birthyear, .id="age") %>%
    mutate(year=birthyear + age - 1)
  
  evol <- pop_expanded %>% 
    filter(year>=0) %>%
    group_by(year, .drop=FALSE) %>%
    summarise(mean_pheno = mean(z), mean_bv = mean(a))
  
  year_range <- min(pop_expanded$year):max(pop_expanded$year)
  not_sampled <-  year_range[!year_range %in% evol$year]
  if(length(not_sampled) > 0)
  {
    evol <- rbind(evol, data.frame(year = not_sampled, mean_pheno=NA, mean_bv=NA))
  }
  
  if(!is.null(maxyear)){
    if(max(evol$year)<maxyear){
      evol <- rbind(evol, data.frame(year=(max(evol$year)+1):maxyear, mean_pheno=NA, mean_bv=NA))
    }
  }
  return(evol)
}#end produce_pop_data()




#' Produce individual data
#'
#' @param pop 
#'
#' @return
#' @export
#'
#' @examples
produce_ind_data <- function(pop)
{
  pop$deathyear[is.na(pop$deathyear)] <- max(pop$birthyear)
  nbrecords <- pop$deathyear - pop$birthyear #don't count juv when they are created
  cumrecords <- c(1,cumsum(nbrecords))
  inddata <- data.frame(obs = 1: (sum(nbrecords)), year=0, id=0, sex="F", z=0.00, birthyear=0, a=0.00, age=0, surv=1, mature=0, repro=0) 
  
  for (i in 1:nrow(pop))
  {
    t_rows <- cumrecords[i]:(cumrecords[i+1]-1)
    inddata[t_rows ,c("id", "sex", "z", "birthyear", "a")] <- pop[i, c("id", "sex", "z", "birthyear", "a")]
    inddata[t_rows, "age"] <- 0:(length(t_rows)-1)
    inddata[t_rows, "year"] <- inddata[t_rows, "age"] + inddata[t_rows, "birthyear"]
    inddata[t_rows[length(t_rows)], "surv"] <- pop[i,"alive"]
    
    if(!is.na(pop$maturation[i]))
    {
      maturation_age <- pop$maturation[i]-pop$birthyear[i]
      inddata[t_rows, "mature"] <- inddata[t_rows,"age"] >= maturation_age
      
      offsp <- which(pop$id[i] == pop$dam | pop$id[i] == pop$sire)
      if(length(offsp)>0)
      {
        repro_calendar <- table(pop[offsp,"birthyear"])
        for (j in 1:length(repro_calendar))
        {
          inddata[t_rows,][inddata[t_rows, "year"] == as.numeric(names(repro_calendar))[j],"repro"] <- repro_calendar[j]
        }
      }
    }
    
  }
  
  return(inddata)
}#end produce_ind_data()

