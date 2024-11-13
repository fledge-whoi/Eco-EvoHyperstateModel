library(tidyverse)
library(ggplot2)
curve(expr = exp(log(11) + x*0.5), from = -4, to = 25)

# Parameterize the life cycles of 5 distinct species.
Data <- list(
  c(11,0.0365,0.0965,0.9),
  c(5, 0.2, 0.25, 0.572),
  c(1,0.8,0.385,0.4),
  c(0.3, 0.93, 0.505, 0.3), 
  c(0.232, 0.95, 0.80, 0.07)
) # (F, SA, SJ, Y)

# Quantitative genetics information
Vp <- 1 # Phenotypic variance
Va <- 0.05 # Heritability value of the phenotypic trait

# Set selective pressure (slope of the relationship between vital rate and phenotype)
beta <- 0.01

Fresponses <- function(Va=0.1, Vp=1, beta=0.15, Data){
nspecies <- length(Data) # number of species
h2 <- Va / Vp # Additive genetic variance (h2 = Va/Vp)
Ve <- Vp - Va # Environmental variance (Vp = Va + Ve)
# Get demographic descriptive info
DEMOout <- matrix(0, nspecies, 3) # outputs lambda, R0, Tgen
S <- matrix(0, nspecies, 5) # sensitivity d lambda/ d theta

dfout <- expand.grid(unit=c("year", "generation"), rate=c("f", "sa", "sj", "m"), selection= c("B", "theta", "R0"), species=1:5)
dfout$Response <- NA

# Main loop for species
for (ispecies in 1:nspecies) {
  # Extract vital rates of the species
  F  <- Data[[ispecies]][1] # offspring produced per year that survive the 1st year
  SA <- Data[[ispecies]][2]
  SJ <- Data[[ispecies]][3]
  Y  <- Data[[ispecies]][4]
  f  <- F / SJ # offspring produced per year
  
  # Projection matrix
  A <- matrix(c(SJ * (1 - Y), f * SJ, SJ * Y, SA), nrow = 2, byrow = TRUE)
  
  # Population growth rate / fitness
  lambda <- (SJ * (1 - Y) + SA) / 2 + sqrt((SJ * (1 - Y) - SA)^2 + 4 * f * SJ^2 * Y) / 2
  
  # Eigenvectors to calculate sensitivity
  v <- c(Y * SJ, lambda - SJ * (1 - Y))
  w <- c(f * SJ, lambda - SJ * (1 - Y))
  
  # Descriptive demographic information
  R0 <- Y * f * SJ^2 / ((1 - SA) * (1 - SJ * (1 - Y)))
  Tgen <- lambda * sum(v * w) / (v[1] * w[2] * f * SJ)
  
  # Sensitivity matrix
  Sensitivity <- outer(v, w) / sum(v * w)
  
  # Calculate sensitivities for each parameter
  theta <- c(f, SJ, SA, Y)
  PDlv <- array(0, dim = c(2, 2, 4))
  
  # Partial derivatives
  PDlv[1, 2, 1] <- SJ
  PDlv[1, 1, 2] <- 1 - Y
  PDlv[1, 2, 2] <- F
  PDlv[2, 1, 2] <- Y
  PDlv[2, 2, 3] <- 1
  PDlv[1, 1, 4] <- -SJ
  PDlv[2, 1, 4] <- SJ
  
  # Calculate sensitivities S(ispecies, j)
  for (j in 1:length(theta)) {
    S[ispecies, j] <- sum(PDlv[,,j] * Sensitivity)
  }
  # Calculate sensitivity of lambda to F
  temp <- matrix(0, 2, 2)
  temp[1, 2] <- 1
  S[ispecies, 5] <- sum(temp * Sensitivity)
  
  # Store outputs in DEMOout and RAout
  DEMOout[ispecies, ] <- c(lambda, R0, Tgen)
  
  # fecundity
  # per time
  dfout$R[dfout$unit=="year" & dfout$rate=="f" & dfout$selection=="B" & dfout$species==ispecies]  <- h2 * Vp * beta / Tgen
  dfout$R[dfout$unit=="year" & dfout$rate=="f" & dfout$selection=="theta" & dfout$species==ispecies] <- h2 * Vp / Tgen
  dfout$R[dfout$unit=="year" & dfout$rate=="f" & dfout$selection=="R0" & dfout$species==ispecies] <- h2 * Vp / Tgen
  # per generation
  dfout$R[dfout$unit=="generation" & dfout$rate=="f" & dfout$selection=="B" & dfout$species==ispecies] <- h2 * Vp * beta
  dfout$R[dfout$unit=="generation" & dfout$rate=="f" & dfout$selection=="theta" & dfout$species==ispecies] <- h2 * Vp
  dfout$R[dfout$unit=="generation" & dfout$rate=="f" & dfout$selection=="R0" & dfout$species==ispecies] <- h2 * Vp
  
  # juvenile survival
  # per time
  dfout$R[dfout$unit=="year" & dfout$rate=="sj" & dfout$selection=="B" & dfout$species==ispecies] <- h2 * Vp * beta * (1 / Tgen) * (1 + (lambda / (lambda - SJ * (1 - Y)))) * (1 - SJ)
  dfout$R[dfout$unit=="year" & dfout$rate=="sj" & dfout$selection=="theta" & dfout$species==ispecies] <- h2 * Vp * (1 / Tgen) * (1 + (lambda / (lambda - SJ * (1 - Y))))
  dfout$R[dfout$unit=="year" & dfout$rate=="sj" & dfout$selection=="R0" & dfout$species==ispecies] <-  h2 * Vp * (1 / Tgen) * ((1 + (lambda / (lambda - SJ * (1 - Y)))) / (1 + (1 / (1 - SJ * (1 - Y)))))
  # per generation
  dfout$R[dfout$unit=="generation" & dfout$rate=="sj" & dfout$selection=="B" & dfout$species==ispecies] <- h2 * Vp * beta * (1 + (lambda / (lambda - SJ * (1 - Y)))) * (1 - SJ)
  dfout$R[dfout$unit=="generation" & dfout$rate=="sj" & dfout$selection=="theta" & dfout$species==ispecies] <- h2 * Vp * (1 + (lambda / (lambda - SJ * (1 - Y))))
  dfout$R[dfout$unit=="generation" & dfout$rate=="sj" & dfout$selection=="R0" & dfout$species==ispecies] <- h2 * Vp * ((1 + (lambda / (lambda - SJ * (1 - Y)))) / (1 + (1 / (1 - SJ * (1 - Y)))))
  
  # maturation rate
  # per time
  dfout$R[dfout$unit=="year" & dfout$rate=="m" & dfout$selection=="B" & dfout$species==ispecies] <- h2 * Vp * beta * (1 / Tgen) * (((lambda - SJ) * (1 - Y)) / (lambda - (SJ * (1 - Y))))
  dfout$R[dfout$unit=="year" & dfout$rate=="m" & dfout$selection=="theta" & dfout$species==ispecies] <- h2 * Vp * (1 / Tgen) * (lambda - SJ) / (lambda - SJ * (1 - Y))
  dfout$R[dfout$unit=="year" & dfout$rate=="m" & dfout$selection=="R0" & dfout$species==ispecies] <-h2 * Vp * (1 / Tgen) * ((lambda - SJ) * (1 - SJ * (1 - Y)) / ((1 - SJ) * (lambda - SJ * (1 - Y))))
  # per generation
  dfout$R[dfout$unit=="generation" & dfout$rate=="m" & dfout$selection=="B" & dfout$species==ispecies] <- h2 * Vp * beta * (((lambda - SJ) * (1 - Y)) / (lambda - (SJ * (1 - Y))))
  dfout$R[dfout$unit=="generation" & dfout$rate=="m" & dfout$selection=="theta" & dfout$species==ispecies] <- h2 * Vp * (lambda - SJ) / (lambda - SJ * (1 - Y))
  dfout$R[dfout$unit=="generation" & dfout$rate=="m" & dfout$selection=="R0" & dfout$species==ispecies] <- h2 * Vp * ((lambda - SJ) * (1 - SJ * (1 - Y)) / ((1 - SJ) * (lambda - SJ * (1 - Y))))
  
  # adult survival
  # per time
  dfout$R[dfout$unit=="year" & dfout$rate=="sa" & dfout$selection=="B" & dfout$species==ispecies] <- h2 * Vp * beta * (1 / Tgen) * ((SA * (1 - SA)) / (lambda - SA))
  dfout$R[dfout$unit=="year" & dfout$rate=="sa" & dfout$selection=="theta" & dfout$species==ispecies] <- h2 * Vp * (1 / Tgen) * (SA / (lambda - SA))
  dfout$R[dfout$unit=="year" & dfout$rate=="sa" & dfout$selection=="R0" & dfout$species==ispecies] <- h2 * Vp * (1 / Tgen) * ((1 - SA) / (lambda - SA))
  # per generation
  dfout$R[dfout$unit=="generation" & dfout$rate=="sa" & dfout$selection=="B" & dfout$species==ispecies] <- h2 * Vp * beta * ((SA * (1 - SA)) / (lambda - SA))
  dfout$R[dfout$unit=="generation" & dfout$rate=="sa" & dfout$selection=="theta" & dfout$species==ispecies] <- h2 * Vp * (SA / (lambda - SA))
  dfout$R[dfout$unit=="generation" & dfout$rate=="sa" & dfout$selection=="R0" & dfout$species==ispecies] <- h2 * Vp * ((1 - SA) / (lambda - SA))
}
return(dfout)
}#end Fresponses

Vas <- c(0.05, 0.2, 0.8)
betas <- c(0.01, 0.15, 0.5)
  
resall <- data.frame()

for (i in 1:length(Vas))
  for(j in 1:length(betas))
  {
    out <- Fresponses(Va=Vas[i], beta = betas[j], Data = Data)
    out$Va <- Vas[i]
    out$beta <- betas[j]
    if(i==1 & j==1){
      resall <- out
    }else{
    resall <- rbind(resall, out)}
  }
resall$rate <- factor(x = resall$rate, levels = c("f", "sj", "sa", "m"))

resall %>% filter(unit=="year", selection=="B", Va==0.05, beta==0.01)

resall %>% filter(unit=="year", selection=="B") %>% ggplot(aes(x=species, y=100*R, fill=rate))  + 
  geom_col(position = position_dodge2()) +
  facet_wrap(~Va + beta, scales = "free", nrow = 3, labeller = "label_both") +
  theme_bw()


library(ggplot2)
resall %>% 
  filter(Va==0.2, beta==0.15) %>%
  ggplot( aes(x=species, y=R, fill=rate)) + 
  geom_col(position = position_dodge2()) + facet_wrap(~selection +unit, scales = "free", nrow = 3) +
  theme_bw()

resall %>% 
  filter(Va==0.2, beta==0.15, unit == "year", selection == "B") %>% mutate(Response=100*R) %>% write_csv(file = "rates.csv")



resall$species <- as.factor(resall$species)

resall %>% 
  filter(Va==0.2, beta==0.5, selection=="B") %>%
  ggplot( aes(x=species, y=R, fill=rate)) + 
  geom_col(position = position_dodge2()) + 
  scale_x_discrete(labels=c("1" = "grass", "2" = "rodent",
                            "3" = "passerine", "4" = "deer", "5" = "albatross")) +
  scale_y_continuous(name = "Genetic change (in phenotypic standard deviation)") + 
  scale_fill_discrete(name = "Vital rate selected", labels = c("fertility", "juvenile survival", "adult survival", "maturation"), type =
                        wesanderson::wes_palette("Moonrise2"))+
  facet_wrap(~unit, scales = "free", nrow = 3, 
             labeller = labeller(unit = c(year = "per year adaptation rate", generation="per generation adaptation rate"))) +
  theme_bw()

ggsave(filename = "compadaptrates.png", width = 10, height = 8)

scale_fill_manual() +c("red", "#d8b365", "#f5f515", "#5ab4ac")
c("#2F4F4F", "#EEB422", "#836FFF", "#CD2626")