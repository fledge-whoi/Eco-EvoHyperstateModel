#' Emat
#'
#' Create a square matrix with a 1 and all other values being 0
#'
#' @param i row coordinate of the 1
#'
#' @param j column coordinate of the 1
#'
#' @param u dimension of the squared matrix
#'
#' @examples
#' Emat(2,3,4)
#'
#' @export
Emat <- function(i,j,u){
  E <- matrix(data=0, nrow=u, ncol=u)
  E[i,j]=1
  return(E)
}#end Emat()


#' BD_proj_mat
#' 
#' Create a block-diagonal matrix from an array
#' 
#' @param B array
#' 
#' @example 
#' 
#' B <- array(runif(100, min=1, max=8), dim=c(2,2,5,5))
#' image(BD_proj_mat(B=B))
#' 
#' @export
BD_proj_mat <- function(B){
  siz=dim(B)[-1]
  s=prod(siz);                    #size of the expected block-diagonal matrix
  sk=s/siz[1];                    #number of blocks on the diagonal
  siz=siz[2:3];                 #maximal value of each index i_1,...,i_r
  A=matrix(data = 0, nrow = s, ncol=s)
  
  for (i in 1:sk)
  {
    r = ((i-1) %% siz[1]) + 1
    c = floor((i-1) / siz[1]) + 1
    A=A+kronecker(Emat(i,i,sk), B[ , , r, c] )
  }
  return(A)
}#end BD_proj_mat()

#' speye
#' 
#' Create a square identity sparse matrix
#' 
#' @param N Matrix dimension (number of rows, or equivalently numbers of column)
#' 
#' @export
speye <- function(N){
    require(Matrix)
    return(sparseMatrix(i=(1:N),j=(1:N), x=rep(1,N) ) )
  }#end speye()



#' GVLE
#' 
#' Function to get the probability density distribution of breeding values from the normal distribution given the additive genetic variance
#' 
#' @param zz breeding values for which to compute the density
#' 
#' @param Va additive genetic variance
#' 
#' @details
#' There is not a factor 2 because the variance we want is the Mendelian segregation variance.
#' 
#' 
#' @export
GVLE <- function(zz, Va){
  exp(-(zz^2)/Va)/sqrt(pi*Va)
}


#' MPM
#' 
#' Run a eco-evo deterministic matrix model
#' 
#' @param Beta numeric, strength of selection on the scale of linear predictors
#' 
#' @param h2 numeric on [0,1], heritability of the trait under selection
#' 
#' @param xp numeric, the range of phenotypic values considered
#' 
#' @param timeA numeric, number of time steps of adaptation. The simulations has two more time steps.
#' 
#' @param g number of phenotypic and genetic classes considered
#' 
#' @param Vp phenotypic variance
#' 
#' @param nind number of individuals to simulate, purely a matter of scale as there is no stochasticity
#' 
#' @param  Fe Fertility
#'
#' @param SA Adult survival rate
#' 
#' @param SJ Juvenile survival rate
#' 
#' @param S0 New born survival rate
#' 
#' @param Y maturation rate
#' 
#' @param selectedrate the vital rate under selection
#' 
#' @param density_tolerance Threshold to trigger warnings about the distribution of breeding values being too narrow or too wide
#'
#' @param enforce_initial_lambda1 Boolean, default FALSE. Should the initial asymptotic population growth be exactly 1.
#' 
#' @param stabilise_population_structure Boolean, default TRUE. Should the simulation starts with stable stage structure.
#' 
#' 
#' @export
MPM <- function(Beta = 0.15, h2 = 0.2, xp=4 ,timeA=100, g = 20, Vp = 1,  nind = 100, 
                Fe=0.2286, SA=0.95, SJ=0.8,  S0=0.8, Y=0.07, selectedrate="Fe", selectmodel="log", 
                census = "prebreeding", density_tolerance = 10^(-5), enforce_initial_lambda1=FALSE, 
                stabilise_population_structure=TRUE){
  require(pracma) # for interpolation
  require(gsignal) # for 2D convolution
  require(purrr) # for reformatting output
  require(tidyr) # for reformatting output
  
  if(h2<=0 | h2>=1) {stop("h2 must be strictly greater than 0 and strictly smaller than 1")}
  if(h2 < ((xp/g)^2)) {warning("g is probably too low given h2 and xp. There are not enough breeding value classes to house Va below (xp/g)^2. Check estimates of Va and Vp in the output to see the realised h2. You may increase xp or decrease g.")}
  
  if(enforce_initial_lambda1){
    Fe <- (1-SA)*(1-(1-Y)*SJ)/(Y*SJ) ##### TO GET LAMBDA 1
    message(paste("I have adjusted Fe to", Fe, "because the enforce_initial_lambda1 is set to TRUE."))
  }
  
  w = 2 # number of stages
  m = 3 # three dimensions: stage, breeding value, phenotype
  nvital = 4 # number of vital rates (juvenile survival, maturation, adult survival, fertility)
  Va = h2*Vp # Additive genetic variance
  Ve = Vp - Va # Get the environmental variance (Vp = Va + Ve)
  
  # Create the phenotypic bins for the matrix
  min_pheno    = -xp 
  max_pheno    = xp 
  min_breeding = -xp 
  max_breeding = xp 
  mean_pheno   = 0 
  
  # Get the categories of phenotypes
  edges_e    = seq(from=min_pheno, to=max_pheno, length.out=g+1)
  binwidth_e = edges_e[2] - edges_e[1] # Bin width
  midpoint_e = edges_e[1:(length(edges_e)-1)] +  (binwidth_e/2) # Mid phenotype in each bin
  
  # Get the categories of breeding values
  edges_a    = edges_e[(edges_e<=max_breeding)&(edges_e>=min_breeding)]
  binwidth_a = edges_a[2] - edges_a[1] # Bin width
  midpoint_a = edges_a[1:(length(edges_a)-1)] +  (binwidth_a/2) # Mid phenotype in each bin
  
  b = length(midpoint_a) # this is the number of classes for breeding value
  midpoint2_a = seq(2*midpoint_a[1], 2*midpoint_a[length(midpoint_a)], by=binwidth_a)
  
  Gle <- GVLE(0.5*(midpoint2_a), Va = Va)
  if(Gle[1]>density_tolerance | Gle[length(Gle)] > density_tolerance) {warning("The initial distribution of breeding values reaches the boundaries of the possible range. You should probably decrease h2, decrease Vp, or increase xp.")}
  if((sum(Gle) - (1/(0.5 * binwidth_a))) > sqrt(density_tolerance)) {warning("The initial distribution of breeding values is too narrow given the coarsness of phenotypic classes; genetic change may not be computed accurately. You should probably increase h2, increase Vp, increase g or decrease xp. ")}
  
  # B matrix
  B = speye(b*w*g)

  # P matrix
  P = speye(b*w*g)

  # Simulate an initial dataframe for the vector of population size at t initial
  nind    = 100 # number of individuals to simulate
  
  N_ve = t(exp((-outer(midpoint_e, midpoint_a, FUN="-")^2 / (2*Ve) )) / colSums(exp(-outer(midpoint_e, midpoint_a, FUN="-")^2 / (2*Ve))))
  
  nbv  = exp(-midpoint_a^2/(2*Va))/sum(exp(-midpoint_a^2/(2*Va))) # Get the proportion of individuals within each class of breeding value
  Nbv  = nind*nbv # Multiply by the total number of individuals in the population to get the number of individuals in each class of breeding value
  nz   = N_ve %*% nbv * binwidth_a
  Nz   = nind*nz
  NN   = t(t(N_ve) * Nbv)
  
  ninit = aperm(array(data = NN/w, dim = c(g,g,w)), perm=c(3,2,1))
  Ninit = t(t(as.vector(ninit))) # here we compute initial population distribution without accounting for stable population structure. We will correct that later
  
  nit = timeA + 2 # number of years
  
  M_EBVS = vector(length = nit)
  M_EBVSJ = vector(length = nit)
  M_EBVSA = vector(length = nit)
  f_Z   = vector(length = nit)
  SJ_Z   = vector(length = nit)
  SA_Z   = vector(length = nit)
  Y_Z  = vector(length = nit)
  Va_hat = vector(length = nit)
  Vp_hat = vector(length = nit)
  
  A <- matrix(data=c(SJ*(1-Y), Fe, SJ*Y, SA), nrow=2, byrow = TRUE) #Leslie matrix
  f <- Fe/SJ
  eigenA <- eigen(A)
  LAMBDA0 <- lambda1_mean <- eigenA$values[1]
  # stable age-stage distribution, w = right eigenvector
  w_vec_mean <- eigenA$vector[,1] / sum(eigenA$vector[,1])
  v_vec_mean <- eigen(t(A))$vector[,1]  / sum(eigen(t(A))$vector[,1])
  S_lambda <- (v_vec_mean %*% t(w_vec_mean)) / as.vector(t(v_vec_mean)%*%w_vec_mean)
  T_GEN <- lambda1_mean*(t(v_vec_mean) %*% w_vec_mean)/(t(v_vec_mean) %*% matrix(c(0,Fe,0,0), nrow=2, byrow = TRUE) %*% w_vec_mean)
  
  # Get Theta (parameters)
  S0_pheno <- rep(S0, times=b)
  SJ_pheno <- rep(SJ, times=b)
  Y_pheno  <- rep(Y, times=b)
  SA_pheno <- rep(SA, times=b)
  f_pheno  <- rep(f, times=b)
  #F_pheno  <- f_pheno*SA_pheno; #POST
  F_pheno  <- f_pheno*SJ_pheno;
  
  if (selectedrate=="Fe"){
    if (selectmodel=="log"){
      # Poisson
      beta <- Beta
      f_pheno  <- exp(log(f) + beta*midpoint_e)
      dS <- Beta*f
    }else{
      if (selectmodel =="logit"){
        beta <- Beta;
        if( (f<=0) | (f>=1)) {stop("When using a logit model for selection the fertility rate (Fe/Sj) must be between 0 and 1 (excluded)")}
        f_pheno  <- 1/(1+exp(-(logit(f) + beta*midpoint_e)))
        dS <- Beta*f*(1-f)
      }
    }
    
    S_LAMBDA <- SJ*S_lambda[1,2]
    F_pheno  <- f_pheno*SJ_pheno #PRE
  }else{
    if(selectedrate=="SJ"){
      beta <- Beta
      SJ_pheno <- 1/(1+exp(-(logit(SJ) + beta*midpoint_e)))
      S_LAMBDA <- (1-Y)*S_lambda[1,1]+Y*S_lambda[2,1]+f*S_lambda[1,2]
      dS <- Beta*SJ*(1-SJ)
      F_pheno <- f_pheno*SJ_pheno
    }else{
      if(selectedrate=="SA"){
        beta <- Beta
        SA_pheno <- 1/(1+exp(-(logit(SA) + beta*midpoint_e)))
        #%F_pheno <- f_pheno*SA_pheno #POST
        #S_LAMBDA <- S_lambda[2,2]+f*S_lambda[1,2] # POST
        S_LAMBDA <- S_lambda[2,2] # PRE
        dS <- Beta*SA*(1-SA)
      }else{
        if(selectedrate=="Y"){
          beta <- Beta
          Y_pheno <- 1/(1+exp(-(logit(Y) + beta*midpoint_e)))
          S_LAMBDA <- SJ*(S_lambda[2,1]-S_lambda[1,1])
          dS <- Beta*Y*(1-Y)
        }else{stop("selectedrate not recognized")}
      }
    }
  }
  
  V_ADAPT <- S_LAMBDA*dS*Va
  
  # Get the matrices U, B, P for Utilde and R, H, M for Ftilde
  # U matrix
  idr = c(1,2,2)
  idc = c(1,1,2)
  if (selectedrate=="Fe"){
    Ujk = c(SJ*(1-Y),
            SJ*Y,
            SA)
    Ui = sparseMatrix(i = idr,j = idc, x = Ujk)
  }
  nId = length(idr)
  Idr = rep(0, times=nId*b*g)
  Idc = rep(0, times=nId*b*g)
  Uu = rep(0, times=nId*b*g)
  
  for (iu in 1:(b*g))
  {
    Ik = ((iu-1)*nId+1):(iu*nId)
    Ptype <- arrayInd(.dim = c(b,g), ind = iu)[,2]
    if (selectedrate=="SJ"){
      Ujk = c(SJ_pheno[Ptype]*(1-Y),SJ_pheno[Ptype]*Y, SA)
    }else{
      if (selectedrate == "SA"){
        Ujk = c(SJ*(1-Y),SJ*Y, SA_pheno[Ptype] )
      }else{
        if(selectedrate == "Y"){
          Ujk = c(SJ*(1-Y_pheno[Ptype]),SJ*Y_pheno[Ptype], SA)
        }
      }
    }
    Ui = sparseMatrix(i = idr,j = idc, x = Ujk)
    Idr[Ik] <- idr+(iu-1)*w
    Idc[Ik] <- idc+(iu-1)*w
    Uu[Ik] <- Ui[Ui!=0]
  }
  U = sparseMatrix(i = Idr, j = Idc, x = Uu, dims = c(b*w*g,b*w*g))
  
  
  #%% R matrix
  idr = 1
  idc = 2
  Rjk = F_pheno[1]
  Ri = sparseMatrix(i = idr,j = idc,x = Rjk)
  nId = length(idr)
  Idr = rep(0, times=nId*b*g)
  Idc = rep(0, times=nId*b*g)
  Rr   = rep(0,times=nId*b*g)
  for (ir in 1:(b*g))
  {
    Ik = ((ir-1)*nId+1):ir*nId;
    if (selectedrate %in% c("Fe", "SJ"))
    {
      Ptype = arrayInd(ind = ir, .dim = c(b,g))[,2]
      Rjk = F_pheno[Ptype]
      Ri = sparseMatrix(i = idr, j = idc, x = Rjk)
    }
    Idr[Ik] = idr+(ir-1)*w
    Idc[Ik] = idc+(ir-1)*w
    Rr[Ik] = Ri[Ri!=0]
  }
  R = sparseMatrix(i = Idr,j = Idc, x = Rr, dims = c(b*w*g,b*w*g))

  #%% H matrix
  #% To construct the H matrix, we need the Ga matrix
  Vle <- Va
  
  
  #%% M matrix
  #% To construct the M matrix, we need the Ge matrix
  #% For each breeding value, what is the potential distribution of offspring phenotype? - Here we add noise (Ve) to the breeding value to obtain offspring phenotype
  emat <- t(exp(-(outer(midpoint_e, midpoint_a, FUN="-"))^2 / (2*Ve) ) * binwidth_a / apply(exp(-(outer(midpoint_e, midpoint_a, FUN="-"))^2/(2*Ve))*binwidth_e, 1, sum))
  
  # The matrix M, which is the transmission of maternal phenotype to offspring phenotype goes through the breeding values only
   nId <-  g*g
  Idr <-   rep(0, times=nId*b)
  Idc <-  rep(0, times=nId*b)
  Mm  <-  rep(0, times=nId*b)
  inew <-  0
  for (i in 1:b)
  {
    iold <- inew
    idr <- which(emat[,i]!=0)
    Mi_n <- emat[idr,i]
    inew <-  length(idr)
    idc <- rep(1:g, each=inew)
    idr <- rep(idr, times=g)
    Ik <- ((i-1)*nId+1):(((i-1)*nId+1)+inew*g-1)
    Idr[Ik] <- idr+((i-1)*w*g)
    Idc[Ik] <- idc+((i-1)*w*g)
    Mm[Ik] <- rep(Mi_n, times=g)
  }
  Idr <- Idr[Idr != 0]
  Idc <- Idc[Idc != 0]
  Mm <- Mm[Mm != 0]
  M = sparseMatrix(i = Idr, j = Idc, x = Mm, dims = c(b*w*g,b*w*g))
  
  #%% Compute N(t) correcting for stable population structure
  if(stabilise_population_structure){
  nninit =array(data = aperm(NN, c(2,1)), dim = c(1,b,g))
  ninit  = outer(w_vec_mean, nninit, "*")
  Ninit =  t(t(as.vector(ninit)))
  }
  
  ninit = aperm(array(data = NN/w, dim = c(g,g,w)), perm=c(3,2,1))
  
  N <- matrix(0, nrow = w*b*g, ncol=nit)
  N[,1] <- Ninit
  NN  <- array(0, c(w,b,g,nit))
  NN[ , , ,1] <- ninit
  
  for (i in 2:nit)
  {
    # Live individual transitions
    # Process #1.1 = Transition of live individuals between states
    NU <- U %*% N[,i-1]
    
    # Rearrange for process #1.2
    nu <- array(NU,c(w,b,g))
    nu_r <- aperm(nu,c(2,1,3)) # bwg
    nu_r <- as.vector(nu_r)
    
    # Process #1.2 = Transition of live individuals between breeding value categories (assumed fixed)
    NB <- B %*% nu_r
    
    # Rearrange for process #1.3
    nb <- array(NB, c(b,w,g))
    nb_r <- aperm(nb, c(3,1,2))
    nb_r <- as.vector(nb_r)
    
    # Process #1.3 = Transition of live individuals between phenotypic categories (assumed fixed)
    NP <- P %*% nb_r
    
    # Rearrange back to initial configuration (stage within breeding value within phenotype)
    np <- array(NP, c(g,b,w))
    np_r <- aperm(np, c(3,2,1))
    np_r <- as.vector(np_r)
    NU_final <- np_r #% Should be the same as NU in the absence of specified transitions between breeding values and phenotypes.
    
    # Reproduction
    # Process #2.1 = Offspring production
    NF <- R %*% N[,i-1]
    
    # Rearrange for process #2.2
    nf <- array(NF,c(w,b,g))
    n <- aperm(nf[1, , ,  drop = TRUE], c(2,1)) # % g,b,1  This is a matrix of offspring produced from maternal phenotypes (rows) and breeding values (columns)
    nn <- apply(n, 2, sum)/ sum(n*binwidth_a)
    
    # Process #2.2 = Transmission of breeding  value
    Gle <- GVLE(0.5*(midpoint2_a), Va = Va)

    H1 <- conv2(nn, Gle * binwidth_a)  #identical with matlab
    H2  <- conv2(H1*binwidth_a, t(n)) # important to transpose n!
    bb2 <- seq(from=2*midpoint_a[1], to = 2*midpoint_a[length(midpoint_a)], by=(binwidth_a/2))
    
    H3 <- apply(H2, 2, function(x) interp1(bb2, x, midpoint_a))
    nf[1,,] <- H3 * sum(n) / sum(H3) # normalise to expected number of individuals
    
    
    # Rearrange for process #2.3
    NH <- as.vector(aperm(nf,c(3,1,2)))
    
    # Process # 2.3 = Attribution of offspring phenotype
    NM <- M%*%NH
    
    # Rearrange back to initial configuration (state within breeding value within phenotype)
    nm <- aperm(array(NM,c(g,w,b)), c(2,3,1))
    nm_r <- as.vector(nm)
    NF_final <- nm_r
    
    # Add all processes to get total population size
    N[,i] <- NF_final+NU_final
  } # end   for (i in 2:nit)
  
  popsize <- apply(N, 2, sum)
  
  NN <- array(N,c(w,b,g,nit)) # pop structured by stage * breeding values; phenotype, time
  ebvs <- apply(NN, c(2,4), sum)
  
  M_EBVS <- apply(ebvs,2, function(x) sum(x*midpoint_a)/sum(x))
  
  sumebvs <- apply(ebvs,2, sum)
  m_ebvs2 <- t(apply(ebvs, 1, function(x) x / sumebvs))
  
  ebvsJ <- aperm(array(apply(NN[1, , , ],c(1,3),sum), c(1,b,1,nit)),c(2,4,1,3))
  M_EBVSJ <- apply(ebvsJ, 2, function(x) sum(x*midpoint_a)/sum(x))
  sumebvsJ <- apply(ebvsJ,2, sum)
  m_ebvsJ2 <- t(apply(ebvsJ, 1, function(x) x / sumebvsJ)) #distribution breeding value
  
  ebvsA <- aperm(array(apply(NN[2, , , ],c(1,3),sum), c(1,b,1,nit)),c(2,4,1,3))
  M_EBVSA <- apply(ebvsA, 2, function(x) sum(x*midpoint_a)/sum(x))# mean breeding value
  sumebvsA <- apply(ebvsA,2, sum)
  m_ebvsA2 <- t(apply(ebvsA, 1, function(x) x / sumebvsA)) #distribution breeding value
  
  V <- mean(M_EBVS[2:timeA]-M_EBVS[1:timeA-1]) # Ebvs   = cumsum(ebvs*binwidth_a,1);
  
  pheno = apply(NN, c(3,4), sum) #frequency of phenotypic classes through time
  
  f_Z <- apply(pheno, 2, function(x) sum(x*f_pheno)/sum(x)) # average vital rate through time
  SJ_Z <- apply(pheno, 2, function(x) sum(x*SJ_pheno)/sum(x))
  SA_Z <- apply(pheno, 2, function(x) sum(x*SA_pheno)/sum(x))
  Y_Z <- apply(pheno, 2, function(x) sum(x*Y_pheno)/sum(x))
  Fe_Z <- f_Z * SJ_Z
  
  lambda_Z <- vector(length = length(f_Z))
  T_GEN_Z <- vector(length = length(f_Z))
  age_juv_Z <- vector(length = length(f_Z)) 
  for (i in 1:length(f_Z))
  {
    A <- matrix(data=c(SJ_Z[i]*(1-Y_Z[i]), Fe_Z[i], SJ_Z[i]*Y_Z[i], SA_Z[i]), nrow=2, byrow = TRUE) #Leslie matrix  
    eigenA <- eigen(A)
    lambda_Z[i] <- eigenA$values[1]
    
    w_vec_mean <- eigenA$vector[,1] / sum(eigenA$vector[,1]) # stable age-stage distribution, w = right eigenvector
    age_juv_Z[i] <- w_vec_mean[1]
    
    v_vec_mean <- eigen(t(A))$vector[,1]  / sum(eigen(t(A))$vector[,1])
    S_lambda <- (v_vec_mean %*% t(w_vec_mean)) / as.vector(t(v_vec_mean)%*%w_vec_mean)
    T_GEN_Z[i] <- lambda_Z[i]*(t(v_vec_mean) %*% w_vec_mean)/(t(v_vec_mean) %*% matrix(c(0,Fe,0,0), nrow=2, byrow = TRUE) %*% w_vec_mean)
  }
  
  MZ <- apply(pheno, 2, function(x) sum(x*midpoint_e)/sum(x)) # average phenotype through time
  Va_hat <- apply(ebvs * outer(X = midpoint_a, Y = M_EBVS, FUN = "-")^2, 2, function(x) sum(x)) / apply(ebvs, 2, sum)
  Vp_hat <- apply(pheno * outer(X = midpoint_e, Y = MZ, FUN = "-")^2, 2, function(x) sum(x)) / apply(pheno, 2, sum)
  
  properties_per_year <- data.frame(year = 1:length(f_Z), 
                                    Fe=Fe_Z,
                                    f=f_Z,
                                    SA=SA_Z,
                                    SJ=SJ_Z,
                                    Y=Y_Z,
                                    Z=MZ,
                                    Va=Va_hat,
                                    Vp=Vp_hat,
                                    M_EBVSA=M_EBVSA,
                                    M_EBVSJ=M_EBVSJ,
                                    M_EBVS=M_EBVS, 
                                    T_GEN=T_GEN_Z,
                                    lambda=lambda_Z,
                                    age_juv=age_juv_Z,
                                    N=popsize)
  
  colnames(pheno) <- 1:ncol(pheno)
  pheno <- as_tibble(pheno) 
  pheno$phenotypic_class <- 1:nrow(pheno)
  pheno$phenotypic_midvalue <- midpoint_e
  pheno <- pheno %>% pivot_longer(cols = 1:(ncol(pheno)-2), names_to = "year", values_to = "density", names_transform = list(year = as.integer))
  
  output <- list(T_GEN=T_GEN, V_ADAPT=V_ADAPT, LAMBDA0=LAMBDA0, V=V, 
                 properties_per_year=properties_per_year,
                 pheno_distribution=pheno, 
                 midpoint_phenotypic_class=midpoint_e)
  
  return(output)
}#end MPM()



#' MPM
#' 
#' Run a eco-evo deterministic matrix model
#' 
#' @param Beta numeric, strength of selection on the scale of linear predictors
#' 
#' @param h2 numeric on [0,1], heritability of the trait under selection
#' 
#' @param xp numeric, the range of phenotypic values considered
#' 
#' @param timeA numeric, number of time steps of adaptation. The simulations has two more time steps.
#' 
#' @param g number of phenotypic and genetic classes considered
#' 
#' @param Vp phenotypic variance
#' 
#' @param nind number of individuals to simulate, purely a matter of scale as there is no stochasticity
#' 
#' @param  Fe Fertility
#'
#' @param SA Adult survival rate
#' 
#' @param SJ Juvenile survival rate
#' 
#' @param S0 New born survival rate
#' 
#' @param Y maturation rate
#' 
#' @param selectedrate the vital rate under selection
#' 
#' @param density_tolerance Threshold to trigger warnings about the distribution of breeding values being too narrow or too wide
#' 
#' 
#' @export

MPM_archive <- function(Beta = 0.15, h2 = 0.2, xp=4 ,timeA=100, g = 20, Vp = 1,  nind = 100, 
                Fe=0.2286, SA=0.95, SJ=0.8,  S0=0.8, Y=0.07, selectedrate="Fe", selectmodel="log", 
                census = "prebreeding", density_tolerance = 10^(-5)){
  require(pracma) # for interpolation
  require(gsignal) # for 2D convolution
  require(purrr) # for reformatting output
  require(tidyr) # for reformatting output
  
  if(h2<=0 | h2>=1) {stop("h2 must be strictly greater than 0 and strictly smaller than 1")}
  
  w = 2 # number of stages
  m = 3 # three dimensions: stage, breeding value, phenotype
  nvital = 4 # number of vital rates (juvenile survival, maturation, adult survival, fertility)
  Va = h2*Vp # Additive genetic variance
  Ve = Vp - Va # Get the environmental variance (Vp = Va + Ve)
  
  # Create the phenotypic bins for the matrix
  min_pheno    = -xp 
  max_pheno    = xp 
  min_breeding = -xp 
  max_breeding = xp 
  mean_pheno   = 0 
  
  # Get the categories of phenotypes
  edges_e    = seq(from=min_pheno, to=max_pheno, length.out=g+1)
  binwidth_e = edges_e[2] - edges_e[1] # Bin width
  midpoint_e = edges_e[1:(length(edges_e)-1)] +  (binwidth_e/2) # Mid phenotype in each bin
  
  # Get the categories of breeding values
  edges_a    = edges_e[(edges_e<=max_breeding)&(edges_e>=min_breeding)]
  binwidth_a = edges_a[2] - edges_a[1] # Bin width
  midpoint_a = edges_a[1:(length(edges_a)-1)] +  (binwidth_a/2) # Mid phenotype in each bin
  
  b = length(midpoint_a) # this is the number of classes for breeding value
  midpoint2_a = seq(2*midpoint_a[1], 2*midpoint_a[length(midpoint_a)], by=binwidth_a)
  
  Gle <- GVLE(0.5*(midpoint2_a), Va = Va)
  if(Gle[1]>density_tolerance | Gle[length(Gle)] > density_tolerance) {warning("The initial distribution of breeding values reaches the boundaries of the possible range. You should probably decrease h2, decrease Vp, or increase xp.")}
  if((sum(Gle) - (1/(0.5 * binwidth_a))) > sqrt(density_tolerance)) {warning("The initial distribution of breeding values is too narrow given the coarsness of phenotypic classes; genetic change may not be computed accurately. You should probably increase h2, increase Vp, increase g or decrease xp. ")}
  
  # B matrix
  Bik = array(data=0, dim = c(b,b,w,g)) # initiate empty cell array
  Bik[1:g,1:g, , ] = diag(b)
  B = BD_proj_mat(Bik);
  
  # P matrix
  Pij = array(data=0, dim=c(g,g,w,b)) # initiate empty cell array
  Pij[1:g,1:g, , ] = diag(g)
  P = BD_proj_mat(Pij);
  
  # Simulate an initial dataframe for the vector of population size at t initial
  nind    = 100 # number of individuals to simulate
  
  N_ve = t(exp((-outer(midpoint_e, midpoint_a, FUN="-")^2 / (2*Ve) )) / colSums(exp(-outer(midpoint_e, midpoint_a, FUN="-")^2 / (2*Ve))))
  
  nbv  = exp(-midpoint_a^2/(2*Va))/sum(exp(-midpoint_a^2/(2*Va))) # Get the proportion of individuals within each class of breeding value
  Nbv  = nind*nbv # Multiply by the total number of individuals in the population to get the number of individuals in each class of breeding value
  nz   = N_ve %*% nbv * binwidth_a
  Nz   = nind*nz
  NN   = t(t(N_ve) * Nbv)
  
  ninit = aperm(array(data = NN/w, dim = c(g,g,w)), perm=c(3,2,1))
  Ninit = t(t(as.vector(ninit)))
  
  
  nit = timeA + 2 # number of years
  
  M_EBVS = vector(length = nit)
  M_EBVSJ = vector(length = nit)
  M_EBVSA = vector(length = nit)
  f_Z   = vector(length = nit)
  SJ_Z   = vector(length = nit)
  SA_Z   = vector(length = nit)
  Y_Z  = vector(length = nit)
  Va_hat = vector(length = nit)
  Vp_hat = vector(length = nit)
  
  A <- matrix(data=c(SJ*(1-Y), Fe, SJ*Y, SA), nrow=2, byrow = TRUE) #Leslie matrix
  f <- Fe/SJ
  eigenA <- eigen(A)
  LAMBDA0 <- lambda1_mean <- eigenA$values[1]
  # stable age-stage distribution, w = right eigenvector
  w_vec_mean <- eigenA$vector[,1] / sum(eigenA$vector[,1])
  v_vec_mean <- eigen(t(A))$vector[,1]  / sum(eigen(t(A))$vector[,1])
  S_lambda <- (v_vec_mean %*% t(w_vec_mean)) / as.vector(t(v_vec_mean)%*%w_vec_mean)
  T_GEN <- lambda1_mean*(t(v_vec_mean) %*% w_vec_mean)/(t(v_vec_mean) %*% matrix(c(0,Fe,0,0), nrow=2, byrow = TRUE) %*% w_vec_mean)
  
  # Get Theta (parameters)
  S0_pheno <- rep(S0, times=b)
  SJ_pheno <- rep(SJ, times=b)
  Y_pheno  <- rep(Y, times=b)
  SA_pheno <- rep(SA, times=b)
  f_pheno  <- rep(f, times=b)
  #F_pheno  <- f_pheno*SA_pheno; #POST
  F_pheno  <- f_pheno*SJ_pheno;
  
  if (selectedrate=="Fe"){
    if (selectmodel=="log"){
      # Poisson
      beta <- Beta
      f_pheno  <- exp(log(f) + beta*midpoint_e)
      dS <- Beta*f
    }else{
      if (selectmodel =="logit"){
        beta <- Beta;
        if( (f<=0) | (f>=1)) {stop("When using a logit model for selection the fertility rate (Fe/Sj) must be between 0 and 1 (excluded)")}
        f_pheno  <- 1/(1+exp(-(logit(f) + beta*midpoint_e)))
        dS <- Beta*f*(1-f)
      }
    }
    
    S_LAMBDA <- SJ*S_lambda[1,2]
    F_pheno  <- f_pheno*SJ_pheno #PRE
  }else{
    if(selectedrate=="SJ"){
      beta <- Beta
      SJ_pheno <- 1/(1+exp(-(logit(SJ) + beta*midpoint_e)))
      S_LAMBDA <- (1-Y)*S_lambda[1,1]+Y*S_lambda[2,1]+f*S_lambda[1,2]
      dS <- Beta*SJ*(1-SJ)
      F_pheno <- f_pheno*SJ_pheno
    }else{
      if(selectedrate=="SA"){
        beta <- Beta
        SA_pheno <- 1/(1+exp(-(logit(SA) + beta*midpoint_e)))
        #%F_pheno <- f_pheno*SA_pheno #POST
        #S_LAMBDA <- S_lambda[2,2]+f*S_lambda[1,2] # POST
        S_LAMBDA <- S_lambda[2,2] # PRE
        dS <- Beta*SA*(1-SA)
      }else{
        if(selectedrate=="Y"){
          beta <- Beta
          Y_pheno <- 1/(1+exp(-(logit(Y) + beta*midpoint_e)))
          S_LAMBDA <- SJ*(S_lambda[2,1]-S_lambda[1,1])
          dS <- Beta*Y*(1-Y)
        }else{stop("selectedrate not recognized")}
      }
    }
  }
  
  V_ADAPT <- S_LAMBDA*dS*Va
  
  # Get the matrices U, B, P for Utilde and R, H, M for Ftilde
  # U matrix
  Ujk <- array(data=0, dim=c(w,w,b,g)) # initiate empty cell array
  Ujk[1,1,,] <- matrix(SJ_pheno*(1-Y_pheno), nrow = b, ncol = g, byrow = FALSE)
  Ujk[2,1,,] <- matrix(SJ_pheno*Y_pheno, nrow = b, ncol = g, byrow = FALSE) 
  Ujk[2,2,,] <- matrix(SA_pheno, nrow = b, ncol = g, byrow = FALSE) 
  U = BD_proj_mat(Ujk)
  
  #%% R matrix
  Rjk <- array(data=0, dim=c(w,w,b,g)) 
  Rjk[1,2,,] <- matrix(F_pheno, nrow = b, ncol = g, byrow = FALSE) 
  R <- BD_proj_mat(Rjk)
  
  #%% H matrix
  #% To construct the H matrix, we need the Ga matrix
  Vle <- Va
  
  
  #%% M matrix
  #% To construct the M matrix, we need the Ge matrix
  #% For each breeding value, what is the potential distribution of offspring phenotype? - Here we add noise (Ve) to the breeding value to obtain offspring phenotype
  emat <- t(exp(-(outer(midpoint_e, midpoint_a, FUN="-"))^2 / (2*Ve) ) * binwidth_a / apply(exp(-(outer(midpoint_e, midpoint_a, FUN="-"))^2/(2*Ve))*binwidth_e, 1, sum))
  
  # The matrix M, which is the transmission of maternal phenotype to offspring phenotype goes through the breeding values only
  Mij <- array(data=0, dim=c(g,g,w,b))
  mij <- kronecker(array(data=1, dim = c(1,1,g)), emat)
  mij <- aperm(mij, c(1,3,2))
  Mij[,,1,] <- mij
  M <- BD_proj_mat(Mij)
  
  #%% Compute N(t)
  N <- matrix(0, nrow = w*b*g, ncol=nit)
  N[,1] <- Ninit
  NN  <- array(0, c(w,b,g,nit))
  NN[ , , ,1] <- ninit
  
  for (i in 2:nit)
  {
    # Live individual transitions
    # Process #1.1 = Transition of live individuals between states
    NU <- U %*% N[,i-1]
    
    # Rearrange for process #1.2
    nu <- array(NU,c(w,b,g))
    nu_r <- aperm(nu,c(2,1,3)) # bwg
    nu_r <- as.vector(nu_r)
    
    # Process #1.2 = Transition of live individuals between breeding value categories (assumed fixed)
    NB <- B %*% nu_r
    
    # Rearrange for process #1.3
    nb <- array(NB, c(b,w,g))
    nb_r <- aperm(nb, c(3,1,2))
    nb_r <- as.vector(nb_r)
    
    # Process #1.3 = Transition of live individuals between phenotypic categories (assumed fixed)
    NP <- P %*% nb_r
    
    # Rearrange back to initial configuration (stage within breeding value within phenotype)
    np <- array(NP, c(g,b,w))
    np_r <- aperm(np, c(3,2,1))
    np_r <- as.vector(np_r)
    NU_final <- np_r #% Should be the same as NU in the absence of specified transitions between breeding values and phenotypes.
    
    # Reproduction
    # Process #2.1 = Offspring production
    NF <- R %*% N[,i-1]
    
    # Rearrange for process #2.2
    nf <- array(NF,c(w,b,g))
    n <- aperm(nf[1, , ,  drop = TRUE], c(2,1)) # % g,b,1  This is a matrix of offspring produced from maternal phenotypes (rows) and breeding values (columns)
    nn <- apply(n, 2, sum)/ sum(n*binwidth_a)
    
    # Process #2.2 = Transmission of breeding  value
    Gle <- GVLE(0.5*(midpoint2_a), Va = Va)
    Gle <- Gle / (sum(Gle) * 0.5*binwidth_a)  # normalise density
    
    H1 <- conv2(nn, Gle * binwidth_a)  #identical with matlab
    H2  <- conv2(H1*binwidth_a, t(n)) # important to transpose n!
    bb2 <- seq(from=2*midpoint_a[1], to = 2*midpoint_a[length(midpoint_a)], by=(binwidth_a/2))
    
    H3 <- apply(H2, 2, function(x) interp1(bb2, x, midpoint_a))
    nf[1,,] <- H3
    
    # Rearrange for process #2.3
    NH <- as.vector(aperm(nf,c(3,1,2)))
    
    # Process # 2.3 = Attribution of offspring phenotype
    NM <- M%*%NH
    
    # Rearrange back to initial configuration (state within breeding value within phenotype)
    nm <- aperm(array(NM,c(g,w,b)), c(2,3,1))
    nm_r <- as.vector(nm)
    NF_final <- nm_r
    
    # Add all processes to get total population size
    N[,i] <- NF_final+NU_final
  } # end   for (i in 2:nit)
  
  popsize <- apply(N, 2, sum)
  
  NN <- array(N,c(w,b,g,nit)) # pop structured by stage * breeding values; phenotype, time
  ebvs <- apply(NN*binwidth_a, c(2,4), sum)
  M_EBVS <- apply(ebvs,2, function(x) sum(x*midpoint_a)/sum(x))
  
  sumebvs <- apply(ebvs*binwidth_a,2, sum)
  m_ebvs2 <- t(apply(ebvs, 1, function(x) x / sumebvs))
  
  ebvsJ <- aperm(array(apply(NN[1, , , ]*binwidth_a,c(1,3),sum), c(1,b,1,nit)),c(2,4,1,3))
  M_EBVSJ <- apply(ebvsJ, 2, function(x) sum(x*midpoint_a)/sum(x))
  sumebvsJ <- apply(ebvsJ*binwidth_a,2, sum)
  m_ebvsJ2 <- t(apply(ebvsJ, 1, function(x) x / sumebvsJ)) #distribution breeding value
  
  ebvsA <- aperm(array(apply(NN[2, , , ]*binwidth_a,c(1,3),sum), c(1,b,1,nit)),c(2,4,1,3))
  M_EBVSA <- apply(ebvsA, 2, function(x) sum(x*midpoint_a)/sum(x))# mean breeding value
  sumebvsA <- apply(ebvsA*binwidth_a,2, sum)
  m_ebvsA2 <- t(apply(ebvsA, 1, function(x) x / sumebvsA)) #distribution breeding value
  
  V <- mean(M_EBVS[2:timeA]-M_EBVS[1:timeA-1]) # Ebvs   = cumsum(ebvs*binwidth_a,1);
  
  pheno = apply(NN, c(3,4), sum) #frequency of phenotypic classes through time
  
  f_Z <- apply(pheno, 2, function(x) sum(x*f_pheno)/sum(x)) # average vital rate through time
  SJ_Z <- apply(pheno, 2, function(x) sum(x*SJ_pheno)/sum(x))
  SA_Z <- apply(pheno, 2, function(x) sum(x*SA_pheno)/sum(x))
  Y_Z <- apply(pheno, 2, function(x) sum(x*Y_pheno)/sum(x))
  Fe_Z <- f_Z * SJ_Z
  
  lambda_Z <- vector(length = length(f_Z))
  T_GEN_Z <- vector(length = length(f_Z))
  age_juv_Z <- vector(length = length(f_Z)) 
  for (i in 1:length(f_Z))
  {
    A <- matrix(data=c(SJ_Z[i]*(1-Y_Z[i]), Fe_Z[i], SJ_Z[i]*Y_Z[i], SA_Z[i]), nrow=2, byrow = TRUE) #Leslie matrix  
    eigenA <- eigen(A)
    lambda_Z[i] <- eigenA$values[1]
    
    w_vec_mean <- eigenA$vector[,1] / sum(eigenA$vector[,1]) # stable age-stage distribution, w = right eigenvector
    age_juv_Z[i] <- w_vec_mean[1]
    
    v_vec_mean <- eigen(t(A))$vector[,1]  / sum(eigen(t(A))$vector[,1])
    S_lambda <- (v_vec_mean %*% t(w_vec_mean)) / as.vector(t(v_vec_mean)%*%w_vec_mean)
    T_GEN_Z[i] <- lambda_Z[i]*(t(v_vec_mean) %*% w_vec_mean)/(t(v_vec_mean) %*% matrix(c(0,Fe,0,0), nrow=2, byrow = TRUE) %*% w_vec_mean)
  }
  
  MZ <- apply(pheno, 2, function(x) sum(x*midpoint_e)/sum(x)) # average phenotype through time
  Va_hat <- apply(ebvs * outer(X = midpoint_a, Y = M_EBVS, FUN = "-")^2, 2, function(x) sum(x)) / apply(ebvs, 2, sum)
  Vp_hat <- apply(pheno * outer(X = midpoint_e, Y = MZ, FUN = "-")^2, 2, function(x) sum(x)) / apply(pheno, 2, sum)
  
  properties_per_year <- data.frame(year = 1:length(f_Z), 
                                    Fe=Fe_Z,
                                    f=f_Z,
                                    SA=SA_Z,
                                    SJ=SJ_Z,
                                    Y=Y_Z,
                                    Z=MZ,
                                    Va=Va_hat,
                                    Vp=Vp_hat,
                                    M_EBVSA=M_EBVSA,
                                    M_EBVSJ=M_EBVSJ,
                                    M_EBVS=M_EBVS, 
                                    T_GEN=T_GEN_Z,
                                    lambda=lambda_Z,
                                    age_juv=age_juv_Z,
                                    N=popsize)
  
  colnames(pheno) <- 1:ncol(pheno)
  pheno <- as_tibble(pheno) 
  pheno$phenotypic_class <- 1:nrow(pheno)
  pheno$phenotypic_midvalue <- midpoint_e
  pheno <- pheno %>% pivot_longer(cols = 1:(ncol(pheno)-2), names_to = "year", values_to = "density", names_transform = list(year = as.integer))
  
  output <- list(T_GEN=T_GEN, V_ADAPT=V_ADAPT, LAMBDA0=LAMBDA0, V=V, 
                 properties_per_year=properties_per_year,
                 pheno_distribution=pheno, 
                 midpoint_phenotypic_class=midpoint_e)
  
  return(output)
}#end MPM()


