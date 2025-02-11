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
#' Created a block-diagonal matrix from an array
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
#' 
#' @export

MPM <- function(Beta = 0.15, h2 = 0.2, xp=4 ,timeA=100, g = 40, Vp = 1,  nind = 100, 
                Fe=1, SA=0.3, SJ=0.1,  S0=0.1, Y=0.8, selectedrate="Fe", selectmodel="logit", 
                census = "prebreeding"){
  require(pracma) # for interpolation
  require(gsignal) # for 2D convolution
  
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
  
  N_ve = t(exp((-outer(midpoint_e, midpoint_a, FUN="-")^2 / (2*Ve) )) / colSums(exp(-outer(midpoint_e, midpoint_a, FUN="-")^2 / (2*Ve)) * binwidth_e))
  
  nbv  = exp(-midpoint_a^2/(2*Va))/sum(exp(-midpoint_a^2/(2*Va))*binwidth_a) # Get the proportion of individuals within each class of breeding value
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
    Gle <- GVLE(0.5*(midpoint2_a))
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
  
  pheno = apply(NN*binwidth_e, c(3,4), sum) #frequency of phenotypic classes through time
  
  f_Z <- apply(pheno, 2, function(x) sum(x*f_pheno)/sum(x)) # average vital rate through time
  SJ_Z <- apply(pheno, 2, function(x) sum(x*SJ_pheno)/sum(x))
  SA_Z <- apply(pheno, 2, function(x) sum(x*SA_pheno)/sum(x))
  Y_Z <- apply(pheno, 2, function(x) sum(x*Y_pheno)/sum(x))
  MZ <- apply(pheno, 2, function(x) sum(x*midpoint_e)/sum(x))
  Va_hat <- apply(ebvs * outer(X = midpoint_a, Y = M_EBVS, FUN = "-")^2, 2, function(x) sum(x)) / apply(ebvs, 2, sum)
  Vp_hat <- apply(pheno * outer(X = midpoint_e, Y = MZ, FUN = "-")^2, 2, function(x) sum(x)) / apply(pheno, 2, sum)

  output <- list(T_GEN=T_GEN, V_ADAPT=V_ADAPT, LAMBDA0=LAMBDA0, V=V, 
                 Vp_hat=Vp_hat, Va_hat=Va_hat, Y_Z=Y_Z, SA_Z=SA_Z, SJ_Z=SJ_Z, f_Z=f_Z, MZ=MZ, 
                 M_EBVSA=M_EBVSA, M_EBVSJ=M_EBVSJ, M_EBVS=M_EBVS) 
  
  return(output)
}#end MPM()


