# Fix the typo on Line 670 and 745, have not tested the performance yet
# This is the version that runs on GPU.
# We set a seed at the very beginning and then generate the data. 
# The jobID will be used as the seed.
walltime = 4 # time limit of the job
time_begin = Sys.time()
library(mvtnorm)
library(tmvtnorm)
library(sde)
library(fChange)
library(invgamma)
library(reshape2)
library(ggplot2)
library(scpt)
library(splines)
library(matrixcalc)

# \phi in the data generation, which is the scalar of spatial distance
phis_real = 2
# \rho in the data generation, which controls the strength of change signal
signal = 1
# if plots of locations, a, b, functional data, CUSUM for all locations
# need to be generated and saved, set pic to TRUE.
pic = FALSE 
# if this is the first time to generate the fake data, set first to TRUE. 
# If not, will read from the directory specified directly.
first = TRUE

iters = 20000 # number of the iterations for MCMC
D = 21 # L, number of Fourier basis functions used to generate the functional data
ns = 50 # N, number of spatial locations in the alternative region
ns_null = 5 # N_0, number of loctaions without changepoint (true null)
Nt = 50 # T, number of curves for each location
M = 365 # number of observations on each curve, only used for Gromenko et al. (2017)
thin_step = 10 # stepsize when thinning
job<-as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) # get the job id and set as the seed

basis = create.fourier.basis(rangeval = c(0, 1), nbasis = D)
ns_alt = ns - ns_null # N_a, number of true alternative locations in the alternative region
nt = Nt-1 # number of time points that are fit in the model, i.e. 1,..,T-1.
n=ns*nt # N(T-1), dimension of the vector e in Stage I

cat("Settings:phi=",phis_real, "signal=",signal,"\n")
cat("Seed", job, "\n")
# Specify the directory of storing the results of accuracy, coverage, and length of intervals
summary_dir = paste("~/summary")
# Specify the directory of storing the details of each simulation
main_dir = paste("~/", phis_real, "_", signal, "_", job, sep="")
read_dir = main_dir
result_dir = main_dir
dir.create(file.path(result_dir), showWarnings = FALSE) # create folder if not exist
setwd(file.path(result_dir))

set.seed(913)
################################## generate locations ################################## 
xx <- runif(ns, min=0, max=10) # x coordinates of the locations
yy <- runif(ns, min=0, max=10) # y coordinates of the locations
location = cbind(xx,yy) # coordiates of the locations

ds = as.matrix(dist(location))
cor_matrix = exp(-ds/phis_real)

null_loc0 = sample(1:ns, 1)
null_dist = ds[, null_loc0]
index_null = order(null_dist)[1:ns_null] # A cluster of locations that are considered as null locations
index_alt = (1:ns)[-index_null]

cor_matrix_alt = cor_matrix[-index_null, -index_null]
if(first){
  ################################## Generate changepoint ################################## 
  k_star_o =  rtmvnorm(1, mean = rep(0.5, ns_alt),
                       sigma = 1*cor_matrix_alt,
                       lower=rep(0.15, length = ns_alt),
                       upper=rep(0.85, length = ns_alt),
                       algorithm="gibbs", burn.in.samples=100)
  
  k_star = rep(1, ns) # Changepoint for all locations. For those without changepoint, k_star=1.
  k_star[index_alt] = round(k_star_o*Nt)/Nt
  # Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('red','blue'))
  loc_col <- rbPal(10)[as.numeric(cut(k_star[index_alt],breaks = 10))]
  if(pic){
    # Plot the changepoints
    #setwd(result_dir)
    png("locations.png")
    plot(xx[index_alt],yy[index_alt],pch = 16,col = loc_col, xlab = "x", ylab = "y", xlim = c(0, 10.3), ylim = c(0, 10.3), 
         main = "Changepoint Plot")
    points(xx[index_null], yy[index_null])
    text(xx+0.3, yy+0.3, 1:ns, cex = 0.6)
    dev.off()
  }
  ################################## Generate change functions ################################## 
  change_f = matrix(0,nrow = D,ncol = ns) # The coefficients of change functions. 
  D0 = (D-1)/2 # The number of pairs of sin and cos basis functions needed
  Ds = rep(1:D0, each = 2)
  sigma_delay = c(1, 1/Ds^3)  # variance of \eta_l in the generation of change function
  mean_delay = signal*c(1, 1/Ds^2)  # mean of \eta_l in the generation of change function
  
  # Generate the coefficients of change functions for all true alternative locations
  for(j in 1:D){
    change_f[j, index_alt] = rmvnorm(1,mean = rep(mean_delay[j], ns_alt), sigma = 0.1*sigma_delay[j] * cor_matrix_alt)
  }
  if(pic){
    # Plot of the change functions at all locations
    png("changefunctions.png")
    for(i in 1:ns){
      if(i ==1){plot(fd(change_f[,i], basis), ylim = c(-4,6))
      }else{
        lines(fd(change_f[,i], basis), col = i)
      }
    }
    dev.off()
    
    # Plot of the absolute value of the coefficients of change functions at all locations
    png("changefunc_abscoef.png")
    for(i in 1:ns){
      if(i ==1){
        plot(1:D, abs(change_f[,i]), main = "abs coefs of change functions",type = "l", xlab = "basis", ylab= "coefficient", ylim = c(0,max(change_f)))
      }else{
        lines(1:D, abs(change_f[,i]), col = i, type = "l")
      }
    }
    dev.off()
  }
  ################################## Generate functional data ##################################
  # During the genetaion, get the Y_{T,k}(s) process and reasonable initial values
  #setwd(result_dir)
  S_sim = matrix(0, ncol = ns, nrow = Nt+1)
  a_expect = b_expect = k_expect = a0_expect = b0_expect = beta_expect = beta_old = rep(0, ns)
  x0 = (0:Nt)/Nt # scaled time points, i.e. 0, 1/T, 2/T, ... (T-1)/T, 1
  
  funcdatpoint = function(x, D, coef){
    # x is the x value interested, and in this function, we will get the data point on the
    # functional data with D Fourier basis and coef as the coeffiences
    result = 0
    for(i in 1:D){
      if(i == 1) result = result + coef[i]
      else {if(i %% 2 == 0) result = result + coef[i] * sqrt(2)  * sin(floor(i/2)*2*pi*x)
      else result = result + coef[i] * sqrt(2)  * cos(floor(i/2)*2*pi*x)}
    }
    return(result)
  }
  # a simple example of illustrating how the function works
  # change.coef = 1:D/10
  # plot(fd(change.coef, basis))
  # lines(seq(0, 1, length = 365), funcdatpoint(seq(0, 1, length = 365), 21, change.coef), col = "red")
  approxFourier = function(dat, pts = 50) {
    t = seq(0, 1, length.out = pts)
    fmat = matrix(0, pts, ncol(dat))
    for(f in 1:ncol(dat)) {
      fmat[,f] = fmat[,f] + dat[1, f]
      
      for(i in 1:((nrow(dat)-1)/2)) {
        fmat[,f] = fmat[,f] + dat[2*i, f]*sin(i*2*pi*t)*sqrt(2)
      }
      for(i in 1:((nrow(dat)-1)/2)) {
        fmat[,f] = fmat[,f] + dat[2*i + 1, f]*cos(i*2*pi*t)*sqrt(2)
      }
    }
    fmat
  }
  
  # \xi_l, Generate the coefficients for error functions at all locations.
  function.coef = array(0, dim = c(D,Nt,ns))
  for(i in 1:D){
    for(j in 1:Nt){
      fiid.coef = rmvnorm(1,mean = rep(0,ns),sigma = 0.5 * sigma_delay[i] * cor_matrix)
      for(k in 1:ns){
        function.coef[i,j,k] = fiid.coef[k]
      }
    }
  }  
  
  #####################################################################
  ################### GKR method Gromenko et al. (2017)################
  #####################################################################
  DATA3 <- array(0, dim=c(M, ns, Nt))
  simpson_r<-function(series){
    len<-length(series);
    sum<-0;
    sum<-.C("simpson", as.double(series), as.double(len-1), as.double(sum))[[3]];
    return(sum);
  }
  for(i in 1:ns){
    # cat(i)
    change.coef = function.coef[,,i]
    if(i %in% index_alt){
      change.coef[, (Nt*k_star[i]+1):Nt] = change.coef[,(Nt*k_star[i]+1):Nt]+change_f[,i]}
    for(m in 1:Nt){
      DATA3[, i, m] = funcdatpoint(seq(0, 1, length = M), D, change.coef[, m])}
  }
  loc = location
  L_per_curve<-dim(DATA3)[1] # number of observations per curve
  N<-dim(DATA3)[2] # number of spatial locations
  NY<-dim(DATA3)[3] # number of curves per location
  L<-L_per_curve*NY # total temporal observations for each location
  DATA.DIFF3<-DATA3[,,2:NY]-DATA3[,,1:(NY-1)] # differenced series in (C.1)
  
  c.hat.pr <- .Call('get_svar',PACKAGE='scpt',as.matrix(loc*pi/180), DATA.DIFF3, dim(DATA.DIFF3))
  dd<-c.hat.pr[[1]]
  H.hat<-c.hat.pr[[3]]
  dd.v<-as.vector(dd)
  H.hat.v<-as.vector(H.hat)
  #removing zeroes
  H.hat.vn<-H.hat.v[dd.v!=0]
  dd.vn<-dd.v[dd.v!=0]
  
  #noparametric BSpline estimation
  fm1 <- lm(H.hat.vn ~ bs(dd.vn, df = 9))
  a<-quantile(dd.vn,.75)
  b<-quantile(dd.vn,.95)
  slope<-predict(fm1, data.frame(dd.vn = a))/(b-a)
  
  H.final<-matrix(1,N,N)
  sp.cov<-predict(fm1, data.frame(dd.vn = dd.vn))
  for(i in 1:length(dd.vn)){
    H.final[dd.vn[i]==dd]<-sp.cov[i]
    if(dd.vn[i] > a & dd.vn[i] < b){
      H.final[dd.vn[i]==dd] <- slope*(b-a)-(dd.vn[i]-a)*slope;
    }
    if(dd.vn[i] > b){
      H.final[dd.vn[i]==dd] <- 0;
    }
  }
  
  H.vec<-as.vector(H.final)
  H.vec.n<-H.vec[dd!=0]
  dd.n<-dd[dd!=0]
  cov.dt<-cbind(dd.n, H.vec.n)
  cov.dt<-cov.dt[order(cov.dt[,1]),]
  
  #caclulating weights
  c.hm<-matrix(0,N,N)
  diag(c.hm)<-c.hat.pr[[2]]
  rho.hat.f<-c.hm%*%H.final%*%c.hm
  one<-matrix(1, ncol=1, nrow=N)
  w<-solve(rho.hat.f^2)%*%one
  w<-w/sum(w[,1])
  
  #temporal covariance
  tvar<-.Call('get_tvar', PACKAGE='scpt',w[,1], c.hat.pr[[2]], DATA.DIFF3, dim(DATA.DIFF3))
  cat("The temporal covariance is positive definite", is.positive.definite(tvar[[1]]))
  FPC<-eigen(tvar[[1]], symmetric=T, only.values = FALSE, EISPACK = FALSE)
  val<-FPC$values[FPC$values>0]
  val<-val/sum(val)
  
  cvar<-val[1]
  for(i in 2:length(val))
    cvar<-c(cvar, cvar[i-1]+val[i])
  
  #NFPC<-length(val[cvar<.78])
  NFPC<-length(val)
  v.hat<-c()
  for(i in 1:NFPC)
    v.hat<-cbind(v.hat, FPC$vectors[,i]/sqrt(simpson_r(FPC$vectors[,i]^2)))
  
  rho_e<-eigen(rho.hat.f)
  rho.hat.f.clip<-(rho_e$vectors[,rho_e$values>0])%*%diag(rho_e$values[rho_e$values>0])%*%t(rho_e$vectors[,rho_e$values>0])
  suppressWarnings(Lb<-matrix(t(chol(rho.hat.f.clip,pivot=TRUE)),nrow=N))
  
  nnn<-NFPC
  aaa<-.Call('get_stat',PACKAGE='scpt',DATA3, dim(DATA3), as.matrix(v.hat[,1:nnn]), w[,1], val[1:nnn], nnn)
  r2<-aaa[[2]][-1,]
  if(pic){
    png("Gromenko's Method.png")
    plot(1:dim(r2)[2], r2[dim(r2)[1],], type = "l", xlab = "Number of curves", ylab = "Statistic curve")
    dev.off()
  }
  k_star_gro = which.max(r2[dim(r2)[1],])/dim(r2)[2]
  
  ##################################################################
  #################### FF Method Aue et al.(2018) ##################
  ##################################################################
  pvalue = change_est = trs = SNR = rep(0, ns)
  confs = matrix(0, nrow = ns, ncol =2)
  
  for(i in 1:ns){
    # cat(i)
    change.coef = function.coef[,,i]
    if(i %in% index_alt){
      change.coef[, (Nt*k_star[i]+1):Nt] = change.coef[,(Nt*k_star[i]+1):Nt]+change_f[,i]}
    fiid_c = fd(change.coef, basis)
    h = round(opt_bandwidth(fiid_c, "BT", "PR")$hat_h_opt)
    # Plot the functional data
    if(pic){
      png(paste("FuncDat", i, ".png", sep=""))
      if(i %in% index_alt){
        args = list(approxFourier(change.coef[, 1:(Nt*k_star[i])]), 
                    approxFourier(change.coef[, (Nt*k_star[i]+1):Nt]))
      }else{
        args = list(approxFourier(change.coef))
      }
      fmats = lapply(seq_along(args), function(x) cbind(melt(as.matrix(args[[x]])), x))
      fmats = do.call("rbind", fmats)
      fmats[["col"]] = as.factor(fmats[["x"]])
      fmats[["Var2"]] = paste0(fmats[["Var2"]], fmats[["col"]])
      
      plt = ggplot() +
        geom_line(data = fmats, aes(x = Var1,
                                    y = value,
                                    group = Var2,
                                    col = col),
                  alpha = 0.9,
                  size = 0.5) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(legend.position="none")
      
      print(plt)
      dev.off()
      
    }
    
    resultiid = change_FF(fiid_c,h=h)
    confint = Conf_int(fiid_c, h=h) # default level 0.05
    confs[i, 1] = confint[1]/Nt
    confs[i, 2] = confint[3]/Nt
    pvalue[i] = resultiid$pvalue
    change_est[i] = resultiid$change
    fdata = center.fd(fiid_c) ; basis = fdata$basis; samp = fdata$coefs;N = ncol(samp)
    S_sim[2,i] = sum((samp[, 1]- (1/N) * rowSums(samp[, 1:N]))^2)/N
    for (j in (2:N)) {
      S_sim[j+1, i] = sum((rowSums(samp[, 1:j]) - (j/N) * rowSums(samp[, 1:N]))^2)/N
    }
    k_star_guess = min(which(S_sim[,i] == max(S_sim[,i])))-1 # estimate from the FF method
    k_expect[i] = k_star_guess/Nt
    beta_y = max(S_sim[,i])
    beta_x = (which.max(S_sim[,i])-1)/Nt
    beta_old[i] = beta_y/(beta_x*(beta_x-1))
    
    LongRunC = LongRun(fdobj = fdata, h = h)
    lambda = LongRunC$e_val
    phi = LongRunC$e_fun$coefs
    dat.b = fiid_c[1:k_star_guess]
    dat.a = fiid_c[(k_star_guess + 1):N]
    mean.b = rowMeans(dat.b$coefs)
    mean.a = rowMeans(dat.a$coefs)
    delta_basis = mean.a - mean.b # estimate the change function
    
    a_expect[i] = 2 * sum(lambda^2)
    for(m in 1:D){
      b_expect[i] = b_expect[i]+4*lambda[m]*(sum(phi[,m]*delta_basis))^2
    }
    a0_expect[i] = sum(lambda)
    b0_expect[i] = sum(delta_basis^2)
    
    trs[i] = sum(diag(LongRunC$covm))
    SNR[i] = k_star[i]*(1-k_star[i])*b0_expect[i]/trs[i]
    m = beta_old[i]*((k_expect[i]-1)*x0+(x0-k_expect[i])*(ifelse(x0>k_expect[i],1,0)))
    v = a_expect[i]*x0^2*(1-x0)^2+b_expect[i]*Nt*(k_expect[i])^2*x0*(1-x0)^3*(ifelse(x0>k_expect[i],1,0))+
      b_expect[i]*Nt*(1-k_expect[i])^2*x0^3*(1-x0)*(ifelse(x0<=k_expect[i],1,0))
    betas = seq(-1, floor(beta_old[i]), length = (abs(floor(beta_old[i]))-1)*10+1)
    res = rep(0, length(betas))
    for(j in 1:length(betas)){
      means = betas[j]*((k_expect[i]-1)*x0+(x0-k_expect[i])*(ifelse(x0>k_expect[i],1,0)))
      res[j] = sum(((S_sim[,i]-means)[-c(1, length(means))]/sqrt(v[-c(1, length(means))]))^2)
    }
    beta_expect[i] = betas[which.min(res)]
    if(pic){
      png(paste("CUSUM", i,".png", sep=""))
      plot(0:N, S_sim[,i], xlab = "time", ylab = "CUSUM", type = "l", main = paste("Location", i, "with p-value", pvalue[i]))
      abline(v = k_star[i]*N, col = "red")
      dev.off()
    }
  }
  
  # cat("p-values(derived from change_FF function) when h = Nt^", alpha, "\n")
  cat(pvalue, "\n")
  
  a_col <- rbPal(10)[as.numeric(cut(a_expect,breaks = 10))]
  if(pic){
    png("a Plot.png")
    plot(xx,yy,pch = 16,col = a_col, xlab = "x", ylab = "y", xlim = c(0, 10.3), ylim = c(0, 10.3), 
         main = "Estimated a Plot")
    text(xx+0.3, yy+0.3, round(a_expect, 2), cex = 0.6)
    dev.off()
    
    png("b Plot.png")
    b_col <- rbPal(10)[as.numeric(cut(b_expect,breaks = 10))]
    plot(xx,yy,pch = 16,col = b_col, xlab = "x", ylab = "y", xlim = c(0, 10.3), ylim = c(0, 10.3), 
         main = "Estimated b Plot")
    text(xx+0.3, yy+0.3, round(b_expect, 2), cex = 0.6)
    dev.off()
    
    cat("Correlation between estimates of a and b:", cor(a_expect, b_expect),"\n")
  }
  save(S_sim, file = "S_sim.RData")
  save(k_star, file = "k_star.RData")
  save(index_alt, file = "index_alt.RData")
  save(beta_expect, file = "beta_expect.RData")
  save(beta_old, file = "beta_old.RData")
  save(a_expect, file = "a_expect.RData")
  save(b_expect, file = "b_expect.RData")
  save(k_expect, file = "k_expect.RData")
  save(a0_expect, file = "a0_expect.RData")
  save(b0_expect, file = "b0_expect.RData")
  save(pvalue, file = "pvalue.RData")
  save(trs, file = "trs.RData")
  save(SNR, file = "SNR.RData")
  save(k_star_gro, file = "k_star_gro.RData")
  write.csv(confs, file = "confs.csv")
}else{
  setwd(read_dir)
  load("S_sim.RData")
  load("k_star.RData")
  load("index_alt.RData")
  load("a_expect.RData")
  load("b_expect.RData")
  load("beta_expect.RData")
  load("beta_old.RData")
  load("k_expect.RData")
  load("a0_expect.RData")
  load("b0_expect.RData")
  load("SNR.RData")
  load("trs.RData")
  load("pvalue.RData")
  load("k_star_gro.RData")
  confs = read.csv("confs.csv")[,-1]
  
  #setwd(result_dir)
}

#####################################################################
################### Bayesian Hierarchical Modeling ##################
#####################################################################

y = as.vector(t(S_sim[-c(1,nrow(S_sim)),]))
# for the model, we ignore the case where t=0 and t=1
t = rep((1:nt)/Nt, each = ns)
dt = as.matrix(dist((1:nt)/Nt))

######################## Specify hyperpriors ########################

# sigma_._a for the shape parameter and sigma_._b for the rate parameter
mu_beta_r = rep(0, ns); tau_beta = 3
mu_c_r = rep(0, ns); tau_c = 3
mu_a_r=0; tau_a = 3
mu_b_r = rep(0, ns); tau_b = 3

sigma_beta_a = 0.1; sigma_beta_b = 0.1 
sigma_c_a = 0.1; sigma_c_b = 0.1
sigma_a_a = 0.1; sigma_a_b = 0.1
sigma_b_a = 0.1; sigma_b_b = 0.1

tau_phi = 0.5; tau_phis = 0.5; tau_phit = 0.1

####################### Define keepers ########################

keep.beta = keep.b = keep.c = matrix(0, ncol = ns, nrow = iters)
keep.mu_beta = keep.mu_c = keep.mu_b = matrix(0, ncol = ns, nrow = iters)

keep.other = matrix(0, ncol = 9, nrow = iters)
colnames(keep.other) = c("mu_a","log(a)","sigma2_beta","sigma2_c", "sigma2_a", "sigma2_b","phi","phis","phit")

###################### initial values ########################
mu_beta = rep(mean(log(-beta_expect)), ns); mu_c = mu_c_r; mu_a = mean(log(a_expect)); mu_b = log(b_expect)
prebeta = log(-beta_expect); prec = qnorm(k_expect); prea = mean(log(a_expect)); preb = log(b_expect) 
sigma2_beta = 1; sigma2_c = 1; sigma2_a = 0.5; sigma2_b = 1
phi = phis_real; phis = 2; phit = 0.2

# log likelihood for the hyperpriors
l_phi = log(dexp(phi, rate=tau_phi))
l_phis = log(dexp(phis, rate=tau_phis))
l_phit = log(dexp(phit, rate=tau_phit))

r_mu_beta = mu_beta-mu_beta_r; l_mu_beta = -ns*log(tau_beta)-1/(2*tau_beta^2)*crossprod(matrix(r_mu_beta, ncol = 1),matrix(r_mu_beta, ncol = 1))
r_mu_c = mu_c-mu_c_r; l_mu_c = -ns*log(tau_c)-1/(2*tau_c^2)*crossprod(matrix(r_mu_c, ncol = 1), matrix(r_mu_c, ncol = 1))
r_mu_a = mu_a - mu_a_r; l_mu_a = -log(tau_a)-1/(2*tau_a^2)*r_mu_a^2
r_mu_b = mu_b-mu_b_r; l_mu_b = -ns*log(tau_b)-1/(2*tau_b^2)*crossprod(matrix(r_mu_b, ncol = 1), matrix(r_mu_b, ncol = 1))

# log likelihood for process
Sigma = exp(-ds/phi)
Sigma_inv = solve(Sigma)
Sigma_det = det(Sigma_inv)
r_beta = prebeta-mu_beta; r_a = prea-mu_a; r_b = preb-mu_b; r_c = prec-mu_c; 

r_beta_sum = crossprod(matrix(r_beta, ncol=1),Sigma_inv%*%matrix(r_beta,ncol=1));
r_a_sum = r_a^2
r_b_sum = crossprod(matrix(r_b, ncol=1),Sigma_inv%*%matrix(r_b,ncol=1))
r_c_sum = crossprod(matrix(r_c, ncol=1),Sigma_inv%*%matrix(r_c,ncol=1))

l_beta = -ns/2*log(sigma2_beta)+1/2*log(Sigma_det)-1/(2*sigma2_beta)*r_beta_sum
l_a = -1/2*log(sigma2_a)-1/(2*sigma2_a)*r_a_sum
l_b = -ns/2*log(sigma2_b)+1/2*log(Sigma_det)-1/(2*sigma2_b)*r_b_sum
l_c = -ns/2*log(sigma2_c)+1/2*log(Sigma_det)-1/(2*sigma2_c)*r_c_sum

# calculate corrrelation matrix
rs <- abs(ds)/phis
rt <- abs(dt)/phit
CORt <- exp(-rt)
diag(CORt) <- 1
CORs <- exp(-rs)
diag(CORs) <- 1

beta = -exp(prebeta); a = rep(exp(prea), ns); b = exp(preb); c = pnorm(prec)
mu = rep(beta,nt)*((rep(c,nt)-1)*t+(t-rep(c,nt))*(ifelse(t>rep(c,nt),1,0)))
R <- y-mu

tau2 <- rep(a,nt)*t^2*(1-t)^2+rep(b,nt)*Nt*(rep(c,nt))^2*t*(1-t)^3*(ifelse(t>rep(c,nt),1,0))+
  rep(b,nt)*Nt*(1-rep(c,nt))^2*t^3*(1-t)*(ifelse(t<=rep(c,nt),1,0))
tauV.inv <- matrix(1/sqrt(tau2),nrow=1)
tauM <- crossprod(tauV.inv, tauV.inv)
CORt.inv <- solve(CORt)
CORs.inv <- solve(CORs)
kro.sol <- kronecker(CORt.inv, CORs.inv)
Siginv  <- tauM * kro.sol
CORt.logdet <- log(det(CORt))
CORs.logdet <- log(det(CORs))
logSigdet <- sum(log(tau2))+ns*CORt.logdet+nt*CORs.logdet
l_y <- -1/2*logSigdet-1/2*crossprod(R,Siginv%*%R)

curll_mu_beta <- l_beta+l_mu_beta
curll_mu_c <- l_c+l_mu_c

curll_beta <- l_y+l_beta
curll_a <- l_y+l_a
curll_b <- l_y+l_b
curll_c <- l_y+l_c

curll_phi <- l_beta+l_b+l_c+l_phi
curll_phis <- l_y+l_phis
curll_phit <- l_y+l_phit

############################## Run MCMC ##############################

t01 = Sys.time()
library(mvtnorm)

acc=matrix(0,ncol = 11, nrow = 2)
colnames(acc) = c("mu_beta","mu_c","mu_a","mu_b","log(-beta)","log(a)","log(b)","Phi^{-1}(c)","phi","phis","phit")
for(i in 1:iters){
  cat("interation:",i,"\n")
  t0 = Sys.time()
  ############################ mu_beta:
  acc[1,1] <- acc[1,1]+1
  curll_mu_beta <- l_beta+l_mu_beta
  
  canmu_beta <- rmvnorm(1,mu_beta,0.01*diag(ns))
  canr_mu_beta = canmu_beta-mu_beta_r 
  canl_mu_beta = -ns*log(tau_beta)-1/(2*tau_beta^2)*crossprod(matrix(canr_mu_beta,ncol=1),matrix(canr_mu_beta,ncol=1))
  canr_beta = prebeta-canmu_beta
  canr_beta_sum = crossprod(matrix(canr_beta,ncol=1),Sigma_inv%*%matrix(canr_beta,ncol=1))
  canl_beta = -ns/2*log(sigma2_beta)+1/2*log(Sigma_det)-1/(2*sigma2_beta)*canr_beta_sum
  canll_mu_beta <- canl_beta+canl_mu_beta
  MH <- canll_mu_beta-curll_mu_beta
  
  if(runif(1)<exp(MH)){
    acc[2,1] <-acc[2,1]+1
    mu_beta  <- canmu_beta; r_mu_beta <-canr_mu_beta; l_mu_beta<-canl_mu_beta
    r_beta  <- canr_beta; r_beta_sum<- canr_beta_sum;l_beta<-canl_beta
  }
  
  ############################ mu_c:
  acc[1,2] <- acc[1,2]+1
  curll_mu_c <- l_c+l_mu_c
  
  canmu_c <- rmvnorm(1,mu_c,0.01*diag(ns))
  canr_mu_c = canmu_c-mu_c_r
  canl_mu_c = -ns*log(tau_c)-1/(2*tau_c^2)*crossprod(matrix(canr_mu_c,ncol=1), matrix(canr_mu_c,ncol=1))
  canr_c = prec-canmu_c; 
  canr_c_sum = crossprod(matrix(canr_c,ncol=1),Sigma_inv%*%matrix(canr_c,ncol=1));
  canl_c = -ns/2*log(sigma2_c)+1/2*log(Sigma_det)-1/(2*sigma2_c)*canr_c_sum
  canll_mu_c <- canl_c+canl_mu_c
  MH <- canll_mu_c-curll_mu_c
  
  if(runif(1)<exp(MH)){
    acc[2,2] <-acc[2,2]+1
    mu_c  <- canmu_c; r_mu_c <-canr_mu_c; l_mu_c<-canl_mu_c
    r_c  <- canr_c; r_c_sum<-canr_c_sum; l_c<-canl_c
  }
  ############################# mu_a:
  acc[1,3] <- acc[1,3]+1
  curll_mu_a <- l_a+l_mu_a
  
  canmu_a <- rnorm(1,mu_a,4)
  canr_mu_a = canmu_a-mu_a_r
  canl_mu_a = -log(tau_a)-1/(2*tau_a^2)*canr_mu_a^2
  canr_a = prea-canmu_a; 
  canr_a_sum = canr_a^2
  canl_a = -1/2*log(sigma2_a)-1/(2*sigma2_a)*canr_a_sum
  canll_mu_a <- canl_a+canl_mu_a
  MH <- canll_mu_a-curll_mu_a
  
  if(runif(1)<exp(MH)){
    acc[2,3] <-acc[2,3]+1
    mu_a  <- canmu_a; r_mu_a <-canr_mu_a; l_mu_a<-canl_mu_a
    r_a  <- canr_a; r_a_sum<-canr_a_sum; l_a<-canl_a
  }
  ############################ mu_b:
  acc[1,4] <- acc[1,4]+1
  curll_mu_b <- l_b+l_mu_b
  
  canmu_b <- rmvnorm(1,mu_b,0.000005*diag(ns))
  canr_mu_b = canmu_b-mu_b_r
  canl_mu_b = -ns*log(tau_b)-1/(2*tau_b^2)*crossprod(matrix(canr_mu_b,ncol=1), matrix(canr_mu_b,ncol=1))
  canr_b = preb-canmu_b; 
  canr_b_sum = crossprod(matrix(canr_b,ncol=1),Sigma_inv%*%matrix(canr_b,ncol=1));
  canl_b = -ns/2*log(sigma2_b)+1/2*log(Sigma_det)-1/(2*sigma2_b)*canr_b_sum
  canll_mu_b <- canl_b+canl_mu_b
  MH <- canll_mu_b-curll_mu_b
  
  if(runif(1)<exp(MH)){
    acc[2,4] <-acc[2,4]+1
    mu_b  <- canmu_b; r_mu_b <-canr_mu_b; l_mu_b<-canl_mu_b
    r_b  <- canr_b; r_b_sum<-canr_b_sum; l_b<-canl_b
  }
  ############################ log(-beta): prebeta
  acc[1,5] <- acc[1,5]+1
  curll_beta <- l_y+l_beta
  
  sigma2_move_beta = 0.001
  canprebeta <- rmvnorm(1,prebeta,sigma2_move_beta*diag(ns))
  
  canr_beta = canprebeta-mu_beta
  canr_beta_sum = crossprod(matrix(canr_beta,ncol=1),Sigma_inv%*%matrix(canr_beta,ncol=1))
  canl_beta = -ns/2*log(sigma2_beta)+1/2*log(Sigma_det)-1/(2*sigma2_beta)*canr_beta_sum
  canbeta = -exp(canprebeta)
  canmu = rep(canbeta,nt)*((rep(c,nt)-1)*t+(t-rep(c,nt))*(ifelse(t>rep(c,nt),1,0)))
  canR <- y-canmu
  
  canl_y <- -1/2*logSigdet-1/2*crossprod(canR,Siginv%*%canR)
  canll_beta <- canl_y+canl_beta
  
  MH <- canll_beta-curll_beta
  
  if(runif(1)<exp(MH)){
    acc[2,5] <-acc[2,5]+1
    prebeta <- canprebeta; r_beta <- canr_beta; r_beta_sum<-canr_beta_sum; l_beta <- canl_beta; beta <- canbeta; 
    mu <- canmu; R <- canR; l_y <- canl_y
  }
  
  ############################ log(a): prea
  acc[1,6] <- acc[1,6]+1
  curll_a <- l_y+l_a
  
  sigma2_move_a = 0.01
  canprea <- rnorm(1, mean = prea, sd = sqrt(sigma2_move_a))
  
  canr_a = canprea-mu_a
  canr_a_sum = canr_a^2
  canl_a = -1/2*log(sigma2_a)-1/(2*sigma2_a)*canr_a_sum
  
  cana = rep(exp(canprea), ns)
  cantau2 <- (rep(cana,nt)*t^2*(1-t)^2 + rep(b,nt)*Nt*rep(c,nt)^2*t*(1-t)^3) * ifelse(t>rep(c,nt),1,0) +
    (rep(cana,nt)*t^2*(1-t)^2 + rep(b,nt)*Nt*(1-rep(c,nt))^2*t^3*(1-t)) * ifelse(t<=rep(c,nt),1,0)
  
  cantauV.inv <- matrix(1/sqrt(cantau2),nrow=1)
  cantauM <- crossprod(cantauV.inv, cantauV.inv)
  canSiginv  <- cantauM * kro.sol
  canlogSigdet <- sum(log(cantau2))+ns*CORt.logdet+nt*CORs.logdet
  canl_y <- -1/2*canlogSigdet-1/2*crossprod(R,canSiginv%*%R)
  canll_a <- canl_y+canl_a
  MH <- canll_a-curll_a
  
  if(runif(1)<exp(MH)){
    acc[2,6] <-acc[2,6]+1
    prea <- canprea; r_a <- canr_a; r_a_sum<-canr_a_sum; l_a <- canl_a;a<-cana;tau2 <- cantau2; tauV.inv <- cantauV.inv;
    tauM<-cantauM; 
    Siginv <- canSiginv; logSigdet <- canlogSigdet; l_y <- canl_y
  }
  ############################ log(b): preb
  acc[1,7] <- acc[1,7]+1
  curll_b <- l_y+l_b
  
  sigma2_move_b = 0.000004
  canpreb <- rmvnorm(1,preb,sigma2_move_b*diag(ns))
  canr_b = canpreb-mu_b
  canr_b_sum = crossprod(matrix(canr_b,ncol=1),Sigma_inv%*%matrix(canr_b,ncol=1))
  canl_b = -ns/2*log(sigma2_b)+1/2*log(Sigma_det)-1/(2*sigma2_b)*canr_b_sum
  canb = exp(canpreb)
  cantau2 <- (rep(a,nt)*t^2*(1-t)^2+rep(canb,nt)*Nt*rep(c,nt)^2*t*(1-t)^3)*ifelse(t>rep(c,nt),1,0)+
    (rep(a,nt)*t^2*(1-t)^2+rep(canb,nt)*Nt*(1-rep(c,nt))^2*t^3*(1-t))*ifelse(t<=rep(c,nt),1,0)
  cantauV.inv <- matrix(1/sqrt(cantau2),nrow=1)
  cantauM <- crossprod(cantauV.inv, cantauV.inv)
  canSiginv  <- cantauM * kro.sol
  canlogSigdet <- sum(log(cantau2))+ns*CORt.logdet+nt*CORs.logdet
  canl_y <- -1/2*canlogSigdet-1/2*crossprod(R,canSiginv%*%R)
  
  canll_b <- canl_y+canl_b
  MH <- canll_b-curll_b
  
  if(runif(1)<exp(MH)){
    acc[2,7] <-acc[2,7]+1
    preb <- canpreb; r_b <- canr_b; r_b_sum<-canr_b_sum; l_b <- canl_b; b<- canb; tau2 <- cantau2; tauV.inv <- cantauV.inv;
    tauM <- cantauM;
    Siginv <- canSiginv; logSigdet <- canlogSigdet; l_y <- canl_y
  }
  
  ############################ Phi^{-1}(c): prec
  acc[1,8] <- acc[1,8]+1
  curll_c <- l_y+l_c
  
  sigma2_move_c = 0.00001
  canprec <- rmvnorm(1,prec,sigma2_move_c*diag(ns))
  
  canr_c = canprec-mu_c
  canr_c_sum = crossprod(matrix(canr_c,ncol=1),Sigma_inv%*%matrix(canr_c,ncol=1))
  canl_c = -ns/2*log(sigma2_c)+1/2*log(Sigma_det)-1/(2*sigma2_c)*canr_c_sum
  canc <- pnorm(canprec)
  canmu = rep(beta,nt)*((rep(canc,nt)-1)*t+(t-rep(canc,nt))*(ifelse(t>rep(canc,nt),1,0)))
  canR <- y-canmu
  cantau2 <- (rep(a,nt)*t^2*(1-t)^2+rep(b,nt)*Nt*rep(canc,nt)^2*t*(1-t)^3)*ifelse(t>rep(canc,nt),1,0)+
    (rep(a,nt)*t^2*(1-t)^2+rep(b,nt)*Nt*(1-rep(canc,nt))^2*t^3*(1-t))*ifelse(t<=rep(canc,nt),1,0)
  cantauV.inv <- matrix(1/sqrt(cantau2),nrow=1)
  cantauM <- crossprod(cantauV.inv, cantauV.inv)
  canSiginv  <- cantauM * kro.sol
  canlogSigdet <- sum(log(cantau2))+ns*CORt.logdet+nt*CORs.logdet
  canl_y <- -1/2*canlogSigdet-1/2*crossprod(canR,canSiginv%*%canR)
  
  canll_c <- canl_y+canl_c
  MH <- canll_c-curll_c
  
  if(runif(1)<exp(MH)){
    acc[2,8] <-acc[2,8]+1
    prec <- canprec; r_c <- canr_c; r_c_sum<-canr_c_sum; l_c <- canl_c; c <- canc;
    mu <- canmu; R <- canR; tau2 <- cantau2; tauV.inv <- cantauV.inv;
    tauM<- cantauM;
    Siginv <- canSiginv; logSigdet <- canlogSigdet; l_y <- canl_y
  }
  
  ############################ phi
  acc[1,9] <- acc[1,9]+1
  curll_phi <- l_beta+l_b+l_c+l_phi
  
  sigma2_move_phi = 0.0004
  canphi <- rnorm(1, mean = phi, sd = sqrt(sigma2_move_phi))
  while(canphi<=0) {canphi <- rnorm(1, mean = phi, sd = sqrt(sigma2_move_phi))}
  canl_phi = log(dexp(canphi, rate=tau_phi))
  canSigma = exp(-ds/canphi)
  canSigma_inv = solve(canSigma)
  canSigma_det = det(canSigma_inv)
  canr_beta_sum = crossprod(matrix(r_beta,ncol=1),canSigma_inv%*%matrix(r_beta,ncol=1))
  canr_b_sum = crossprod(matrix(r_b,ncol=1),canSigma_inv%*%matrix(r_b,ncol=1))
  canr_c_sum = crossprod(matrix(r_c,ncol=1),canSigma_inv%*%matrix(r_c,ncol=1))
  
  canl_beta = -ns/2*log(sigma2_beta)+1/2*log(canSigma_det)-1/(2*sigma2_beta)*canr_beta_sum
  canl_b = -ns/2*log(sigma2_b)+1/2*log(canSigma_det)-1/(2*sigma2_b)*canr_b_sum
  canl_c = -ns/2*log(sigma2_c)+1/2*log(canSigma_det)-1/(2*sigma2_c)*canr_c_sum
  canll_phi <- canl_beta+canl_b+canl_c+canl_phi
  MH <- canll_phi-curll_phi-log(1-pnorm(-canphi/sqrt(sigma2_move_phi)))+log(1-pnorm(-phi/sqrt(sigma2_move_phi)))
  
  if(runif(1)<exp(MH)){
    acc[2,9] <-acc[2,9]+1
    phi <- canphi; Sigma <- canSigma; l_phi <- canl_phi; Sigma_inv <- canSigma_inv; Sigma_det<- canSigma_det;
    l_beta <- canl_beta; l_b <- canl_b; l_c <- canl_c
    r_beta_sum <- canr_beta_sum ; r_b_sum<-canr_b_sum;r_c_sum<-canr_c_sum
  }
  
  ############################ phis
  acc[1,10] <- acc[1,10]+1
  curll_phis <- l_y+l_phis
  
  sigma2_move_phis = 0.1
  canphis <- rnorm(1, mean = phis, sd = sqrt(sigma2_move_phis))
  while(canphis<=0){canphis <- rnorm(1, mean = phis, sd = sqrt(sigma2_move_phis))}
  
  canl_phis = log(dexp(canphis, rate=tau_phis))
  canrs <- abs(ds)/canphis
  canCORs <- exp(-canrs)
  diag(canCORs) <- 1
  canCORs.inv <- solve(canCORs)
  cankro.sol <- kronecker(CORt.inv, canCORs.inv)
  canSiginv  <- tauM * cankro.sol
  canCORs.logdet <- log(det(canCORs))
  canlogSigdet <- sum(log(tau2))+ns*CORt.logdet+nt*canCORs.logdet
  canl_y <- -1/2*canlogSigdet-1/2*crossprod(R,canSiginv%*%R)
  
  canll_phis <- canl_y+canl_phis
  
  MH <- canll_phis-curll_phis-log(1-pnorm(-canphis/sqrt(sigma2_move_phis)))+log(1-pnorm(-phis/sqrt(sigma2_move_phis)))
  if(runif(1)<exp(MH)){
    acc[2,10] <-acc[2,10]+1
    phis <- canphis; l_phis <- canl_phis; rs <- canrs; CORs <- canCORs; 
    CORs.inv <-canCORs.inv; kro.sol <- cankro.sol; 
    Siginv <- canSiginv; CORs.logdet <- canCORs.logdet; logSigdet <- canlogSigdet; l_y<- canl_y
  }
  
  ############################ phit
  acc[1,11] <- acc[1,11]+1
  curll_phit <- l_y+l_phit
  
  sigma2_move_phit = 0.01
  canphit <- rnorm(1, mean = phit, sd = sqrt(sigma2_move_phit))
  while(canphit<=0){canphit <- rnorm(1, mean = phit, sd = sqrt(sigma2_move_phit))}
  
  canl_phit = log(dexp(canphit, rate=tau_phit))
  canrt <- abs(dt)/canphit
  canCORt <- exp(-canrt)
  diag(canCORt) <- 1
  canCORt.inv <- solve(canCORt)
  cankro.sol <- kronecker(canCORt.inv, CORs.inv)
  canSiginv  <- tauM * cankro.sol
  canCORt.logdet <- log(det(canCORt))
  canlogSigdet <- sum(log(tau2))+ns*canCORt.logdet+nt*CORs.logdet
  canl_y <- -1/2*canlogSigdet-1/2*crossprod(R,canSiginv%*%R)
  
  canll_phit <- canl_y+canl_phit
  
  MH <- canll_phit-curll_phit-log(1-pnorm(-canphit/sqrt(sigma2_move_phit)))+log(1-pnorm(-phit/sqrt(sigma2_move_phit)))
  if(runif(1)<exp(MH)){
    acc[2,11] <-acc[2,11]+1
    phit <- canphit; l_phit <- canl_phit; rt <- canrt; CORt <- canCORt; 
    CORt.inv <- canCORt.inv; kro.sol <- cankro.sol; 
    Siginv <- canSiginv; CORt.logdet <- canCORt.logdet; logSigdet <- canlogSigdet; l_y<- canl_y
  }
  ############# sigmas (Gibbs) #################
  
  sigma2_beta <- rinvgamma(1, shape = sigma_beta_a+ns/2, rate=sigma_beta_b+r_beta_sum/2)
  sigma2_c <- rinvgamma(1, shape = sigma_c_a+ns/2, rate=sigma_c_b+r_c_sum/2)
  sigma2_a <- rinvgamma(1, shape = sigma_a_a+1/2, rate=sigma_a_b+r_a_sum/2)
  sigma2_b <- rinvgamma(1, shape = sigma_b_a+ns/2, rate=sigma_b_b+r_b_sum/2)
  
  ##############################################:
  #####        KEEP TRACK OF STUFF       #######:
  ##############################################:
  keep.beta[i,] = prebeta;
  keep.b[i,] = preb; keep.c[i,] = prec
  keep.mu_beta[i,] = mu_beta; keep.mu_c[i,] = mu_c; keep.mu_b[i,] = mu_b
  keep.other[i,] = c(mu_a, prea, sigma2_beta, sigma2_c,sigma2_a, sigma2_b, phi,phis,phit)
  
  cat(acc[1,], "/", acc[2,],"\n")
  if(walltime-as.numeric(Sys.time()-time_begin, units = "hours")<=15/3600){
    #setwd(result_dir)
    write.csv(keep.beta, "prebeta_temp.csv")
    write.csv(keep.b, "preb_temp.csv")
    write.csv(keep.c, "prec_temp.csv")
    write.csv(keep.other, "other_temp.csv")
    write.csv(keep.mu_beta, "mu_beta_temp.csv")
    write.csv(keep.mu_c, "mu_c_temp.csv")
    write.csv(keep.mu_b, "mu_b_temp.csv")
    write.csv(acc,"acceptance.csv")
  }
  cat(Sys.time()-t0, "\n")
}

##############################################:
########        SAVE RESULTS       ###########:
##############################################:
#setwd(result_dir)
write.csv(keep.beta, "prebeta_temp.csv")
write.csv(keep.b, "preb_temp.csv")
write.csv(keep.c, "prec_temp.csv")
write.csv(keep.other, "other_temp.csv")
write.csv(keep.mu_beta, "mu_beta_temp.csv")
write.csv(keep.mu_c, "mu_c_temp.csv")
write.csv(keep.mu_b, "mu_b_temp.csv")
write.csv(acc,"acceptance.csv")
cat(acc[2,]/acc[1,],"\n")
acc[2,]/acc[1,]
Sys.time()-t01
r_thin = seq(iters*3/4+1, iters, thin_step) 
c.mean_alt = apply(pnorm(as.matrix(keep.c[r_thin,index_alt])), 2, mean)
c.q1_alt = apply(pnorm(as.matrix(keep.c[r_thin,index_alt])), 2, function(x) quantile(x, probs=0.05/2))
c.q2_alt = apply(pnorm(as.matrix(keep.c[r_thin,index_alt])), 2, function(x) quantile(x, probs=1-0.05/2))
logST = log(c.q2_alt - c.q1_alt)
logFF = log(confs[index_alt,2] - confs[index_alt,1])

results = matrix(0, ncol = 17, nrow = 1)
results[1] = sum(k_star[index_alt]>=c.q1_alt & k_star[index_alt]<=c.q2_alt)/ns_alt
results[2] = sum(k_star[index_alt]>=confs[index_alt,1] & k_star[index_alt]<=confs[index_alt,2])/ns_alt
results[3] = mean((c.mean_alt-k_star[index_alt])^2)
results[4] = mean((k_expect[index_alt]-k_star[index_alt])^2)
results[5] = mean((k_star_gro-k_star[index_alt])^2)
results[6] = sqrt(mean((c.mean_alt-k_star[index_alt])^2))
results[7] = sqrt(mean((k_expect[index_alt]-k_star[index_alt])^2))
results[8] = sqrt(mean((k_star_gro-k_star[index_alt])^2))
results[9] = mean(abs(c.mean_alt-k_star[index_alt]))
results[10] = mean(abs(k_expect[index_alt]-k_star[index_alt]))
results[11] = mean(abs(k_star_gro-k_star[index_alt]))
results[12] = k_star_gro
# For the following ones, we round our estimate to the closest integer and save the corresponding
# MSE, RMSE, MAE
c.mean_alt_round = round(c.mean_alt * Nt)/Nt
results[13] = mean((c.mean_alt_round-k_star[index_alt])^2)
results[14] = sqrt(mean((c.mean_alt_round-k_star[index_alt])^2))
results[15] = mean(abs(c.mean_alt_round -k_star[index_alt]))
results[16] = mean(logST)
results[17] = mean(logFF)
                        
colnames(results) = c("Coverage ST", "Coverage FF", "MSE ST", "MSE FF","MSE G", "RMSE ST", "RMSE FF",
                      "RMSE G", "MAE ST", "MAE FF", "MAE G", "Estimate G", "MSE rST", "RMSE rST", 
                      "MAE rST", "mean(logST)", "mean(logFF)")
write.csv(results,"results.csv")
save(logST, file = "logST.RData")
save(logFF, file = "logFF.RData")

setwd(summary_dir)
write.csv(results,paste("results_", phis_real, "_", signal, "_", job, ".csv", sep = ""))
save(logST, file = paste("logST", phis_real, "_", signal, "_", job, ".RData", sep = ""))
save(logFF, file = paste("logFF", phis_real, "_", signal, "_", job, ".RData", sep = ""))
#setwd(result_dir)

##############################################:
########        SAVE CI IMAGES       #########:
##############################################:
result_new = matrix(0, ncol = 5, nrow = 0)
result_new = as.data.frame(result_new)
colnames(result_new) = c("location", "estimate", "lower", "upper", "Method")

for(i in 1:length(index_alt)){
  result_add = data.frame(matrix(0, ncol = 5, nrow = 3))
  colnames(result_add) = c("location", "estimate", "lower", "upper", "Method")
  
  result_add[, 1] = index_alt[i]
  result_add[, 2] = c(k_star[index_alt[i]], k_expect[index_alt[i]], c.mean_alt[i])
  result_add[, 3] = c(k_star[index_alt[i]], confs[index_alt[i],1], c.q1_alt[i])
  result_add[, 4] = c(k_star[index_alt[i]], confs[index_alt[i],2], c.q2_alt[i])
  result_add[, 5] = factor(c("True", "FFChange","ST_Model"))
  result_new = rbind(result_new, result_add)
}

png(paste("0.05CI_", phis_real, "_", signal, "_", job, ".png", sep = ""))
ggplot() +
  geom_point(data = result_new, aes(location, estimate, col = Method, shape = Method, size = Method)) +
  ggtitle(paste("Confidence Interval for alternative locations with Level 0.05 with coverage", 
                round(sum(k_star[index_alt]>=c.q1_alt & k_star[index_alt]<=c.q2_alt)/ns_alt,3)
                , "and Aue:", round(sum(k_star[index_alt]>=confs[index_alt,1] & k_star[index_alt]<=confs[index_alt,2])/ns_alt,3))
  ) +
  xlab("Location")+
  ylab("Changepoint")+
  geom_errorbar(data = result_new[result_new$Method!="True",], aes(location, estimate, ymin=lower, ymax=upper, col = Method), size = 1.01)+
  scale_color_manual(values = c('blue4', 'orangered2', 'black'))+
  scale_shape_manual(values = c(4,4,20))+
  scale_size_manual(values = c(3,3,3.1))
dev.off()
