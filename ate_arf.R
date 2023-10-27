# Call library
library("lme4")       # For multilevel models
library("estimatr")

# This is an R function that calculates the ATE-ARF estimator developed in
# Schochet (Statistics in Medicine, 2024). The associated documentation, ate-arf.pdf, discusses the
# function inputs and how to run the program.

# The ate_arf function

ate_arf <- function(data_csv,y_var,trt_var,x_logit,x_wls,clus_var,wgt,marg,se_ipw,df_wls,out_regr,out_est) {

# Read in csv data

mydat1 <- read.csv(data_csv)

sum_dat1  <- summary(mydat1)
nobs_dat1 <- NROW(mydat1)

# Define defaults

if ((wgt == '0') | (wgt == "")) {
    wgt <- c(0)
    }

if ((marg == '1') | (marg == "")) {
    marg <- c(1)
    } else if (marg == '0') {
      marg <- c(0)
      }

if ((se_ipw == '1') | (se_ipw == "")) {
  se_ipw <- c(1)
  } else if (se_ipw == '0') {
    se_ipw <- c(0)
    }

if (out_regr == "") {
  out_regr  <- c('ate_arf log.txt')
  }

if (out_est == "") {
  out_est  <- c('ate_arf estimates.txt')
  }

# Add the weight onto mydat

if (wgt == 0) {
  mydat1['ww'] <- 1
  } else {
      mydat1['ww'] <- mydat1[wgt]
      }

# Add on treatment variable and create control variable

mydat1['trt_var']   <- mydat1[trt_var]
mydat1['cnt_var']   <- 1 - mydat1$trt_var

# Add on cluster variable

mydat1['clusv'] <- mydat1[clus_var]

mobs_dat1 <- NROW(unique(mydat1$clusv))

# Keep only cases with treatment = 0,1 and >0 wgts

mydat1 <- mydat1[which(((mydat1$trt_var==0)  | (mydat1$trt_var==1))),]
mydat1 <- mydat1[which(mydat1$ww > 0),]

# Check whether x_wls has elements

hasx <- 0
if (any(is.element(x_wls,'0')) == FALSE | any(is.element(x_wls,0)) == FALSE) {
  hasx <- 1
  }

#
# Conduct a complete case analysis so add on other X and Y variables
#

compl  <- c('trt_var','cnt_var','ww','clusv')

# Add on x_logit variables

compla <- c("")
for (i in 1:NROW(x_logit)) {
  vname <- paste0("zzxl", i)
  mydat1[[vname]] <- as.numeric(unlist(mydat1[x_logit[[i]]]))

  if (i == 1) {
    compla <- c('zzxl1')
    } else {
      compla <- c(compla,paste0("zzxl", i))
      }
    }

compl <- c(compl,compla)

# Add x_wls if specified

if (hasx == 1) {

  compla <- c("")
  for (i in 1:NROW(x_wls)) {
    vname <- paste0("zzx", i)
    mydat1[[vname]] <- as.numeric(unlist(mydat1[x_wls[[i]]]))

    if (i == 1) {
      compla <- c('zzx1')
      } else {
        compla <- c(compla,paste0("zzx", i))
        }
      }

  compl <- c(compl,compla)
  }

# DO OVER Y

for (iy in 1:NROW(y_var)) {

yname <- y_var[[iy]]
vname <- paste0("zzy", iy)
mydat1[[vname]] <- as.numeric(unlist(mydat1[y_var[[iy]]]))

compla <- c(paste0("zzy", iy))

complg <- c(compl,compla)

###
# Take complete cases only
###

mydat <- mydat1[complete.cases(mydat1[ , complg]), ]

nobs_dat <- NROW(mydat)
mobs_dat <- NROW(unique(mydat$clusv))

# Normalize the ww weights to nt and nc

mydat['ttt'] <- mydat[trt_var]
mydat['ccc'] <- 1 - mydat[trt_var]

normt <- sum(mydat$ww[mydat$ttt == 1])
n_wwt <- NROW(mydat$ww[mydat$ttt == 1])

mydat$ww[mydat$ttt == 1] <- mydat$ww[mydat$ttt == 1]*n_wwt/normt

normc <- sum(mydat$ww[mydat$ttt == 0])
n_wwc <- NROW(mydat$ww[mydat$ttt == 0])

mydat$ww[mydat$ttt == 0] <- mydat$ww[mydat$ttt == 0]*n_wwc/normc

# Create X covariates for logit regressions

xldat <- mydat[x_logit]

# Check whether x_wls has elements and create X covariates for WLS regression

if (hasx == 1) {
  xdat <- mydat[x_wls]
  }

###
# RUN LOGIT MODEL FOR THE CONTROL GROUP RECRUITS VERSUS TREATMENT GROUP RECRUITS
###

# Create ppc dependent variable for control logit model

mydat$ppc <- mydat$cnt_var

# Create logit data for control model

cldat      <- cbind(mydat['ppc'],xldat)
cldat      <- cbind(cldat,mydat['ww'])
cldat_clus <- cbind(cldat,mydat['clusv'])

# Create formulas for control logit model

formc <- as.formula(paste("ppc ~ ", paste(colnames(xldat), collapse = " + ")))
formc_clus <- as.formula(paste("ppc ~ ", paste(colnames(xldat), collapse = " + "),"+ (1|clusv)"))

# Fit the logistic regression model for the control group

if (marg == 1) {
  logc <- suppressWarnings(glm(formc,
              family="binomial",
              weights=ww,
              data=cldat))
  predc <- predict(logc, newdata = cldat, type = 'response')

  } else if ((marg == 0) | (marg == '0')) {
    logc <- suppressWarnings(glmer(formc_clus,
                  family="binomial",
                  weights=ww,
                  nAGQ=10,
                  data=cldat_clus))
    predc <- predict(logc, newdata = cldat, re.form=NA, type = 'response')
    }

#view model summary

logc_res <- summary(logc)

# Add e(x) (predicted probs) back to the main dataset, calculate (1-e(x)), and calculate residuals

mydat$ex    <- predc
mydat$ex1   <- (1 - predc)

mydat$resid_lc <- mydat$cnt_var - mydat$ex

# Construct ipw weights

ccc <- mydat$cnt_var
ex    <- predc
ex1   <- (1 - predc)

mydat['wgt_arf'] <- ccc  + (1-ccc)*(ex/ex1)

# Multiply weight by sample weight by ww

mydat$wgt_arf <- mydat$wgt_arf*mydat$ww

# Run WLS regressions

# Create WLS formula

if (hasx == 1) {
  form_wls <- as.formula(paste("yyy ~ ttt +", paste(colnames(xdat), collapse = " + ")))
  } else if (hasx == 0) {
    form_wls <- as.formula(paste("yyy ~ ttt"))
    }

# Create WLS data

mydat['yyy'] <- mydat[yname]
wlsdat <- cbind(mydat['yyy'],mydat['ttt'])
wlsdat <- cbind(wlsdat,mydat['ccc'])
wlsdat <- cbind(wlsdat,mydat['ww'])
wlsdat <- cbind(wlsdat,mydat['clusv'])

# Calculate mean and std of yyy for the control group

yyc <- wlsdat$yyy[wlsdat$ttt==0]

mean_yc <- mean(yyc)
std_yc  <- sd(yyc)

if (hasx == 1) {
  wlsdat <- cbind(wlsdat,xdat)
  }

wlsdat_arf <- cbind(wlsdat,mydat['wgt_arf'])

###
# RUN WLS REGRESSIONS
###

###
# DO THE REGRESSIONS IF DO NOT WANT ADJ IPW STANDARD ERRORS - USE LM_ROBUST EXISTING PACKAGES
###

if (se_ipw == 0) {

sand_arf   <- lm_robust(form_wls, data = wlsdat_arf, weights = wgt_arf, clusters = clusv,se_type = "stata")

kk <- 0

if (hasx == 1) {
  kk         <- length(sand_arf$coefficients)        # Number of coefficients
}

ate_arf_res     <- summary(sand_arf)
imp_arf         <- sand_arf$coefficients[2]
sand_se_arf     <- sand_arf$std.error[2]
ts_arf          <- sand_arf$statistic[2]
p_arf           <- sand_arf$p.value[2]
nobs_arf        <- sand_arf$nobs
mclus_arf       <- sand_arf$nclusters

} # end if (se_ipw == 0)

###
# DO THE REGRESSIONS IF WANT ADJ IPW STANDARD ERRORS
###

if (se_ipw == 1) {

ate_arf <- lm(form_wls, data = wlsdat_arf, weights = wgt_arf)

ate_arf_res <- summary(ate_arf)
imp_arf     <- ate_arf$coefficients[2]
resid_arf   <- ate_arf$residuals
nobs_arf    <- NROW(wlsdat_arf)

kk_arf      <- length(ate_arf$coefficients)

#print(ate_arf_res)

nx <- 0
if (hasx == 1) {
  ncoeff <- length(ate_arf$coefficients)     # Number of coefficients
  betax  <- ate_arf$coefficients[3:ncoeff]   # Pull off x coefficients
  nx     <- NROW(betax)                      # Number of covariates
}

###
# VARIANCE CALCULATIONS FOR ADJ ESTIMATORS
###

if (hasx == 1) {
  xdat1 <- cbind(mydat['clusv'],xdat)   # Add clus to xdat
  }

xldat1  <- cbind(mydat['clusv'],xldat)  # Add clus to xldat

# Create u vectors with WLS residuals, eta vector with logit residuals, ex, and add clusv and same for weights

u_arf    <- cbind(mydat['clusv'],resid_arf)
eta_arf  <- cbind(mydat['clusv'],mydat['resid_lc'])
wgt_arf  <- cbind(mydat['clusv'],mydat['wgt_arf'])
e_arf    <- cbind(mydat['clusv'],mydat['ex'])

nk <- NROW(x_logit)              # Number of covariates in logit model

mclus     <- length(unique(mydat$clusv))  # Number of clusters
mclus_arf <- mclus

# INITIALIZE VARIABLES FOR VARIANCES FOR ATE-ARF IPW ESTIMATOR

# ATE-ARF Variables to calculate WLS GEE

phij   <- matrix(0,(nx+2),1)
gammaj <- matrix(0,(nx+2),(nx+2))

delta <- gammaj
gamma <- gammaj

# ATE-ARF Logit

phij_l   <- matrix(0,(nk+1),1)
gammaj_l <- matrix(0,(nk+1),(nk+1))

delta_l <- gammaj_l
gamma_l <- gammaj_l

# ATE-ARF Variables to calculate Off-Diagonal WLS and Logit GEE - the B Matrix in Paper

gammaj_b <- matrix(0,(nx+2),(nk+1))
gamma_b  <- gammaj_b

# ATE-ARF Variables for all pieces

nv <- nx + nk + 3

phij_all <- matrix(0,nv,1)
phi_all  <- phij_all

gammaj_all <- matrix(0,nv,nv)
gamma_all  <- gammaj_all

delta_all <- gammaj_all

###
# DO OVER ALL CLUSTERS
###

for (j in 1:mclus) {

  # CALCULATE GEE VARIANCE PIECES FOR ATE_ARF IPW ESTIMATOR

  # ATE-ARF WLS variance

  tj <- mydat$trt_var[mydat$clusv == j]
  cj <- mydat$cnt_var[mydat$clusv == j]
  qj <- cj

  onesj <- ifelse(tj == 1,1,1)

  if (nx > 0) {
    xj  <- xdat1[xdat1$clusv == j,2:(nx+1)]
    xj  <- as.matrix(xj)
    }

  uj  <- u_arf[u_arf$clusv == j,2]

  wj  <- wgt_arf[wgt_arf$clusv == j,2]

  wuj <- wj*uj
  wtj <- wj*tj
  wcj <- wj*cj

  if (nx > 0) {
    wxj <- sqrt(wj)*xj
    }

  phij[1,1] <- t(tj) %*% wuj
  phij[2,1] <- t(cj) %*% wuj

  if (nx>0) {
    phij[3:(2+nx),1] <- t(xj) %*% wuj
    }

  deltaj <- phij %*% t(phij)
  delta  <- delta + deltaj

  gammaj[1,1]        <- t(tj) %*% wj
  gammaj[1,2]        <- 0
  gammaj[2,1]        <- 0
  gammaj[2,2]        <- t(cj) %*% wj

  if (nx>0) {
    gammaj[1,3:(2+nx)] <- t(wtj) %*% xj
    gammaj[2,3:(2+nx)] <- t(wcj) %*% xj
    gammaj[3:(2+nx),1] <- t(xj) %*% wtj
    gammaj[3:(2+nx),2] <- t(xj) %*% wcj
    gammaj[3:(2+nx),3:(2+nx)] <- t(wxj) %*% wxj
    }

  gamma <- gamma + gammaj

  # ATE_ARF B Matrix

  ej     <- e_arf[e_arf$clusv == j,2]
  etaj   <- eta_arf[eta_arf$clusv == j,2]

  xlogitj <- xldat1[xldat1$clusv == j,2:(nk+1)]
  xlogitj  <- as.matrix(xlogitj)

  nrue    <- uj*wj
  tnrue   <- tj*uj*wj
  if (nx>0) {
    tnruex  <- tnrue*xj
  }

  gammaj_b[1,1]        <- -t(tj) %*% nrue
  gammaj_b[1,2:(nk+1)] <- -t(tnrue) %*% xlogitj

  if (nx>0) {
    gammaj_b[3:(nx+2),1]        <- -t(xj) %*% tnrue
    gammaj_b[3:(nx+2),2:(nk+1)] <- -t(tnruex) %*% xlogitj
  }

  gamma_b <- gamma_b + gammaj_b

  # ATE-ARF logit variance

  eej   <- ej*(1-ej)
  eexj  <- sqrt(eej)*xlogitj

  vvv <- t(onesj) %*% etaj

  phij_l[1,1]        <- t(onesj) %*% etaj
  phij_l[2:(nk+1),1] <- t(xlogitj) %*% etaj

  deltaj_l <- phij_l %*% t(phij_l)
  delta_l  <- delta_l + deltaj_l

  gammaj_l[1,1]        <- t(onesj) %*% eej
  gammaj_l[1,2:(nk+1)] <- t(eej) %*% xlogitj
  gammaj_l[2:(nk+1),1] <- t(xlogitj) %*% eej
  gammaj_l[2:(nk+1),2:(nk+1)] <- t(eexj) %*% eexj

  gamma_l <- gamma_l + gammaj_l

  # ATE-ARF The whole shebang by putting A, B, and C matrixes together

  phij_all[1:(nx+2),1]         <- phij
  phij_all[(nx+3):(nx+nk+3),1] <- phij_l

  deltaj_all <- phij_all %*% t(phij_all)

  delta_all  <- delta_all + deltaj_all

  gammaj_all[1:(nx+2),1:(nx+2)] <- gammaj

  gammaj_all[1:(nx+2),(nx+3):nv]<- gammaj_b

  gammaj_all[(nx+3):nv,(nx+3):nv] <- gammaj_l

  gamma_all <- gamma_all + gammaj_all

  }  # For j

###
# CALCULATE FULL GEE SEs
###

# ATE-ARF WLS

delta_hat <- delta/mclus
gamma_hat <- gamma/mclus

ig <- solve(gamma_hat)     # Solve is the inverse function in R

# WLS All variables

sand_vara <- ig %*% delta_hat %*% t(ig)
sand_vara <- sand_vara/mclus
sand_sea  <- sqrt(diag(sand_vara))

# ATE-ARF WLS Impacts by differencing T and C mean estimates

lamda <- matrix(0,1,(nx+2))
lamda[1,1] <-  1
lamda[1,2] <- -1

sand_var <- lamda %*% sand_vara %*% t(lamda)
sand_se_arf  <- sqrt(sand_var)

seadj_arf    <- sqrt((mclus_arf/(mclus_arf-1))*((nobs_arf-1)/(nobs_arf-kk_arf)))
sand_se_arf  <- sand_se_arf*seadj_arf

# Logits - All Variables

delta_hat_l <- delta_l/mclus
gamma_hat_l <- gamma_l/mclus

ig_l <- solve(gamma_hat_l)

sand_vara_l <- ig_l %*% delta_hat_l %*% t(ig_l)
sand_vara_l <- sand_vara_l/mclus
sand_sea_l  <- sqrt(diag(sand_vara_l))

# ATE-ARF Whole shebang

delta_hat_all <- delta_all/mclus
gamma_hat_all <- gamma_all/mclus

ig_all = solve(gamma_hat_all)

# ATE-ARF Shebang All variables

sand_vara_all <- ig_all %*% delta_hat_all %*% t(ig_all)
sand_vara_all <- sand_vara_all/mclus
sand_sea_all  <- sqrt(diag(sand_vara_all))

# ATE-ARF Shebang Impacts by differencing T and C mean estimates

lamda_all = matrix(0,1,nv)
lamda_all[1,1] <-  1
lamda_all[1,2] <- -1

sand_var_all <- lamda_all %*% sand_vara_all %*% t(lamda_all)
sand_se_all_arf  <- sqrt(sand_var_all)

sand_se_all_arf  <- sand_se_all_arf*seadj_arf

# Conduct significance testing

if (df_wls == 1) {
  df_arf  <- mclus - nx - 2
}
else if (df_wls == 0) {
  df_arf  <- mclus - nx - nk - 3
}

ts_arf <- abs(imp_arf)/sand_se_arf
p_arf  <- 2*pt(ts_arf,  df=df_arf,  lower.tail=FALSE)

ts_all_arf <- abs(imp_arf)/sand_se_all_arf
p_all_arf  <- 2*pt(ts_all_arf,  df=df_arf,  lower.tail=FALSE)

} # END IF (se_ipw == 1)

mean_yt <- mean_yc + imp_arf

imp_arf_es <- imp_arf / std_yc

Sum_data  <- sum_dat1
Persons   <- nobs_dat1
Clusters  <- mobs_dat1

orig_samp  <- cbind(Persons,Clusters)
orig_samp1 <- as.data.frame(orig_samp)

Persons  <- nobs_dat
Clusters <- mobs_dat

anal_samp  <- cbind(Persons,Clusters)
anal_samp1 <- as.data.frame(anal_samp)

# Print Data Summary and Logit and WLS Regression Results

tits1  <- sprintf("SUMMARY OF INPUT CSV DATA FILE")
tits2  <- sprintf("ORIGINAL SAMPLE SIZES")
tits3  <- sprintf("ANALYSIS SAMPLE SIZES FOR OUTCOME %s AFTER REMOVING MISSING DATA",yname)

titl1  <- sprintf("LOGIT RESULTS PREDICTING CONTROL GROUP RECRUITMENT FOR OUTCOME %s",yname)
titl1a <- sprintf("WITH WEIGHTS %s",wgt)

titw1 <- sprintf("WLS RESULTS FOR THE ATE-ARF ESTIMATOR FOR OUTCOME %s",yname)

blank <- c("")

# PRINT DATA SUMMARY AND LOGIT AND REGRESSION RESULTS

if (iy == 1) {
  d_app <- c('FALSE')      # Do not append the output files
  } else {
    d_app <- c('TRUE')     # Append the output files with new output
    }

if (iy == 1) {
  #cat(blank,tits1,blank, sep="\n")
  #print(Sum_data)
  cat(blank,tits2,blank, sep="\n")
  print(orig_samp1, row.names = F)
}

cat(blank,tits3,blank, sep="\n")
print(anal_samp1, row.names = F)

sink(out_regr,d_app)       # Sink writes out results

if (iy == 1) {
  cat(blank,tits1,blank, sep="\n")
  print(Sum_data)
  cat(blank,tits2,blank, sep="\n")
  print(orig_samp1, row.names = F)
  }

cat(blank,tits3,blank, sep="\n")
print(anal_samp1, row.names = F)

if (wgt == 0) {
  cat(blank,titl1,blank, sep="\n")
  print(logc_res)

  cat(blank,titw1,blank, sep="\n")
  print(ate_arf_res)

} else {
  cat(blank,titl1,titl1a,blank, sep="\n")
  print(logc_res)

  cat(blank,titw1,titl1a,blank, sep="\n")
  print(ate_arf_res)

  }

sink()     # Turns off sink

# Create results for printing

if (se_ipw == 1) {
  # Merge the impact, se, tstat, and pvalue columns
    resg <- do.call("cbind", list(mean_yt, mean_yc, imp_arf, imp_arf_es, sand_se_all_arf, ts_all_arf, p_all_arf,
                                     sand_se_arf, ts_arf, p_arf))
  } else if (se_ipw == 0) {
    resg <- do.call("cbind", list(mean_yt, mean_yc, imp_arf, imp_arf_es, sand_se_arf, ts_arf, p_arf))
    }

# Format
Mean_T  <- sprintf("%.4f",resg[,1])
Mean_C  <- sprintf("%.4f",resg[,2])
Impact  <- sprintf("%.4f",resg[,3])
Effect_size  <- sprintf("%.4f",resg[,4])
Std_err <- sprintf("%.4f",resg[,5])
t_value <- sprintf("%.4f",resg[,6])
p_value <- sprintf("%.4f",resg[,7])

# Create Estimator name variable
Estimator <- c("ATE-ARF")
Estimator <- format(Estimator, justify = "left", width = 9)

# Merge again
resg1 <- cbind(Estimator,Mean_T)
resg1 <- cbind(resg1,Mean_C)
resg1 <- cbind(resg1,Impact)
resg1 <- cbind(resg1,Effect_size)
resg1 <- cbind(resg1,Std_err)
resg1 <- cbind(resg1,t_value)
resg1 <- cbind(resg1,p_value)

# Make it a data frame
resg2 <- as.data.frame(resg1)

# Add significance stars
resg2$signif <- c("")
resg2$signif[which(resg2$p_value <= .01)] <- c('***')
resg2$signif[which((resg2$p_value <= .05) & (resg2$p_value > .01))] <- c('**')
resg2$signif[which((resg2$p_value <= .10) & (resg2$p_value > .05))] <- c('*')

resg2$signif <- format(resg2$signif, justify = "left", width = 5)

# Add blank row for printing
blank <- c("")
resg2 <- rbind(blank,resg2)

# PRINT ATE-ARF ESTIMATES

if (marg == 1) {
  titc1 <- sprintf("Table %s. ATE-ARF GEE ESTIMATION RESULTS FOR THE MARGINAL LOGIT MODEL",iy)
  } else if (marg == 0) {
    titc1 <- sprintf("Table %s. ATE-ARF GEE ESTIMATION RESULTS FOR THE RANDOM EFFECTS LOGIT MODEL",iy)
    }

    titc2 <- sprintf("FOR OUTCOME VARIABLE %s: CLUSTERED DESIGN USING %s",yname,clus_var)

if (se_ipw == 1) {
  titc3 <- c("(Standard errors adjust for estimation error in the IPW weights)")
  } else if (se_ipw == 0) {
    titc3 <- c("(Standard errors ignore estimation error in the IPW weights)")
    }

titc4a <- c("Notes. See Schochet (Statistics in Medicine, 2023) for details")
titc4b <- c("on the ATE-ARF estimand and inverse probability weighting (IPW) estimator.")
titc4c <- c("ATE-ARF = Average treatment effect for the always-recruited in the favored cluster.")

if (wgt == 0) {
  cat(blank,titc1,titc2,titc3,blank, sep="\n")
  print(resg2, row.names = F)
  cat(blank,titc4a,titc4b,titc4c,blank, sep="\n")
} else {
  cat(blank,titc1,titc2,titl1a,titc3,blank, sep="\n")
  print(resg2, row.names = F)
  cat(blank,titc4a,titc4b,titc4c,blank, sep="\n")
}

# Print out results to file using sink()

sink(out_est,d_app)

if (wgt == 0) {
  cat(blank,titc1,titc2,titc3,blank, sep="\n")
  print(resg2, row.names = F)
  cat(blank,titc4a,titc4b,titc4c,blank, sep="\n")
} else {
  cat(blank,titc1,titc2,titl1a,titc3,blank, sep="\n")
  print(resg2, row.names = F)
  cat(blank,titc4a,titc4b,titc4c,blank, sep="\n")
}

sink()

}  # End Y Loop

} # End of function

