#****h* /simulation
#  NAME
#    simulation
#  FUNCTION 
#    Routines to generate simulated bivariate data. None of these are
#    user-facing, they were only used in simulations
#*******

#########################
### Main Simulation Routine
#########################

#****f* simulation/bivrecsim
#  NAME
#    bivrecsim --- conduct a set of simulations
#  FUNCTION
#    Conducts simulations and dumps the output to a file
#  INPUTS
#    Nsim       number of simulations
#    timedep    boolean, whether to use time-depencent covariates
#    dumpfile   file into which to dump the results
#    outputfrailties    boolean, whether to dump frailty estimates
#    type       type of frailties to generate
#    ragged     boolean, whether to use clusters of different sizes
#    censortime     fixed censoring time
#    ...        other parameters matching those of bivrec.agdata
#  OUTPUTS
#    a dump file containing simulation results
#  SYNOPSIS
bivrecsim <- function(Nsim, m = 10, Ji = rep(5, 10), K = 10, Kd = 10, timedep = FALSE,
        dumpfile = "dump", correction = "none", outputfrailties = FALSE,
        computesd = FALSE, type = "lognormal", dispest = "pearson", ragged = FALSE,
        censortime = 100, fixzero = NULL, smooth = FALSE, fullS = TRUE, boot = NULL,
        verbose = 2, maxiter = 200, alternating = FALSE)
#  SOURCE
#
{
        
    # Set global variables
    params <- initsetup(timedep, m, Ji, K, Kd, fixzero, alternating)
    
    # Simulation loop
     for(isim in 1:Nsim){
        
        cat("\nSimulation ", isim, "\n")
        if(ragged){
            Jilocal <- ceiling(runif(m) * max(Ji))
        }else{
            Jilocal <- Ji
        }

        ###########################
        #### Begin initialization of simulation i:
        
        ## Generate frailties, events, deaths, etc. There must be at least one event
        ## and one death in each cluster.
        nclustevents1 <- 0; nclustevents2 <- 0;ngens <- 0
        while(min(nclustevents1) == 0 | min(nclustevents2) == 0){
        
            #### Generate frailties:
            UV <- generatefrailty(type = type, params)
            Ui <- UV$Ui; Vi <- UV$Vi; Uij <- UV$Uij; Vij <- UV$Vij; 
            
            ## Generate a single covariate value (may or may not be constant)
            Zij <- generatecovariate(params)

            ## Compute followup times
            Cij <- generatecensoring(m, Jilocal, params$lambda0c, params$gamweibc)
            
            ## Generate recurrent event times
            Rij1 <- generaterecurrent(m, Jilocal, Zij, Uij, params$beta1,
                params$beta2, params$lambda0, params$gamweib, params$timedep, Cij)
            Rij2 <- generaterecurrent(m, Jilocal, Zij, Vij, params$beta1d,
                params$beta2d, params$lambda0d, params$gamweibd, params$timedep, Cij)
            
            ## Generate death times, based on covariates and frailty
            #Dij <- generatedeath(Zij, Rij, Vij, m, Jilocal)
            
            
            ## Compute the number of recurrent events observed
            nevents1 <- rep(0, sum(Jilocal))
            nevents2 <- rep(0, sum(Jilocal))
            for(ij in 1:sum(Jilocal)) nevents1[ij] <- sum(cumsum(Rij1[ij, ]) < Cij[ij])
            for(ij in 1:sum(Jilocal)) nevents2[ij] <- sum(cumsum(Rij2[ij, ]) < Cij[ij])
            ## Compute the number of events and deaths observed in each cluster
            ij <- rep(1:m, Jilocal)
            for (i in 1:m){
                nclustevents1[i] <- sum(nevents1[ij == i])
                nclustevents2[i] <- sum(nevents2[ij == i])
            }
            ngens <- ngens + 1
            
         #   if(alternating & any(nevents1 == 0)) nclustevents1 <- 0
            
        }
        
        ## Transform generated data into Anderson - Gill format
        agdata <- gen2AG(m, Jilocal, Zij, Rij1, Rij2, Cij, alternating, params$timedep)
            
        ## Clean up initialization
        #rm(Dij, Rij, Xij, nevents, nclustevents, nclustdeath, UV)
        
        ## Fit the model:
        fit <- bivrec(agdata, K, Kd, correction, dispest = dispest,
            computesd = computesd, fixzero = fixzero, smooth = smooth, fullS = fullS,
            verbose = verbose, maxiter = maxiter, alternating = alternating)
        if(!is.null(boot)){
            bootopts <- paste(names(boot), boot, sep = " = ", collapse = ", ")
            bootout <- NULL;
            bootstring <- paste("bootout <- bootrd(agdata, fit, ", bootopts, ")",
                sep = "")
            eval(parse(text = bootstring))
        }
        
        if(!is.null(names(fit))){
            ## DEBUG: Dump stuff to a file
            cat("Sim:\t", isim, "\tbetahat:\t", fit$regressionoutput$betahat, "\n" , 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)
            cat("Sim:\t", isim, "\tbetadhat:\t ", fit$regressionoutput$betadhat, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)
            if(computesd){
            cat("Sim:\t", isim, "\tstdbetahat:\t ", fit$stderr[1], "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)
            cat("Sim:\t", isim, "\tstdbetadhat:\t ", fit$stderr[2], "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)
            }
            if(outputfrailties){ 
            cat("Sim:\t", isim, "\tUi:\t ", Ui, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)  
            cat("Sim:\t", isim, "\tUihat:\t ", fit$frailtyoutput$Uihat, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)  
            cat("Sim:\t", isim, "\tVi:\t ", Vi, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)  
            cat("Sim:\t", isim, "\tVihat:\t ", fit$frailtyoutput$Vihat, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)  
            cat("Sim:\t", isim, "\tUij:\t ", Uij, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)  
            cat("Sim:\t", isim, "\tUijhat:\t ", fit$frailtyoutput$Uijhat, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)  
            cat("Sim:\t", isim, "\tVij:\t ", Vij, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)  
            cat("Sim:\t", isim, "\tVijhat:\t ", fit$frailtyoutput$Vijhat, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)      
            }
            cat("Sim:\t", isim, "\tsigma2hat:\t ", fit$dispparams$sigma2hat, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)      
            cat("Sim:\t", isim, "\tsigma2dhat:\t ", fit$dispparams$sigma2dhat, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)      
            cat("Sim:\t", isim, "\tnu2hat:\t ", fit$dispparams$nu2hat, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)      
            cat("Sim:\t", isim, "\tnu2dhat:\t ", fit$dispparams$nu2dhat, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)    
            cat("Sim:\t", isim, "\tthetahat:\t ", fit$dispparams$thetahat, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)      
 
           # Dump initial values to file
            cat("Sim:\t", isim, "\tbetahat:\t", fit$initial$betahat, "\n" ,
				file = paste(dumpfile, ".init.txt", sep = ""), append = TRUE)
            cat("Sim:\t", isim, "\tbetadhat:\t ", fit$initial$betadhat, "\n", 
				file = paste(dumpfile, ".init.txt", sep = ""), append = TRUE)           
            cat("Sim:\t", isim, "\tsigma2hat:\t ", fit$initial$dispparams$sigma2hat, "\n", 
				file = paste(dumpfile, ".init.txt", sep = ""), append = TRUE)      
            cat("Sim:\t", isim, "\tsigma2dhat:\t ", fit$initial$dispparams$sigma2dhat, "\n", 
				file = paste(dumpfile, ".init.txt", sep = ""), append = TRUE)      
            cat("Sim:\t", isim, "\tnu2hat:\t ", fit$initial$dispparams$nu2hat, "\n", 
				file = paste(dumpfile, ".init.txt", sep = ""), append = TRUE)      
            cat("Sim:\t", isim, "\tnu2dhat:\t ", fit$initial$dispparams$nu2dhat, "\n", 
				file = paste(dumpfile, ".init.txt", sep = ""), append = TRUE)    
            cat("Sim:\t", isim, "\tthetahat:\t ", fit$initial$dispparams$thetahat, "\n", 
				file = paste(dumpfile, ".init.txt", sep = ""), append = TRUE)      
 
 
        }
        if(!is.null(boot)) {
            if(!is.null(names(bootout))){
            cat("Sim:\t", isim, "\tboot.betahat:\t", bootout$mean$betahat, "\n" ,
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)
            cat("Sim:\t", isim, "\tboot.betadhat:\t ", bootout$mean$betadhat, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)
            cat("Sim:\t", isim, "\tboot.stdbetahat:\t ", bootout$sd$betahat, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)
            cat("Sim:\t", isim, "\tboot.stdbetadhat:\t ", bootout$sd$betadhat, "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)
            cat("Sim:\t", isim, "\tboot.sigma2hat:\t ", bootout$mean$dispparams["sigma2hat"], "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)      
            cat("Sim:\t", isim, "\tboot.sigma2dhat:\t ", bootout$mean$dispparams["sigma2dhat"], "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)      
            cat("Sim:\t", isim, "\tboot.nu2hat:\t ", bootout$mean$dispparams["nu2hat"], "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)      
            cat("Sim:\t", isim, "\tboot.nu2dhat:\t ", bootout$mean$dispparams["nu2dhat"], "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)    
            cat("Sim:\t", isim, "\tboot.thetahat:\t ", bootout$mean$dispparams["thetahat"], "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)      
            cat("Sim:\t", isim, "\tboot.stdsigma2hat:\t ", bootout$sd$dispparams["sigma2hat"], "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)      
            cat("Sim:\t", isim, "\tboot.stdsigma2dhat:\t ", bootout$sd$dispparams["sigma2dhat"], "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)      
            cat("Sim:\t", isim, "\tboot.stdnu2hat:\t ", bootout$sd$dispparams["nu2hat"], "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)      
            cat("Sim:\t", isim, "\tboot.stdnu2dhat:\t ", bootout$sd$dispparams["nu2dhat"], "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)    
            cat("Sim:\t", isim, "\tboot.stdthetahat:\t ", bootout$sd$dispparams["thetahat"], "\n", 
				file = paste(dumpfile, ".txt", sep = ""), append = TRUE)      
            }
        }

    }
    
}
#************ bivrecsim 

#****f* simulation/initsetup
#  NAME
#    initsetup  --- intial setup for simulation
#  FUNCTION
#    Creates parameters for simulation
#  SYNOPSIS
initsetup <- function(timedep, m, Ji, K, Kd, fixzero, alternating = F)
#  SOURCE
#
{
    m <- m
    Ji <- Ji
    if(is.null(fixzero)) fixzero <- ""
    # Variance of cluster frailty (rec)
    if(fixzero != "sigma2hat") sigma2 <- .25  else sigma2 <- 0 
    # Variance of cluster frailty (death)
    if(fixzero != "sigma2dhat") sigma2d <- .25 else sigma2d <- 0 
    # Additional variance for individual frailty (rec)
    if(fixzero != "nu2hat") nu2 <- .25 else nu2 <- 0 
    # Additional variance for individual frailty (death)
    if(fixzero != "nu2dhat") nu2d <- .25 else nu2d <- 0 
    # Covariance
    if(fixzero != "thetahat") theta <- .125 else theta <- 0 
    beta1d <- 1 # Coefficient for the death process
    beta1 <- 1 # Coefficient for the recurrent process
    beta2d <- 1 # Coefficient for the death process
    beta2 <- 1 # Coefficient for the recurrent process
    lambda0 <- 10 # Baseline hazard for recurrent events
    lambda0d <- 10 # Baseline hazard for death
    lambda0c <- 1 # Baseline hazard for censoring
    if(alternating) lambda0c <- .25
    gamweib <- 1.8 # Shape parameter for Weibull distribution (rec)
    gamweibd <- 1.8 # Shape parameter for Weibull distribution (death)
    gamweibc <- 1.8 # Shape parameter for Weibull distribution (censoring)
    Z1mean <- 0 # Mean of the Zij covariate (generation)
    Z1var <- .5 # Variance of the Zij covariate (for generation)
    Z2mean <- 0 # Mean of the Zij covariate (generation)
    Z2var <- .5 # Variance of the Zij covariate (for generation)
    maxiter.outer <- 200
    maxiter.inner <- 50
    thresh.outer <- 1e-3
    thresh.inner <- 1e-5
    timedep <- timedep
    K <- K # Number of breaks in the baseline hazard
    Kd <- Kd # Number of breaks in the baseline hazard for death
    params = list(m = m, Ji = Ji, sigma2 = sigma2, sigma2d = sigma2d, 
            nu2 = nu2, nu2d = nu2d, theta = theta,
            beta1d = beta1d, beta2d = beta2d, beta1 = beta1, beta2 = beta2,
            lambda0 = lambda0, lambda0d = lambda0d, lambda0c = lambda0c,
            gamweib = gamweib, gamweibd = gamweibd, gamweibc = gamweibc,
            Z1mean = Z1mean, Z2mean = Z2mean, Z1var = Z1var, Z2var = Z2var,
            maxiter.outer = maxiter.outer, maxiter.inner = maxiter.inner,
            thresh.outer = thresh.outer, thresh.inner = thresh.inner,
            timedep = timedep, K = K, Kd = Kd)
}
#************ initsetup 

#****f* simulation/L
#  NAME
#    L --- stratification function
#  FUNCTION
#    Simulates strata
#  SYNOPSIS
L <- function(i, j, k, Z)
#  SOURCE
#
{
    # Simple stratification: First event, stratum 1, events 2 - 3, stratum 2,
    # after that, stratum 3
    if(k == 1){out <- 1}
    if(k > 1 & k < 4){out <- 2}
    if(k > 3){out <- 3}
    out <- 1
    return(out)
}
#************ L 

#****f* simulation/gen2AG
#  NAME
#    gen2AG
#  FUNCTION
#    Convert generated data into Anderson - Gill format
#  INPUTS
#    m      number of clusters
#    Ji     cluster size
#    Zij    list of generated covariates and times
#    Rij1   recurrent event gap times 1
#    Rij2   recurrent event gap times 2
#    Cij    censoring times
#    alternating  boolean epsisodic data indicator
#    timedep      time-dependent covariate indicator
#  OUTPUTS
#    data frame in format used by bivrec.agdata
#  SYNOPSIS
gen2AG <- function(m, Ji, Zij, Rij1, Rij2, Cij, alternating, timedep){
#  SOURCE
#
    if(timedep){
        Zij1 <- Zij$Zij1
        Zij1times <- Zij$Zij1times
        Zij2 <- matrix(0, dim(Zij1)[1], 100)
    }else{
        if(!is.matrix(Zij) && !is.list(Zij)) {
            Zij1 <- matrix(rep(Zij, 100), length(Zij), 100)
            Zij2 <- matrix(0, dim(Zij1)[1], 100)
        }
        if(is.list(Zij)){
            Zij1 <- matrix(rep(Zij$Zij1, 100), length(Zij$Zij1), 100)
            Zij2 <- matrix(rep(Zij$Zij2, 100), length(Zij$Zij2), 100)
        }    
        Zij1times <- matrix(Inf, dim(Zij1)[1], 100)
    }
	kmax <- 100 # maximum number of events allowed before censoring
    ## Compute the number of observed events for each individual
	nevents1 <- rep(0, sum(Ji))  
    ## Compute the number of observed events for each individual
	nevents2 <- rep(0, sum(Ji))  
	ncovchange <- rep(0, sum(Ji))
    ## observed events
	for(ij in 1:sum(Ji)) nevents1[ij] <- sum(cumsum(Rij1[ij, ]) < Cij[ij]) 
	for(ij in 1:sum(Ji)) nevents2[ij] <- sum(cumsum(Rij2[ij, ]) < Cij[ij]) 
    ## observed covariate changes
	for(ij in 1:sum(Ji)) ncovchange[ij] <- sum(cumsum(Zij1times[ij, ]) < Cij[ij]) 
	outdata.matrix <- matrix(0, sum(nevents1 + nevents2 + ncovchange + 1), 10)
	outrow <- 1
	# Loop runs over all individuals and generates an Anderson - Gill line
	# for each event time that is less than the follow - up time
	# Columns: i: Cluster, j: Individual, k: At risk for which event
    #   start / stop: Cumulative times, Z: Covariate values, Delta / delta: indicators 
    if(!alternating){
        for(i in 1:m){
            for(j in 1:Ji[i]){
                ij <- sum(Ji[0:(i - 1)]) + j
                ThisFollowTime <- Cij[ij]
                CumTimes1 <- cumsum(Rij1[ij, ])
                CumTimes2 <- cumsum(Rij2[ij, ])
                CumTimesZ <- cumsum(Zij1times[ij, ])
                CumTimes1 <- CumTimes1[CumTimes1 < ThisFollowTime]
                CumTimes2 <- CumTimes2[CumTimes2 < ThisFollowTime]
                CumTimesZ <- CumTimesZ[CumTimesZ < ThisFollowTime]
                l1 <- length(CumTimes1)
                l2 <- length(CumTimes2)
                lZ <- length(CumTimesZ)
                thisoutdata <- matrix(0, l1 + l2 + lZ + 1, 10)
                thisoutdata[, 1] <- i;thisoutdata[, 2] <- j;
                thisoutdata[, 6] <- c(CumTimes1, CumTimes2, CumTimesZ, ThisFollowTime)
                thisoutdata[, 7] <- c(rep(1, l1), rep(0, l2 + lZ + 1))
                thisoutdata[, 8] <- c(rep(0, l1), rep(1, l2), rep(0, lZ + 1))
                thisoutdata[, 9] <- c(rep(NA, l1 + l2), Zij1[ij, (lZ > 0) * (1:lZ)],
                    Zij1[ij, min(which(cumsum(Zij1times[ij, ]) > ThisFollowTime))])
                thisoutdata[, 10] <- c(rep(NA, l1 + l2), Zij2[ij, (lZ > 0) * (1:lZ)],
                    Zij2[ij, min(which(cumsum(Zij1times[ij, ]) > ThisFollowTime))])
                thisoutdata <- matrix(thisoutdata[order(thisoutdata[, 6]), ],
                    dim(thisoutdata)[1], dim(thisoutdata)[2])    
                thisoutdata[, 3] <- cumsum(c(1, thisoutdata[, 7] + 
                    thisoutdata[, 8]))[ - (dim(thisoutdata)[1] + 1)]     
                for(ind in (dim(thisoutdata)[1]):1){
                    if(!is.na(thisoutdata[ind, 9])) thisZ1 <- thisoutdata[ind, 9] 
                    else thisoutdata[ind, 9] <- thisZ1
                    if(!is.na(thisoutdata[ind, 10])) thisZ2 <- thisoutdata[ind, 10] 
                    else thisoutdata[ind, 10] <- thisZ2
                    k <- thisoutdata[ind, 3]
                    thisoutdata[ind, 4] <- L(i, j, k, Zij1[ij, k])
                    if(ind > 1) thisoutdata[ind, 5] <- thisoutdata[ind - 1, 6]
                }
                newoutrow <- outrow + dim(thisoutdata)[1]
                outdata.matrix[outrow:(newoutrow - 1), ] <- thisoutdata
                outrow <- newoutrow
            }
        }
    }else{
        for(i in 1:m){
            for(j in 1:Ji[i]){
                ij <- sum(Ji[0:(i - 1)]) + j
                ThisFollowTime <- Cij[ij]
                CumTime <- 0
                k <- 1;eventct <- 1
                while(CumTime < ThisFollowTime){
                    if((CumTime + Rij1[ij, eventct]) < ThisFollowTime){
                        thisoutdata <- c(i, j, k, L(i, j, k, Zij1[ij, eventct]),
                            CumTime, CumTime + Rij1[ij, eventct], 1, 0, 
                            Zij1[ij, eventct], Zij2[ij, eventct])
                        CumTime <- CumTime + Rij1[ij, eventct]
                        outdata.matrix[outrow, ] <- thisoutdata
                        k <- k + 1;outrow <- outrow + 1
                    }else{
                        thisoutdata <- c(i, j, k, L(i, j, k, Zij1[ij, eventct]),
                            CumTime, ThisFollowTime, 0, 0, Zij1[ij, eventct],
                            Zij2[ij, eventct])
                        CumTime <- ThisFollowTime
                        outdata.matrix[outrow, ] <- thisoutdata
                        k <- k + 1;outrow <- outrow + 1
                    }
                    if((CumTime + Rij2[ij, eventct]) < ThisFollowTime){
                        thisoutdata <- c(i, j, k, L(i, j, k, Zij2[ij, eventct]),
                            CumTime, CumTime + Rij2[ij, eventct], 0, 1,
                            Zij1[ij, eventct], Zij2[ij, eventct])
                        CumTime <- CumTime + Rij2[ij, eventct]
                        outdata.matrix[outrow, ] <- thisoutdata
                        k <- k + 1;outrow <- outrow + 1
                    }else{
                        if(CumTime < ThisFollowTime){
                            thisoutdata <- c(i, j, k, L(i, j, k, Zij2[ij, eventct]),
                                CumTime, ThisFollowTime, 0, 0, Zij1[ij, eventct],
                                Zij2[ij, eventct])
                            CumTime <- ThisFollowTime
                            outdata.matrix[outrow, ] <- thisoutdata
                            k <- k + 1;outrow <- outrow + 1
                        }
                    }
                    eventct <- eventct + 1
                }
            }
        }
    
    }
	outdata <- as.data.frame(outdata.matrix[1:(outrow - 1), ])
	colnames(outdata) <- c("i", "j", "k", "r", "start", "stop",
                            "delta", "Delta", "Z1", "Z2")
	return(outdata)
}
#************ gen2AG 

#****f* simulation/generatefrailty
#  NAME
#    generatefrailty --- simulate frailties
#  FUNCTION
#    Generates simulated frailties by different methods depending on the
#    distribution chosen.
#  INPUTS
#    type    distribution string
#    params  cluster information and dispersion parameters
#  OUTPUTS
#    Ui, Uij, Vi, Vij vectors
#  SYNOPSIS
generatefrailty <- function(type = "lognormal", params)
#  SOURCE
#
{
    Ui <- -1;Vi <- -1;Uij <- -1;Vij <- -1;
    m <- params$m;Ji <- params$Ji
    sigma2 <- params$sigma2;sigma2d <- params$sigma2d;
    nu2 <- params$nu2; nu2d <- params$nu2d; theta <- params$theta
    # Occasionally, trivariate reduction or gaussian frailties can be negative, 
    # continue generating until this is no longer the case.
    while(sum(Ui < 0) + sum(Vi < 0) + sum(Uij < 0) + sum(Vij < 0) > 0){
        Ui <- rep(0, m); Vi <- rep(0, m)
        Uij <- rep(0, sum(Ji)); Vij <- rep(0, sum(Ji))
        # Cumulative sum of number of cluster individuals
        Jicum <- c(1, cumsum(Ji) + 1) 
        if(type == "gamma"){
            #! Gamma subject - level frailties are generated by trivariate reduction, as follows:
             for (i in 1:m){
                # Ensure that conditions for trivariate reduction are met
                while(Ui[i] <= theta / nu2d | Vi[i] <= theta / nu2){
                    # Mean 1, variance sigma2
                    Ui[i] <- rgamma(1, shape = 1 / sigma2, scale = sigma2)  
                    # Mean 1, variance sigma2d
                    Vi[i] <- rgamma(1, shape = 1 / sigma2d, scale = sigma2d)  
                }
                Y1 <- rgamma(Ji[i], shape = Ui[i] / nu2 - theta / (nu2 * nu2d), scale = 1)
                Y2 <- rgamma(Ji[i], shape = Vi[i] / nu2d - theta / (nu2 * nu2d), scale = 1)
                Y3 <- rgamma(Ji[i], shape = theta / (nu2 * nu2d), scale = 1)
                if(theta == 0){ Y3 <- 0 }
                Uij[Jicum[i]:(Jicum[i + 1] - 1)] <- nu2 * (Y1 + Y3)
                Vij[Jicum[i]:(Jicum[i + 1] - 1)] <- nu2d * (Y2 + Y3)
            }
        }
        if(type == "gaussian"){
            #!For Gaussian frailties, the cluster - level frailties are generated as iid Normals,
            #!and given these, the subject - level frailties are generated with the covariance matrix
            #!above. The code checks that the covariance matrix is positive definite, which is not
            #!always assured.
            for (i in 1:m){
                eig <- -1
                while(sum(eig < 0) > 0){
                    Ui[i] <- rnorm(1, mean = 1, sd = sqrt(sigma2))
                    Vi[i] <- rnorm(1, mean = 1, sd = sqrt(sigma2d))
                    covmat <- matrix(c(Ui[i] * nu2, theta, theta, Vi[i] * nu2d), 2, 2)
                    eig <- eigen(covmat)$values
                    if(sum(eig < 0) > 0) print("hello")
                }
                UV <- mvrnorm(Ji[i], c(Ui[i], Vi[i]), covmat)
                Uij[Jicum[i]:(Jicum[i + 1] - 1)] <- UV[, 1]
                Vij[Jicum[i]:(Jicum[i + 1] - 1)] <- UV[, 2]
            }
        }
        if(type == "lognormal"){
            for (i in 1:m){
                eig <- -1
                while(sum(eig < 0) > 0){
                    # Mean 1, variance sigma2                
                    Ui[i] <- exp(rnorm(1, log(1 / sqrt(1 + sigma2)), sqrt(log(1 + sigma2))) ) 
                    # Mean 1, variance sigma2d               
                    Vi[i] <- exp(rnorm(1, log(1 / sqrt(1 + sigma2d)), sqrt(log(1 + sigma2d))) ) 
                    muprime <- c(log(Ui[i]^2 / sqrt(Ui[i]^2 + Ui[i] * nu2)),
                        log(Vi[i]^2 / sqrt(Vi[i]^2 + Vi[i] * nu2d)))
                    covmat <- matrix(0, 2, 2)
                    covmat[1, 1] <- log(1 + Ui[i] * nu2 / Ui[i]^2)
                    covmat[1, 2] <- log(1 + theta / Ui[i] / Vi[i])
                    covmat[2, 1] <- log(1 + theta / Ui[i] / Vi[i])
                    covmat[2, 2] <- log(1 + nu2d * Vi[i] / Vi[i]^2)
                    eig <- eigen(covmat)$values
                }
                UV <- matrix(mvrnorm(Ji[i], muprime, covmat), Ji[i], 2)
                Uij[Jicum[i]:(Jicum[i + 1] - 1)] <- exp(UV[, 1])
                Vij[Jicum[i]:(Jicum[i + 1] - 1)] <- exp(UV[, 2])
            }
        }
        if(type == "lognormperm"){
            genlognormal <- function(n, mu, sigma2){
                sigma2prime <- log(1 + sigma2 / mu^2)        
                muprime <- log(mu) - 1 / 2 * sigma2prime
                return(rlnorm(n, meanlog = muprime, sdlog = sqrt(sigma2prime)))            
            }
            correlate <- function(x, y, covar){
                sorttop <- function(x, pct){
                     nsort <- round(length(x) * pct)
                     if(nsort < 2) return(x)
                     if(nsort == length(x)) return(sort(x))
                     return(c(sort(x[1:nsort]), x[(nsort + 1):length(x)]))
                }
                if(length(x) == 1) return(list(x = x, y = y))
                pctopt <- optimize(function(pct, target, x, y) (cov(sorttop(x, pct), sorttop(y, pct)) - target)^2, interval = c(0, 1), target = covar, x = x, y = y)$minimum    
                cat(" ", pctopt)
                return(list(x = sorttop(x, pctopt), y = sorttop(y, pctopt)))
            }
            for(i in 1:m){
                Ui[i] <- exp(rnorm(1, log(1 / sqrt(1 + sigma2)), sqrt(log(1 + sigma2))) ) 
                Vi[i] <- exp(rnorm(1, log(1 / sqrt(1 + sigma2d)), sqrt(log(1 + sigma2d))) ) 
                Uijprop <- genlognormal(Ji[i], Ui[i], Ui[i] * nu2)
                Vijprop <- genlognormal(Ji[i], Vi[i], Vi[i] * nu2d)
                corrUV <- correlate(Uijprop, Vijprop, theta)
                Uij[Jicum[i]:(Jicum[i + 1] - 1)] <- corrUV$x
                Vij[Jicum[i]:(Jicum[i + 1] - 1)] <- corrUV$y            
            }
        }
                   
    }
    return(list(Ui = Ui, Vi = Vi, Uij = Uij, Vij = Vij))
}
#************ generatefrailty 

#****f* simulation/generatecovariate
#  NAME
#    generatecovariate
#  FUNCTION
#    Generates a single covariate, Normal with a given mean and variance.
#    Old code to generate multiple covariates or time - dependent covariates 
#    has been removed.
#  SYNOPSIS
generatecovariate <- function(params)
#  SOURCE
#
{
    m <- params$m;Ji <- params$Ji;timedep <- params$timedep
    Z1mean <- params$Z1mean;Z1var <- params$Z1var
    gamweib <- params$gamweib;lambda0 <- params$lambda0;
    if(timedep){
        Zij1 <- rnorm(sum(Ji) * 100, Z1mean, Z1var);dim(Zij1) <- c(sum(Ji), 100)
        Zij1times <- rweibull(sum(Ji) * 100, shape = gamweib,
            scale = lambda0^(-1 / gamweib));dim(Zij1times) <- c(sum(Ji), 100)
        return(list(Zij1 = Zij1, Zij1times = Zij1times))
    }else{
        Zij1 <- rnorm(sum(Ji), Z1mean, Z1var)
        return(Zij1)
    }
}
#************ generatecovariate 


#****f* simulation/generaterecurrent
#  NAME
#    generaterecurrent
#  FUNCTION
# Generates recurrent event times from a Weibull hazard with baseline lambda_0, shape gamma.
#  INPUTS
#    m      number of clusters
#    Ji     cluster sizes
#    Zij    covariates
#    Uij    frailties
#    beta1  coefficient for covariate 1
#    beta2  coefficient for covariate 2 (if applicable)
#    labmda0    weibull baseline
#    gamweib    weibull shape
#    timedep    boolean to indicate time-dependence
# SYNOPSIS
generaterecurrent <- function(m, Ji, Zij, Uij, beta1, beta2, lambda0, gamweib, timedep, Cij)
{
#  SOURCE
#
    # Convert the provided covariates into two matrices of time - dependent covariates
    if(timedep){
        Zij1 <- Zij$Zij1
        Zij1times <- Zij$Zij1times
        Zij2 <- matrix(0, dim(Zij1)[1], 100)  
        
        Rij <- rep(0, sum(Ji) * 100); dim(Rij) <- c(sum(Ji), 100)
        for(ind in 1:sum(Ji))
        {
            timescum <- cumsum(Zij1times[ind, ])
            timescum[min(which(timescum > Cij[ind]))] <- Inf
            obstimes <- c(timescum[timescum < Cij[ind]], Cij[ind])
            Rijcum <- 0;k <- 1
            while(k < 100){
                
                hazards <- Uij[ind] * lambda0 * gamweib * (obstimes - Rijcum)^(gamweib - 1) *
                    exp(beta1 * Zij1[ind, 1:length(obstimes)] + beta2 *
                    Zij2[ind, 1:length(obstimes)])
                maxhaz <- max(hazards, na.rm = TRUE)
                accept <- FALSE
                # Generate by accept-reject
                n <- 0;Rijprop <- 0
                while(!accept){
                    genunif <- runif(1)
                    Rijprop <- Rijprop - 1 / maxhaz * log(genunif)
                    thistime <- min(which(timescum > (Rijcum + Rijprop)))
                    if(thistime == 101) thistime <- 100
                    thishaz <- Uij[ind] * lambda0 * gamweib * (Rijprop)^(gamweib - 1) *
                        exp(beta1 * Zij1[ind, thistime] + beta2 * Zij2[ind, thistime])
                    accunif <- runif(1)
                    if(accunif * maxhaz < thishaz){
                        Rij[ind, k] <- Rijprop
                        Rijcum <- Rijcum + Rijprop
                        if(Rijcum > Cij[ind]) k <- 100
                        accept <- TRUE
                    }
                    n <- n + 1
                }
            k <- k + 1
            }   
        }
        return(Rij)     
    }else{
        if(!is.matrix(Zij) && !is.list(Zij)) {
            Zij1 <- matrix(rep(Zij, 100), length(Zij), 100)
            Zij2 <- matrix(0, dim(Zij1)[1], 100)
        }
        if(is.list(Zij)){
            Zij1 <- matrix(rep(Zij$Zij1, 100), length(Zij$Zij1), 100)
            Zij2 <- matrix(rep(Zij$Zij2, 100), length(Zij$Zij2), 100)
        }    
        # Initialize output matrix
        Rij <- rep(0, sum(Ji) * 100)  # Gap data for times between events
        dim(Rij) <- c(sum(Ji), 100)
        # Generate interevent gap times by inversion
        for(k in 1:100)  Rij[, k] <- ((-exp(-beta1 * Zij1[, k] - beta2 * Zij2[, k]) *
            log(1 - runif(sum(Ji))) / (lambda0 * Uij))^(1 / gamweib))
        Rij[Rij == 0] <- 1
        return(Rij)
    }
}
#************ generaterecurrent 


#****f* simulation/generatecensoring
#  NAME
#    generatecensoring --- generate censoring times
#  SYNOPSIS
generatecensoring <- function(m, Ji, lambda0c, gamweibc)
#  SOURCE
#
{
    ### Random censoring
    C <- rweibull(sum(Ji), shape = gamweibc, scale = lambda0c^(-1 / gamweibc))    
}
#************ generatecensoring 
