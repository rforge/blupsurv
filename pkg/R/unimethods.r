#****h* /methodsUni
#  NAME
#    methodsUni --- methods for univariate recurrent event processes
#  FUNCTION
#    S3 methods for univariate recurrent events
#*******

#****f* methodsUni/unirec
#  NAME
#    unirec --- S3 generic
#  FUNCTION
#    Dispatches to either unirec.formula, or unirec.data.frame 
#  SYNOPSIS
unirec <- function(x, ...)
{
#  SOURCE
#
    UseMethod("unirec")
}

#************ unirec 

#****f* methodsUni/unirec.data.frame
unirec.data.frame <- function(x, ...)
#  NAME
#    unirec.data.frame --- method for data frames 
#  FUNCTION
#    Method for data frames in in the required format, see bivrec.agdata for detail.
#    It works by adding extra columns to the input and passing it into the bivariate
#    functions.
#  INPUTS
#    x      a data frame
#  OUTPUTS
#           an object of class unirec
#  SYNOPSIS
{
#  SOURCE
#
    call <- match.call()
    if(dim(x)[2] > 7 && all(colnames(x)[1:7] == c("i", "j", "k", "r", "start",
        "stop", "delta"))) {
        x <- x[order(x$i, x$j, x$k), ]
        varnames <- colnames(x)[ - (1:7)]
        # Add some extra empty columns, so it can be passsed to bivrec.agdata
        x <- cbind(x[, 1:7], 0, x[, -(1:7)])
        colnames(x)[ - (1:7)] <- c("Delta", varnames)
        out <- bivrec.agdata(x, fixzero = c("clust2", "subj2", "cov"), ...)
        out$call <- call
        return(out)
    } else {
        stop("input data frame needs to have colnames i, j, k, r, start, stop, delta")
    }   
}
#************ unirec.data.frame 

#****f* methodsUni/unirec.formula
#  NAME
#    unirec.formula --- method for formulas
#  FUNCTION
#    This is the main user-facing method of the package for univariate data
#  INPUTS
#    See package documentation
#  OUTPUTS
#    An object of class unirec. See package documentation
#  SYNOPSIS
unirec.formula <- function(formula, data = parent.frame(), ...)
{
#  SOURCE
#
    # in part based on coxph function
    # Some of the code here is duplicated from bivrec.formula, but I don't see
    # a good way around it
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, sys.parent()))) m$data <- as.data.frame(data)
    m$...<-NULL
    m[[1]] <- as.name("model.frame")
    special <- c("id", "cluster", "strata")
    Terms <- if (missing(data)) terms(formula, special) 
        else terms(formula, special, data = data)    
    m$formula <- Terms
    m <- eval(m, sys.parent())
    n <- nrow(m)
    resp <- model.extract(m, "response")
    if (!is.Surv(resp)) stop("model response must be a Surv object")
    if(dim(resp)[2] != 3) stop("response must be of type Surv(start, stop, status)")
    start <- resp[, 1]
    stop <- resp[, 2]
    delta <- resp[, 3]
    dropx <- NULL
    clusterind <- attr(Terms, "specials")$cluster

    # Cluster handling
    clusternames <- NULL
    if(length(clusterind) > 0){
        cluster <- m[, clusterind]
        if(is.factor(cluster)) clusternames <- levels(cluster)
        if(is.numeric(cluster)) clusternames <- as.character(unique(cluster))
        i <- as.numeric(as.factor(cluster))
        tempc <- untangle.specials(Terms, "cluster", 1:10)
        ord <- attr(Terms, "order")[tempc$terms]
        if (any(ord > 1)) stop("Cluster cannot be used in an interaction")
        dropx <- c(dropx, tempc$terms)
    }else{
        cluster <- rep(1, n)
        i <- rep(1, n)
    }

    # ID handling
    idind <- attr(Terms, "specials")$id
    if(length(idind) == 0) stop("a subject id (unique within clusters) is required")
    else{
        id <- m[, idind]
        subjnames <- paste(cluster, id, sep = ".")
        j <- rep(0, n)
        for(ii in 1:max(i)){
            iids <- unique(id[i == ii])
            jj <- 1
            for(thisid in iids){
                j[i == ii & id == thisid] <- jj
                jj <- jj + 1
            }
        }
        tempi <- untangle.specials(Terms, "id", 1:10)
        ord <- attr(Terms, "order")[tempi$terms]
        if (any(ord > 1)) stop("id cannot be used in an interaction")
        dropx <- c(dropx, tempi$terms)
    }

    # Stratum handling
    stratind <- attr(Terms, "specials")$strata
    if(length(stratind) > 0) 
    {
        strat <- m[, stratind]
        ustrat <- unique(strat)
        stratnames <- as.character(ustrat)
        r <- rep(0, n)
        for(rr in 1:length(ustrat)){
            r[strat == ustrat[rr]] <- rr
        }
        tempr <- untangle.specials(Terms, "strata", 1:10)
        ord <- attr(Terms, "order")[tempr$terms]
        if (any(ord > 1)) stop("strata cannot be used in an interaction")
        dropx <- c(dropx, tempr$terms)
    }else{
        r <- rep(1, n)
        stratnames <- "1"
    }
    
    # Construct data frame to pass to unirec.agdata
    Ki <- table(i * 1e6 + j)
    k <- unlist(as.vector(sapply(Ki, function(x) 1:x)))
    newTerms <- if(length(dropx))  Terms[ - dropx] else Terms
    X <- model.matrix(newTerms, m)
    X <- X[, -1, drop = FALSE]
    agdata <- as.data.frame(cbind(i, j, k, r, start, stop, delta, X))
    sortord <- order(agdata$i, agdata$j)
    agdata <- agdata[sortord, ]
    subjnames <- unique(subjnames[sortord])
    processname <- as.character(call$formula[[2]][[4]])
    fit <- unirec(agdata, ...)
    fit$call <- call
    fit <- cleanunirecoutput(fit, clusternames, subjnames, stratnames, processname)
    return(fit)
}
#************ unirec.formula 


#****f* methodsUni/cleanunirecoutput
#  NAME
#    cleanunirecoutput --- clean output of unirec.agdata
#  FUNCTION
#    Construct the unirec object and 
#    restore the original stratum, process, cluster and subject names that
#    were stripped by unirec.formula
#  INPUTS
#    fit            a fit returned by unirec.agdata
#    clusternames   vector of original cluster IDs
#    subjnames      vector of original subject IDs
#    stratnames     vector of original stratum IDs
#    processnames   vector of original process IDs
#  OUTPUTS
#    an object of class unirec, with properly named components
#  SYNOPSIS
cleanunirecoutput <- function(fit, clusternames, subjnames, stratnames, processname)
#  SOURCE
#
{
    out <- NULL
    out$call <- fit$call
    # regression
    out$regression = list(coefficients = fit$regressionoutput$betahat,
                        loglik = fit$regressionoutput$loglik1)
    # frailties
    clust <- fit$frailtyoutput$Uihat
    subj <- fit$frailtyoutput$Uijhat
    names(clust) <- clusternames
    names(subj) <- subjnames
    out$frailty = list(clust = clust, subj = subj)
    # dispersion
    disp <- list(clust = fit$dispparams$sigma2hat, subj = fit$dispparams$nu2hat)
    out$dispersion <- disp
    # baseline
    breaks <- fit$regressionoutput$as
    basehaz <- fit$regressionoutput$alphars
    rownames(breaks) <- rownames(basehaz) <- stratnames
    out$hazard = list(breaks = breaks, hazard = basehaz)
    # summaries
    summary.reg <- fit$summary.reg[, 1:3]
    colnames(summary.reg) <- c("coef", "sd", "pval")
    summary.disp <- fit$summary.disp[c(1, 3)]
    names(summary.disp) <- c("clust", "subj")
    out$summaries <- list(regression = summary.reg, dispersion = summary.disp)
    attr(out, "processname") <- processname
    class(out) <- "unirec"
    return(out)
}
#************ cleanunirecoutput 

######### Print and summary methods

#****f* methodsUni/print.unirec
#  NAME
#    print.unirec  --- S3 print method
#  FUNCTION
#    Prints unirec objects
#  INPUTS
#    x      object of class unirec
#  SYNOPSIS
print.unirec <- function(x, ...)
#  SOURCE
#
{
    processname <- attr(x, "processname")
    cat("Process ", processname, ":\n", sep = "")
    print(x$regression$coefficients, ...)
    invisible(x)
}
#************ print.unirec 


#****f* methodsUni/summary.unirec
#  NAME
#    summary.unirec --- S3 summary method
#  FUNCTION
#    Computes a summary for the unirec class
#  INPUTS
#    object     an object of class unirec
#    digits     precision to print the output
#  OUTPUTS
#    an object of class summary.unirec
#  SYNOPSIS
summary.unirec <- function(object, digits = 4, ...)
#  SOURCE
#
{
    x <- object
    out <- NULL
    dots <- as.list(substitute(list(...)))[ - 1]
    sreg <- format(x$summaries$regression, digits = digits, ...)
    sdisp <- format(x$summaries$dispersion, digits = digits, ...)
    out <- list(call = x$call, summary.reg = sreg, summary.disp = sdisp)
    processname <- attr(x, "processname")
    attr(out, "processname") <- processname
    class(out) <- "summary.unirec"   
    return(out)
}
#************ summary.unirec 


#****f* methodsUni/print.summary.unirec
#  NAME
#    print.summary.unirec -- S3 method to print summary.unirec objects
#  INPUTS
#    x      an objet of class summary.unirec
#  SYNOPSIS
print.summary.unirec <- function(x, ...)
#  SOURCE
#
{
    processname <- attr(x, "processname")
    dots <- x$dots
    sreg <- x$summary.reg
    sdisp <- x$summary.disp
    colnames(sreg)[1:3] <- paste(processname, c("coef", "sd", "pval"), sep = ".")
    cat("\nSummary of regression coefficients: \n")
    sreg <- cbind(sreg, "")
    colnames(sreg)[4] <- " "
    pval1 <- as.numeric(sreg[, 3])
    sreg[, 4] <- paste(ifelse(pval1<.1, "*", ""), ifelse(pval1<.05, "*", ""), ifelse(pval1<.01, "*", ""), sep = "")
    sreg[pval1 < 1e-4, 3] <- "<1e-4"
    print(sreg, ...)
    cat("\nSummary of dispersion coefficients: \n")
    sdisp <- data.frame(sdisp);colnames(sdisp) <- " "
    rownames(sdisp) <- paste("var", processname, c("clust", "subj"), sep = ".")
    print(sdisp, ...)
    invisible(x)
}
#************ print.summary.unirec 

################# Plot method

#****f* methodsUni/plot.unirec
#  NAME
#    plot.unirec --- S3 plot method for unirec objects
#  FUNCTION
#    User-facing plotting function to plot survivor functions
#  INPUTS
#    see package documentation
#  SYNOPSIS
plot.unirec <- function(x, main = NULL, xscale = 1, hazscale = 1, add = FALSE, ...)
#  SOURCE
#
{
    processname <- attr(x, "processname")
    if(!is.null(main)) processname <- main
    plotsurvivor(x$hazard$breaks * xscale, x$hazard$hazard / xscale * hazscale, main = processname, add = add, ...)
    invisible(x)
}

#************ plot.unirec 
