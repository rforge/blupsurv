################# Method dispatch (bivariate)

bivrec <- function(x, ...)
{
    UseMethod("bivrec")
}

bivrec.data.frame <- function(x, ...)
{
    call <- match.call()
    if(dim(x)[2] > 8 && all(colnames(x)[1:8] == c("i", "j", "k", "r", "start", "stop", "delta", "Delta"))) {
        x <- x[order(x$i, x$j, x$k), ]
        out <- bivrec.agdata(x, ...)
        out$call <- call
        return(out)
    } else {
        if(dim(x)[2] == 8) stop("at least one covariate is required")
        else stop("input data frame needs to have colnames i, j, k, r, start, stop, delta, Delta")
    }
}

bivrec.formula <- function(formula, data = parent.frame(), ...)
{
    # in part based on coxph function
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, sys.parent()))) m$data <- as.data.frame(data)
    m$...<-NULL
    m[[1]] <- as.name("model.frame")
    special <- c("id", "cluster", "strata")
    Terms <- if (missing(data)) terms(formula, special) else terms(formula, special, data = data)    
    m$formula <- Terms
    oldNAoption <- getOption("na.action"); options(na.action = na.fail)
    m <- eval(m, sys.parent())
    options(na.action = oldNAoption)
    n <- nrow(m)
    resp <- model.extract(m, "response")
    if (!is.Surv2(resp)) stop("model response must be a Surv2 object")
    start <- resp[, "start"]
    stop <- resp[, "stop"]
    delta <- resp[, "status1"]
    Delta <- resp[, "status2"]
    processnames <- attr(resp, "processnames")
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
        if(is.factor(id)) levels(id) <- 1:length(levels(id))
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
    
    Ki <- table(i * 1e6 + j)
    k <- unlist(as.vector(sapply(Ki, function(x) 1:x)))
    newTerms <- if(length(dropx))  Terms[ - dropx] else Terms
    X <- model.matrix(newTerms, m)
    X <- X[, -1, drop = FALSE]
    agdata <- as.data.frame(cbind(i, j, k, r, start, stop, delta, Delta, X))
    sortord <- order(agdata$i, agdata$j)
    agdata <- agdata[sortord, ]
    subjnames <- unique(subjnames[sortord])
    check <- checkinput.bivrec(agdata, clusternames, subjnames, stratnames, processnames)
    fit <- bivrec(agdata, ...)
    fit$call <- call
    fit <- cleanbivrecoutput(fit, clusternames, subjnames, stratnames, processnames)
    return(fit)
}


################# Bivariate survival class

Surv2 <- function(start, stop, status1, status2)
{
    l <- length(start)
    if(length(stop) != l || length(status1) != l || length(status2) != l)
        stop("All four vectors must be of the same length")
    if(!all(sort(unique(status1)) == c(0, 1)) || !all(sort(unique(status2)) == c(0, 1)))
        stop("status must be 0 - 1")
    if(!all(stop >= start))
        stop("stop must be greater than start")
    out <- cbind(start, stop, status1, status2)
    attr(out, "processnames") <- as.character(c(match.call()[[4]], match.call()[[5]]))
    class(out) <- "Surv2"
    invisible(out)
}

is.Surv2 <- function(x) inherits(x, "Surv2")

id <- function(x) x

################# Output cleanup

cleanbivrecoutput <- function(fit, clusternames, subjnames, stratnames, processnames)
{
    out <- NULL
    out$call <- fit$call
    # regression
    out$regression = list(coefficients1 = fit$regressionoutput$betahat,
                        coefficients2 = fit$regressionoutput$betadhat,
                        loglik1 = fit$regressionoutput$loglik1,
                        loglik2 = fit$regressionoutput$loglik2)
    # frailties
    clust1 <- fit$frailtyoutput$Uihat
    clust2 <- fit$frailtyoutput$Vihat
    subj1 <- fit$frailtyoutput$Uijhat
    subj2 <- fit$frailtyoutput$Vijhat
    names(clust1) <- names(clust2) <- clusternames
    names(subj1) <- names(subj2) <- subjnames
    out$frailty = list(clust1 = clust1, clust2 = clust2, subj1 = subj1, subj2 = subj2)
    # dispersion
    disp <- fit$dispparams
    names(disp) <- c("clust1", "clust2", "subj1", "subj2", "cov")
    out$dispersion <- disp
    # baseline
    breaks1 <- fit$regressionoutput$as
    breaks2 <- fit$regressionoutput$asd
    basehaz1 <- fit$regressionoutput$alphars
    basehaz2 <- fit$regressionoutput$alpharsd
    rownames(breaks1) <- rownames(breaks2) <- rownames(basehaz1) <- rownames(basehaz2) <- stratnames
    out$hazard = list(breaks1 = breaks1, breaks2 = breaks2, hazard1 = basehaz1, hazard2 = basehaz2)
    # summaries
    summary.reg <- fit$summary.reg
    colnames(summary.reg) <- c("coef1", "sd1", "pval1", "coef2", "sd2", "pval2")
    summary.disp <- fit$summary.disp
    names(summary.disp) <- c("clust1", "clust2", "subj1", "subj2", "cov", "corr")
    out$summaries <- list(regression = summary.reg, dispersion = summary.disp)
    attr(out, "processnames") <- processnames
    class(out) <- "bivrec"
    return(out)
}

################# Print and summary

print.bivrec <- function(x, ...)
{
    processnames <- attr(x, "processnames")
    cat("Process 1 (", processnames[1], "):\n", sep = "")
    print(x$regression$coefficients1, ...)
    cat("Process 2 (", processnames[2], "):\n", sep = "")
    print(x$regression$coefficients2, ...)
    invisible(x)
}

summary.bivrec <- function(object, digits = 4, ...)
{
    x <- object
    out <- NULL
    dots <- as.list(substitute(list(...)))[ - 1]
    sreg <- format(x$summaries$regression, digits = digits, ...)
    sdisp <- format(x$summaries$dispersion, digits = digits, ...)
    out <- list(call = x$call, summary.reg = sreg, summary.disp = sdisp)
    processnames <- attr(x, "processnames")
    attr(out, "processnames") <- processnames
    class(out) <- "summary.bivrec"   
    return(out)
}

print.summary.bivrec <- function(x, ...)
{
    processnames <- attr(x, "processnames")
    dots <- x$dots
    sreg <- x$summary.reg
    sdisp <- x$summary.disp
    colnames(sreg)[1:3] <- paste(processnames[1], c("coef", "sd", "pval"), sep = ".")
    colnames(sreg)[4:6] <- paste(processnames[2], c("coef", "sd", "pval"), sep = ".")
    cat("\nSummary of regression coefficients: \n")
    sreg <- cbind(sreg[1:3], "", sreg[4:6], "")
    colnames(sreg)[c(4, 8)] <- " "
    pval1 <- as.numeric(sreg[, 3])
    sreg[, 4] <- paste(ifelse(pval1<.1, "*", ""), ifelse(pval1<.05, "*", ""), ifelse(pval1<.01, "*", ""), sep = "")
    sreg[!is.na(pval1) & pval1 < 1e-4, 3] <- "<1e-4"
    pval2 <- as.numeric(sreg[, 7])
    sreg[, 8] <- paste(ifelse(pval2<.1, "*", ""), ifelse(pval2<.05, "*", ""), ifelse(pval2<.01, "*", ""), sep = "")
    sreg[!is.na(pval2) & pval2 < 1e-4, 7] <- "<1e-4"
    for(ind in 1:dim(sreg)[2]) {
        fixrow <- grep("NA", sreg[, ind]); 
        if(length(fixrow) > 0) sreg[fixrow, ind] <- ""
    }
    print(sreg, ...)
    cat("\nSummary of dispersion coefficients: \n")
    sdisp <- data.frame(sdisp);colnames(sdisp) <- " "
    rownames(sdisp)[c(1, 3)] <- paste("var", processnames[1], c("clust", "subj"), sep = ".")
    rownames(sdisp)[c(2, 4)] <- paste("var", processnames[2], c("clust", "subj"), sep = ".")
    rownames(sdisp)[c(5, 6)] <- c("covariance", "correlation")
    print(sdisp, ...)
    invisible(x)
}

################# Plot method

plot.bivrec <- function(x, which = c(0, 1, 2), main = NULL, xscale = 1, hazscale = 1, add = FALSE, legend = NULL, ...)
{
    if(length(which) > 1) which <- 0
    which <- c(0, 1, 2)[pmatch(which, c(0, 1, 2))]
    processnames <- attr(x, "processnames")
    if(!is.null(main)) processnames <- c(main[1], ifelse(length(main) > 1, main[2], main[1]))
    if(which == 0) par(mfrow = c(1, 2))
    if(which == 0 | which == 1) plotsurvivor(x$hazard$breaks1 * xscale, x$hazard$hazard1 / xscale * hazscale, main = processnames[1], add = add, legend = legend, ...)
    if(which == 0 | which == 2) plotsurvivor(x$hazard$breaks2 * xscale, x$hazard$hazard2 / xscale * hazscale, main = processnames[2], add = add, legend = legend, ...)
    invisible(x)
}

################# Plotter

plotsurvivor <- function(breaks, hazard, col = 1:10, lty = rep(1, 10), main = "survivor", xlab = "Time",
    ylab = "Survivor function", xlim = NULL, ylim = NULL, type = "s", add = FALSE, legend = NULL, ...)
{
    interlength <- t(apply(breaks, 1, diff))
    survival <- cbind(1, exp(-t(apply(hazard[, -dim(hazard)[2], drop = FALSE] * interlength, 1, cumsum))))
    if(is.null(xlim))   xlim = c(0, max(breaks))
    if(is.null(ylim))    ylim = c(max(0, min(survival-.05)), 1)
    if(!add)  plot(-1e5, -1e5, xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab = ylab, ...)
    for(r in 1:dim(breaks)[1]){
        lines(breaks[r, ], survival[r, ], col = col[r], lty = lty[r], type = type, ...)
    }
    if(!is.null(legend)) legend("bottomleft", rownames(hazard), col = col, lty = lty, ...)
}

################# Input check

checkinput.bivrec <- function(agdata, clusternames, subjnames, stratnames, processnames)
{
    # Check that start - stop times match
    starttimes <- agdata$start[agdata$start > 0]
    stoptimes <- agdata$stop[which(agdata$start > 0) - 1]
    if(!isTRUE(all.equal(starttimes, stoptimes)))
        stop("Interval start times must be equal to the preceding stop times")
    
    # Check that there are events for each cluster
    nclustevents1 <- table(agdata$i, agdata$delta)
    nclustevents2 <- table(agdata$i, agdata$Delta)
    if(!all(nclustevents1[, 2] > 0)){
        badclust <- which(nclustevents1[, 2] <= 0)[1]
        stop(paste("Cluster", clusternames[badclust], "has no events of type", processnames[1]))
    }
    if(!all(nclustevents2[, 2] > 0)){
        badclust <- which(nclustevents2[, 2] <= 0)[1]
        stop(paste("Cluster", clusternames[badclust], "has no events of type", processnames[2]))
    }
    # Check for missing values or NaN
    if(any(sapply(agdata, is.infinite)))
        stop("Infinite values in data")
}

