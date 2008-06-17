########### Method dispatch (univariate)

unirec<-function(x,...)
{
    UseMethod("unirec")
}

unirec.data.frame<-function(x,...)
{
    call<-match.call()
    if(dim(x)[2]>7 && all(colnames(x)[1:7]==c("i","j","k","r","start","stop","delta"))) {
        x<-x[order(x$i,x$j,x$k),]
        varnames<-colnames(x)[-(1:7)]
        x<-cbind(x[,1:7],0,x[,-(1:7)])
        colnames(x)[-(1:7)]<-c("Delta",varnames)
        out<-bivrec.agdata(x,fixzero=c("clust2","subj2","cov"),...)
        out$call<-call
        return(out)
    } else {
        stop("input data frame needs to have colnames i,j,k,r,start,stop,delta")
    }   
}

unirec.formula<-function(formula,data=parent.frame(),...)
{
    # in part based on coxph function
    call<-match.call()
    m<-match.call(expand.dots=FALSE)
    if(is.matrix(eval(m$data,sys.parent()))) m$data<-as.data.frame(data)
    m$...<-NULL
    m[[1]]<-as.name("model.frame")
    special<-c("id","cluster","strata")
    Terms <- if (missing(data)) terms(formula, special) else terms(formula, special, data = data)    
    m$formula<-Terms
    m<-eval(m,sys.parent())
    n<-nrow(m)
    resp <- model.extract(m, "response")
    if (!is.Surv(resp)) stop("model response must be a Surv object")
    if(dim(resp)[2]!=3) stop("response must be of type Surv(start,stop,status)")
    start<-resp[,1]
    stop<-resp[,2]
    delta<-resp[,3]
    dropx<-NULL
    clusterind<-attr(Terms,"specials")$cluster
    # Cluster handling
    clusternames<-NULL
    if(length(clusterind)>0){
        cluster<-m[,clusterind]
        if(is.factor(cluster)) clusternames<-levels(cluster)
        if(is.numeric(cluster)) clusternames<-as.character(unique(cluster))
        i<-as.numeric(as.factor(cluster))
        tempc <- untangle.specials(Terms, "cluster", 1:10)
        ord <- attr(Terms, "order")[tempc$terms]
        if (any(ord > 1)) stop("Cluster cannot be used in an interaction")
        dropx <- c(dropx,tempc$terms)
    }else{
        cluster<-rep(1,n)
        i<-rep(1,n)
    }
    # ID handling
    idind<-attr(Terms,"specials")$id
    if(length(idind)==0) stop("a subject id (unique within clusters) is required")
    else{
        id<-m[,idind]
        subjnames<-paste(cluster,id,sep=".")
        j<-rep(0,n)
        for(ii in 1:max(i)){
            iids<-unique(id[i==ii])
            jj<-1
            for(thisid in iids){
                j[i==ii & id==thisid]<-jj
                jj<-jj+1
            }
        }
        tempi<-untangle.specials(Terms, "id", 1:10)
        ord <- attr(Terms, "order")[tempi$terms]
        if (any(ord > 1)) stop("id cannot be used in an interaction")
        dropx <- c(dropx,tempi$terms)
    }
    # Stratum handling
    stratind<-attr(Terms,"specials")$strata
    if(length(stratind)>0) 
    {
        strat<-m[,stratind]
        ustrat<-unique(strat)
        stratnames<-as.character(ustrat)
        r<-rep(0,n)
        for(rr in 1:length(ustrat)){
            r[strat==ustrat[rr]]<-rr
        }
        tempr<-untangle.specials(Terms, "strata", 1:10)
        ord <- attr(Terms, "order")[tempr$terms]
        if (any(ord > 1)) stop("strata cannot be used in an interaction")
        dropx <- c(dropx,tempr$terms)
    }else{
        r<-rep(1,n)
        stratnames<-"1"
    }
    
    Ki<-table(i*1e6+j)
    k<-unlist(as.vector(sapply(Ki,function(x) 1:x)))
    newTerms <- if(length(dropx))  Terms[-dropx] else Terms
    X <- model.matrix(newTerms, m)
    X<-X[,-1,drop=FALSE]
    agdata<-as.data.frame(cbind(i,j,k,r,start,stop,delta,X))
    sortord<-order(agdata$i,agdata$j)
    agdata<-agdata[sortord,]
    subjnames<-unique(subjnames[sortord])
    processname<-as.character(call$formula[[2]][[4]])
    fit<-unirec(agdata,...)
    fit$call<-call
    fit<-cleanunirecoutput(fit,clusternames,subjnames,stratnames,processname)
    return(fit)
}

################# Output cleanup


cleanunirecoutput<-function(fit,clusternames,subjnames,stratnames,processname)
{
    out<-NULL
    out$call<-fit$call
    # regression
    out$regression=list(coefficients=fit$regressionoutput$betahat,
                        loglik=fit$regressionoutput$loglik1)
    # frailties
    clust<-fit$frailtyoutput$Uihat
    subj<-fit$frailtyoutput$Uijhat
    names(clust)<-clusternames
    names(subj)<-subjnames
    out$frailty=list(clust=clust,subj=subj)
    # dispersion
    disp<-list(clust=fit$dispparams$sigma2hat, subj=fit$dispparams$nu2hat)
    out$dispersion<-disp
    # baseline
    breaks<-fit$regressionoutput$as
    basehaz<-fit$regressionoutput$alphars
    rownames(breaks)<-rownames(basehaz)<-stratnames
    out$hazard=list(breaks=breaks,hazard=basehaz)
    # summaries
    summary.reg<-fit$summary.reg[,1:3]
    colnames(summary.reg)<-c("coef","sd","pval")
    summary.disp<-fit$summary.disp[c(1,3)]
    names(summary.disp)<-c("clust","subj")
    out$summaries<-list(regression=summary.reg,dispersion=summary.disp)
    attr(out,"processname")<-processname
    class(out)<-"unirec"
    return(out)
}

################# Print and summary
print.unirec<-function(x,...)
{
    processname<-attr(x,"processname")
    cat("Process ", processname,":\n",sep="")
    print(x$regression$coefficients,...)
    invisible(x)
}

summary.unirec<-function(object,digits=4,...)
{
    x<-object
    out<-NULL
    dots<-as.list(substitute(list(...)))[-1]
    sreg<-format(x$summaries$regression,digits=digits,...)
    sdisp<-format(x$summaries$dispersion,digits=digits,...)
    out<-list(call=x$call,summary.reg=sreg,summary.disp=sdisp)
    processname<-attr(x,"processname")
    attr(out,"processname")<-processname
    class(out)<-"summary.unirec"   
    return(out)
}

print.summary.unirec<-function(x,...)
{
    processname<-attr(x,"processname")
    dots<-x$dots
    sreg<-x$summary.reg
    sdisp<-x$summary.disp
    colnames(sreg)[1:3]<-paste(processname,c("coef","sd","pval"),sep=".")
    cat("\nSummary of regression coefficients: \n")
    sreg<-cbind(sreg,"")
    colnames(sreg)[4]<-" "
    pval1<-as.numeric(sreg[,3])
    sreg[,4]<-paste(ifelse(pval1<.1,"*",""),ifelse(pval1<.05,"*",""),ifelse(pval1<.01,"*",""),sep="")
    sreg[pval1<1e-4,3]<-"<1e-4"
    print(sreg,...)
    cat("\nSummary of dispersion coefficients: \n")
    sdisp<-data.frame(sdisp);colnames(sdisp)<-" "
    rownames(sdisp)<-paste("var",processname,c("clust","subj"),sep=".")
    print(sdisp,...)
    invisible(x)
}

################# Plot method

plot.unirec<-function(x,main=NULL,xscale=1,hazscale=1,add=FALSE,...)
{
    processname<-attr(x,"processname")
    if(!is.null(main)) processname<-main
    plotsurvivor(x$hazard$breaks*xscale,x$hazard$hazard/xscale*hazscale,main=processname,add=add,...)
    invisible(x)
}

