! Compile string: R CMD SHLIB bivrec.f

!****h* /fortran
! NAME
!   fortran --- FORTRAN routines for fast computation
! FUNCTION
!   FORTRAN routines called by the bivrec.agdata fitter. All of these
!   are also implemented in R in ZZdebug, but with less BLAS magic.
!   All input matrices are in the form of column-major vectors. This
!   will not be mentioned elsewere.
!****** 

!****m* fortran/fmkmrs
! NAME
!   fmkmrs --- terms for profile likelihood
! FUNCTION
!     Compute terms needed in the profile likelihood and gradient.
!     See also the R function makemrs.
! INPUTS
!    betahat    regression parameters (length p)
!    index      d x 6 matrix of indices i,j,k,r,smin,smax 
!    delta      length d vector of event indicators (double)
!    times      length d vector of actual time span for each row (double)
!    Z          d x ncovs matrix of covariates (as a double row vector)
!    as         nr x ns matrix of discretization breakpoints (double) 
!    Uijmat     m x maxj matrix of frailty estimates (double row vector)
!    d          number of rows of Z
!    ncovs      number of covariates
!    nr         number of strata
!    ns         number of breakpoints
!    m          number of clusters
!    maxj       largest cluster size
!    mrs        nr x ns matrix for output
!    mkgr       integer indicator for whether to make gradient
! OUTPUTS
!    mrsgr      nr x ns x ncovs gradient for output
! SYNOPSIS
      subroutine fmkmrs(betahat,index,times,Z,as,Uijmat,d,ncovs,nr,ns,
     $  m,maxj,mrs,mrsgr,mkgr)
! SOURCE
      integer d,ncovs,nr,ns,nk,m,maxj,mkgr
      integer index(d,6)
      double precision betahat(ncovs)
      double precision times(d)
      double precision Z(d,ncovs)
      double precision Uijmat(m,maxj)
      double precision lp,time,pred1,pred2,thisA
      double precision mrs(nr,ns)
      double precision mrsgr(nr,ns,ncovs)
      double precision as(nr,ns)
      
      integer ind,i,j,k,r,s,smax,smin,maxs    
      parameter(icol=1, jcol=2, kcol=3, ircol=4, ismincol=5, ismaxcol=6)
      
!     clear output arrays
      call dzero(mrs,nr*ns)
      call dzero(mrsgr,nr*ns*ncovs)

      maxs=(ns-1)*95/100      
      do 100 ind=1,d
        i=index(ind,icol)
        j=index(ind,jcol)
        k=index(ind,kcol)
        r=index(ind,ircol)
        smin=index(ind,ismincol)
        smax=index(ind,ismaxcol)
        time=times(ind)
        lp=ddot(ncovs,Z(ind,1),d,betahat,1)
!       ignore the last 5% of intervals
!        if(smax.GT.maxs) then
!            smax=maxs
!        endif
        pred1=Uijmat(i,j)*exp(lp)  
        do 110 s=smin,smax
            thisA=0.d0
            call A(thisA,time,as,r,s,nr,ns)
            pred2=pred1*thisA
            mrs(r,s)=mrs(r,s)+pred2
            if(mkgr.EQ.1) then
                call daxpy(ncovs,pred2,Z(ind,1),d,mrsgr(r,s,1),nr*ns)
            endif
 110    continue
                
 100  continue
      end
!******
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!****m* fortran/A
! NAME
!   A --- time in an interval
! FUNCTION
!    Computes the amount of time in interval (r,s) prior to time t
! INPUTS
!    out    output storage
!    time   time of interest
!    as     nr x ns matrix of breakpoints
!    r      stratum
!    s      interval
!    nr         number of strata
!    ns         number of breakpoints
! SYNOPSIS
      subroutine A(out,time,as,r,s,nr,ns)
! SOURCE
      integer nr,ns,r,s
      double precision out,time,as(nr,ns)
      
      s=s+1
      if(s.EQ.0) then
        out=0.d0
      else if(time.LT.as(r,s-1)) then
        out=0.d0
      else if(time.GE.as(r,s)) then
        out=as(r,s)-as(r,s-1)
      else
        out=time-as(r,s-1)
      endif
      s=s-1
      end
!******

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!****m* fortran/fproflik
! NAME
!   fproflik --- compute profile likelihood
! FUNCTION
!    Compute profile likelihood of regression parameters for a single process.
!    See also the R function mkproflik.
! INPUTS
!    betahat    regression parameters (length p)
!    index      d x 6 matrix of indices i,j,k,r,smin,smax 
!    delta      length d vector of event indicators (double)
!    times      length d vector of actual time span for each row (double)
!    Z          d x ncovs matrix of covariates (as a double row vector)
!    as         nr x ns matrix of discretization breakpoints (double) 
!    Uijmat     m x maxj matrix of frailty estimates (double row vector)
!    d          number of rows of Z
!    ncovs      number of covariates
!    nr         number of strata
!    ns         number of breakpoints
!    m          number of clusters
!    maxj       largest cluster size
! OUTPUTS
!    lik        storage for output
! SYNOPSIS
      subroutine fproflik(betahat,index,delta,times,Z,as,Uijmat,d,
     $ ncovs,nr,ns,m,maxj,lik)
! SOURCE
      integer d,ncovs,nr,ns,nk,m,maxj
      integer index(d,6)
      double precision betahat(ncovs)
      double precision delta(d),times(d)
      double precision Z(d,ncovs)
      double precision Uijmat(m,maxj)
      double precision lp,deltat,lik
      double precision mrs(nr,ns)
      double precision mrsgr(nr,ns,ncovs)
      double precision as(nr,ns)
      integer ind,i,j,k,r,s,smax,smin

      parameter(icol=1, jcol=2, kcol=3, ircol=4, ismincol=5, ismaxcol=6)
      
      call fmkmrs(betahat,index,times,Z,as,Uijmat,d,ncovs,nr,ns,
     $  m,maxj,mrs,mrsgr,0)
      
      lik=0.d0
      do 100 ind=1,d
        i=index(ind,icol)
        j=index(ind,jcol)
        k=index(ind,kcol)
        r=index(ind,ircol)
        smin=index(ind,ismincol)
        smax=index(ind,ismaxcol)
        deltat=delta(ind)
        lp=ddot(ncovs,Z(ind,1),d,betahat,1)
        lik=lik+deltat*(log(Uijmat(i,j))-log(mrs(r,smax))+lp)
 100  continue
      end
!******

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!****m* fortran/fprofgr
! NAME
!  fprofgr --- compute profile likelihood gradient
! FUNCTION
!    Compute profile likelihood of regression parameters for a single process.
!    See also the R function mkproflik.
! INPUTS
!    betahat    regression parameters (length p)
!    index      d x 6 matrix of indices i,j,k,r,smin,smax 
!    delta      length d vector of event indicators (double)
!    times      length d vector of actual time span for each row (double)
!    Z          d x ncovs matrix of covariates (as a double row vector)
!    as         nr x ns matrix of discretization breakpoints (double) 
!    Uijmat     m x maxj matrix of frailty estimates (double row vector)
!    d          number of rows of Z
!    ncovs      number of covariates
!    nr         number of strata
!    ns         number of breakpoints
!    m          number of clusters
!    maxj       largest cluster size
! OUTPUTS
!    gr         ncovs vector for output
! SYNOPSIS
      subroutine fprofgr(betahat,index,delta,times,Z,as,Uijmat,d,
     $ ncovs,nr,ns,m,maxj,gr)
! SOURCE
      integer d,ncovs,nr,ns,nk,m,maxj
      integer index(d,6)
      double precision betahat(ncovs)
      double precision delta(d),times(d)
      double precision Z(d,ncovs)
      double precision Uijmat(m,maxj)
      double precision lp,deltat,gr(ncovs)
      double precision mrs(nr,ns)
      double precision mrsgr(nr,ns,ncovs)
      double precision as(nr,ns)
      integer ind,i,j,k,r,s,smax,smin

      parameter(icol=1, jcol=2, kcol=3, ircol=4, ismincol=5, ismaxcol=6)
      
      call fmkmrs(betahat,index,times,Z,as,Uijmat,d,ncovs,nr,ns,
     $  m,maxj,mrs,mrsgr,1)
      
      lik=0.d0
      do 100 ind=1,d
        i=index(ind,icol)
        j=index(ind,jcol)
        k=index(ind,kcol)
        r=index(ind,ircol)
        smin=index(ind,ismincol)
        smax=index(ind,ismaxcol)
        deltat=delta(ind)
        if(deltat.EQ.1.d0) then
            call daxpy(ncovs,1.d0,Z(ind,1),d,gr,1)
            call daxpy(ncovs,-1.d0/mrs(r,smax),mrsgr(r,smax,1),
     $                 nr*ns,gr,1)
        endif
 100  continue
      end
!******

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!****m* fortran/fsmuij
! NAME
!  fsmuij --- compute sums of events and means
! FUNCTION
!    Compute the sums of mu_rijks and delta_rijks for every (i,j). 
!    See also the R function Smuij.
! INPUTS
!    betahat    regression parameters (length p)
!    index      d x 6 matrix of indices i,j,k,r,smin,smax 
!    delta      length d vector of event indicators (double)
!    Z          d x ncovs matrix of covariates (double row vector)
!    alphars    nr x ns matrix of baseline hazards (double) 
!    d          number of rows of Z
!    ncovs      number of covariates
!    nr         number of strata
!    ns         number of breakpoints
!    m          number of clusters
!    maxj       largest cluster size
! OUTPUTS
!    Smu        m x maxj matrix containg sum(mu_rijks) for all (i,j)
!    Sdelta     m x maxj matrix containg sum(delta_rijks) for all (i,j)
! SYNOPSIS
      subroutine fsmuij(betahat,index,delta,Z,alphars,as,
     $ d,ncovs,nr,ns,m,maxj,Smu,Sdelta)
! SOURCE
      integer d,ncovs,nr,ns,nk,m,maxj
      integer index(d,6)
      double precision betahat(ncovs)
      double precision delta(d)
      double precision Z(d,ncovs)
      double precision alphars(nr,ns)
      double precision as(nr,ns)
      double precision lp, deltat,w,mu,power
      double precision Smu(m,maxj),Sdelta(m,maxj)
      integer ind,i,j,k,r,s,smax,smin,maxs
      
      parameter(icol=1, jcol=2, kcol=3, ircol=4, ismincol=5, ismaxcol=6)

      maxs=(ns-1)*95/100
      do 100 ind=1,d
!	call intpr("ind",-1,ind,1)	 
       i=index(ind,icol)
        j=index(ind,jcol)
        k=index(ind,kcol)
        r=index(ind,ircol)
        smin=index(ind,ismincol)
        smax=index(ind,ismaxcol)
        deltat=delta(ind)
        lp=ddot(ncovs,Z(ind,1),d,betahat,1)
!       ignore the last 5% of intervals
!        if(smax.GT.maxs) then
!            smax=maxs
!        endif  
        Sdelta(i,j)=Sdelta(i,j)+deltat
        do 110 s=smin,smax
           mu=exp(lp)*alphars(r,s)*(as(r,s+1)-as(r,s))
           mu=(1.d0-exp(-mu))
           Smu(i,j)=Smu(i,j)+mu
 110    continue   
 100  continue
      end
!******
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!****m* fortran/fsmud2
! NAME
!  fsmuij --- compute squared means and cross-terms
! FUNCTION
!    Compute the sums of mu_rijks * delta_rijks and
!    mu_rijks^2 for every (i,j). 
!    See also the R function Smud2.
! INPUTS
!    betahat    regression parameters (length p)
!    index      d x 6 matrix of indices i,j,k,r,smin,smax 
!    delta      length d vector of event indicators (double)
!    Z          d x ncovs matrix of covariates (double row vector)
!    alphars    nr x ns matrix of baseline hazards (double) 
!    d          number of rows of Z
!    ncovs      number of covariates
!    nr         number of strata
!    ns         number of breakpoints
!    m          number of clusters
!    maxj       largest cluster size
! OUTPUTS
!    Smud       m x maxj matrix containg sum(mu*delta) for all (i,j)
!    Smu2       m x maxj matrix containg sum(mu^2) for all (i,j)
! SYNOPSIS
      subroutine fsmud2(betahat,index,delta,Z,alphars,as,
     $ d,ncovs,nr,ns,Smud,Smu2)
! SOURCE
      integer d,ncovs,nr,ns,nk
      integer index(d,6)
      double precision betahat(ncovs)
      double precision delta(d)
      double precision Z(d,ncovs)
      double precision alphars(nr,ns)
      double precision as(nr,ns)
      double precision lp, deltat,deltatt,w,mu,power
      double precision Smud,Smu2
      integer ind,i,j,k,r,s,smax,smin,maxs
      
      parameter(icol=1, jcol=2, kcol=3, ircol=4, ismincol=5, ismaxcol=6)

      maxs=(ns-1)*95/100
      do 100 ind=1,d
!	call intpr("ind",-1,ind,1)	 
       i=index(ind,icol)
        j=index(ind,jcol)
        k=index(ind,kcol)
        r=index(ind,ircol)
        smin=index(ind,ismincol)
        smax=index(ind,ismaxcol)
        deltat=delta(ind)
        lp=ddot(ncovs,Z(ind,1),d,betahat,1)
!       ignore the last 5% of intervals
!        if(smax.GT.maxs) then
!            smax=maxs
!        endif  
        do 110 s=smin,smax
           mu=exp(lp)*alphars(r,s)*(as(r,s+1)-as(r,s))
           mu=(1.d0-exp(-mu))
           if(s.EQ.smax) then
                deltatt=deltat
           else
                deltatt=0.d0
           endif
           Smud=Smud+(mu-deltatt)**2
           Smu2=Smu2+mu**2
 110    continue   
 100  continue
      end
!******
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
!****m* fortran/fmkalpha2
! NAME
!  fsmuij --- compute baseline hazard
! FUNCTION
!   Compute the MLEs for the baseline hazard parameters alphars, 
!   given estimates of regression parameters andf frailties, for
!   a single recurrent event process.
!   See also the R function makealphars2.
! INPUTS
!    betahat    regression parameters (length p)
!    index      d x 6 matrix of indices i,j,k,r,smin,smax 
!    delta      length d vector of event indicators (double)
!    times      length d vector of event times
!    Z          d x ncovs matrix of covariates (double row vector)
!    alphars    nr x ns matrix of baseline hazards (double) 
!    as         nr x ns matrix of discretization breakpoints (double) 
!    Uijmat     m x maxj matrix of frailty estimates (double row vector)
!    d          number of rows of Z
!    ncovs      number of covariates
!    nr         number of strata
!    ns         number of breakpoints
!    m          number of clusters
!    maxj       largest cluster size
! OUTPUTS
!    alphars updated with baseline hazard MLEs
! SYNOPSIS
      subroutine fmkalpha2(betahat,index,delta,times,Z,alphars,as,
     $ Uijmat,d,ncovs,nr,ns,m,maxj)
! SOURCE
      integer d,ncovs,nr,ns,m,maxj
      integer index(d,6)
      double precision betahat(ncovs)
      double precision delta(d),times(d)
      double precision Z(d,ncovs)
      double precision alphars(nr,ns)
      double precision as(nr,ns)
      double precision Uijmat(m,maxj)
      double precision pred,thisA,time
      double precision Srs(nr,ns), drs(nr,ns)
      integer ind,i,j,k,r,s,smax,smin
      
      parameter(icol=1, jcol=2, kcol=3, ircol=4, ismincol=5, ismaxcol=6)
      
      call dzero(Srs,nr*ns)
      call dzero(drs,nr*ns)

      do 100 ind=1,d
        i=index(ind,icol)
        j=index(ind,jcol)
        k=index(ind,kcol)
        r=index(ind,ircol)
        smin=index(ind,ismincol)
        smax=index(ind,ismaxcol)
        deltat=delta(ind)
        time=times(ind)
        pred=exp(ddot(ncovs,Z(ind,1),d,betahat,1))*Uijmat(i,j)
        drs(r,smax)=drs(r,smax)+deltat
        do 110 s=smin,smax
            thisA=0.d0
            call A(thisA,time,as,r,s,nr,ns)
            Srs(r,s)=Srs(r,s)+pred*thisA
 110    continue
 100  continue

      do 200 r=1,nr
        do 210 s=1,ns
            if(Srs(r,s).EQ.0.d0) then
                alphars(r,s)=100.d0
            else
                alphars(r,s)=drs(r,s)/Srs(r,s)
            endif
 210    continue
 200  continue
      end
!******

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!****m* fortran/fmkfrail
! NAME
!  fmkfrail --- compute BLUP frailty estimates
! FUNCTION
!    Compute cluster- and subject-level frailty BLUPs under the
!    auxiliary Poisson model. 
!    See also the R function fupdatefrailties3.
! INPUTS
!   Too many to list, but analogous to those in other routines in
!   this file. Inputs appended with d refer to process 2.
! OUTPUTS
!    Uihat      vector of cluster-level frailties for process 1
!    Vihat      vector of cluster-level frailties for process 2
!    Uijhat     vector of subject-level frailties for process 1
!    Vijhat     vector of subject-level frailties for process 2
!    p,q,r,s, etc   intermediate values
! SYNOPSIS
      subroutine fmkfrail(index,indexd,delta,deltad,Z,Zd,alphars,
     $ alpharsd,as,asd,betahat,betadhat,m,Ji,jimax,Jicum,jisum,d1,
     $ d2,ncovs1,ncovs2,nr,ns,nsd,sigma2,sigma2d,nu2,nu2d,theta,Uihat,
     $ Vihat,Uijhat,Vijhat,pi,qi,ri,si,piprime,qiprime,riprime,siprime,
     $ pij,qij,rij,sij,pijprime,qijprime,rijprime,sijprime,wi,wij,zij)
! SOURCE
      integer d1,d2,ncovs1,ncovs2,nr,ns,nsd,m,jimax
      integer Ji(m),Jicum(m)
      integer index(d1,6), indexd(d2,6)
      double precision delta(d1), deltad(d2)
      double precision Z(d1,ncovs1), Zd(d2,ncovs2)
      double precision alphars(nr,ns),alpharsd(nr,nsd)
      double precision as(nr,ns), asd(nr,nsd)
      double precision betahat(ncovs1), betadhat(ncovs2)
      double precision sigma2,sigma2d,nu2,nu2d,theta
      double precision Uihat(m),Vihat(m),Uijhat(jisum),Vijhat(jisum)
      double precision pi(m),qi(m),ri(m),si(m),piprime(m),qiprime(m),
     $ riprime(m),siprime(m),pij(m,jimax),qij(m,jimax),rij(m,jimax),
     $ sij(m,jimax),pijprime(m,jimax),qijprime(m,jimax),
     $ rijprime(m,jimax),sijprime(m,jimax)
      double precision wi(m),wij(m,jimax),zij(m,jimax)
      
      integer i,j
      double precision Smu(m,jimax),Seta(m,jimax),Sdelta(m,jimax),
     $ SDeltad(m,jimax)
      double precision Smuij,Setaij,Sdeltaij,Sdeltaijd
      
      call dzero(Smu,m*jimax)
      call dzero(Seta,m*jimax)
      call dzero(Sdelta,m*jimax)
      call dzero(SDeltad,m*jimax)
      
!     call dblepr("Smu",-1,Smu,m*jimax)

      call fsmuij(betahat,index,delta,Z,alphars,as,d1,ncovs1,nr,ns,
     $ m,jimax,Smu,Sdelta)
      call fsmuij(betadhat,indexd,deltad,Zd,alpharsd,asd,d2,ncovs2,nr,
     $ nsd,m,jimax,Seta,SDeltad)    
      
      do 100 i=1,m
!    compute intermediate values
        do 110 j=1,Ji(i)
            Smuij=Smu(i,j)
            Setaij=Seta(i,j)
            Sdeltaij=Sdelta(i,j)
            SDeltaijd=SDeltad(i,j)
            
            wij(i,j)=nu2-theta*theta*Setaij/(1.0d0+nu2d*Setaij)
            zij(i,j)=nu2d-theta*theta*Smuij/(1.0d0+nu2*Smuij)
            pij(i,j)=Smuij/(1.0d0+wij(i,j)*Smuij)
            qij(i,j)=-theta*Smuij*Setaij/((1.0d0+nu2*Smuij)*
     $          (1.0d0+nu2d*Setaij)-theta*theta*Smuij*Setaij)
            rij(i,j)=qij(i,j)
            sij(i,j)=Setaij/(1.0d0+zij(i,j)*Setaij)
            pijprime(i,j)=Sdeltaij/(1.0d0+wij(i,j)*Smuij)
            qijprime(i,j)=-theta*Smuij*SDeltaijd/((1.0d0+nu2*Smuij)*
     $          (1.0d0+nu2d*Setaij)-theta*theta*Smuij*Setaij)
            rijprime(i,j)=-theta*Sdeltaij*Setaij/((1.0d0+nu2*Smuij)*
     $          (1.0d0+nu2d*Setaij)-theta*theta*Smuij*Setaij)   
            sijprime(i,j)=SDeltaijd/(1.0d0+zij(i,j)*Setaij)       
            
            piprime(i)=piprime(i)+pijprime(i,j)
            qiprime(i)=qiprime(i)+qijprime(i,j)
            riprime(i)=riprime(i)+rijprime(i,j)
            siprime(i)=siprime(i)+sijprime(i,j)
            pi(i)=pi(i)+pij(i,j)
            qi(i)=qi(i)+qij(i,j)
            ri(i)=ri(i)+rij(i,j)
            si(i)=si(i)+sij(i,j)            
 110    continue     
!  compute cluster-level frailties
        wi(i)=1.0d0/((1.0d0+sigma2*pi(i))*(1.0d0+sigma2d*si(i))
     $      -sigma2*sigma2d*ri(i)*ri(i))
        Uihat(i)=1.0d0+sigma2*wi(i)*((1.0d0+sigma2d*si(i))*
     $      (piprime(i)-pi(i)+qiprime(i)-qi(i))-sigma2d*ri(i)*
     $      (riprime(i)-ri(i)+siprime(i)-si(i)))
        Vihat(i)=1.0d0+sigma2d*wi(i)*((1.0d0+sigma2*pi(i))*
     $      (riprime(i)-ri(i)+siprime(i)-si(i))-sigma2*ri(i)*
     $      (piprime(i)-pi(i)+qiprime(i)-qi(i))) 
        
        if(Uihat(i).LE.0.01) then 
            Uihat(i)=0.01
        endif
        if(Vihat(i).LE.0.01) then 
            Vihat(i)=0.01
        endif
!  compute subject-level frailties
        do 120 j=1,Ji(i)
            Uijhat(Jicum(i)+j)=Uihat(i)-(nu2*pij(i,j)+theta*rij(i,j))
     $           *Uihat(i)-(nu2*qij(i,j)+theta*sij(i,j))*Vihat(i)
     $          +nu2*(pijprime(i,j)+qijprime(i,j))
     $          +theta*(rijprime(i,j)+sijprime(i,j))
            Vijhat(Jicum(i)+j)=Vihat(i)-(theta*qij(i,j)+nu2d*sij(i,j))
     $          *Vihat(i)-(theta*pij(i,j)+nu2d*rij(i,j))*Uihat(i)
     $          +theta*(pijprime(i,j)+qijprime(i,j))
     $          +nu2d*(rijprime(i,j)+sijprime(i,j))           
            if(Uijhat(Jicum(i)+j).LE.0.01) then 
                Uijhat(Jicum(i)+j)=0.01
            endif
            if(Vijhat(Jicum(i)+j).LE.0.01) then 
                Vijhat(Jicum(i)+j)=0.01
            endif
 120    continue
 100  continue   
       
      end
!******
 
!cccccccccccccccccccccccccccccccccccccccc     

!****m* fortran/fmksens2
! NAME
!  fmksens2 --- compute sensitivity matrix
! FUNCTION
!    Construct the sensitivity matrix at the end of the procedure.
!    This is only for the regression coefficients. The full sensitivity
!    is computed by fmksens2full
!    See also the R function makeSens2.
! INPUTS
!   Too many to list, but analogous to those in other routines in
!   this file. Inputs appended with d refer to process 2.
! OUTPUTS
!   Smat    matrix of size ncovs1+ncovs2 x ncovs2+ncovs2
!           sensitivity matrix
! SYNOPSIS
      subroutine fmksens2(Smat,index1,index2,Z,Zd,alphars,alpharsd,
     $ as,asd,betahat,betadhat,times1,times2,pi,qi,ri,si,wi,wij,zij,
     $ sig2,sig2d,nu2,nu2d,theta,ncovs1,ncovs2,nr,ns,nsd,
     $ d1,d2,m,Ji,maxj)
! SOURCE
     
! output of frailty estimation
      double precision pi(m),qi(m),ri(m),si(m),wi(m),
     $  wij(m,maxj),zij(m,maxj)
! output of dispersion estimates
      double precision sig2,sig2d,nu2,nu2d,theta
! d1,d2: length of data matrix
! m, Ji, maxj: # of clusters, individuals, max individuals
      integer d1,d2,m,Ji(m),maxj,ncovs1,ncovs2,nr,ns,nsd
! Sensitivity matrix (for output)
      double precision Smat(ncovs1+ncovs2,ncovs1+ncovs2)
! matrix of alpharh values
      double precision as(nr,ns), asd(nr,nsd)
      double precision alphars(nr,ns),alpharsd(nr,nsd)
! matrix of coefficient estimates
      double precision betahat(ncovs1),betadhat(ncovs2)
! matrix of indices
      integer index1(d1,6),index2(d2,6)
! covariates: matrix of covariates for each index entry
      double precision Z(d1,ncovs1)
      double precision Zd(d2,ncovs2)
      double precision times1(d1),times2(d2)
! variables to hold the current values
      double precision bZd,bZ,mu,eta
      double precision w,wj,time,thisA
      integer i,ind,r,iind0,iind1,iind0d,iind1d,s
! Constants for column indices
      parameter(icol=1, jcol=2, kcol=3, ircol=4, ismincol=5, ismaxcol=6)
      
      double precision S11(ncovs1,ncovs1)
      double precision S12(ncovs1,ncovs2)
      double precision S13(ncovs2,ncovs2)
      double precision S21(ncovs1),S22(ncovs2)
      double precision S31(ncovs1),S32(ncovs2)      
      double precision S2(ncovs1+ncovs2)
      double precision S3(ncovs1+ncovs2)
      double precision Sii(ncovs1+ncovs2,ncovs1+ncovs2)
!      double precision Sinv(g+gd,g+gd)
      
      double precision Sm(maxj)
      double precision Se(maxj)
      double precision Smx(maxj,ncovs1)
      double precision Sex(maxj,ncovs2)
      
      double precision phi(maxj),Z2(ncovs1),Z2d(ncovs2)
      
! indices for the start and end of the i portion of the data matrix      
      iind0=0
      iind1=0
      iind0d=0
      iind1d=0
      
      call dzero(Smat,(ncovs1+ncovs2)*(ncovs1+ncovs2))
      
! Loop over the values of i
      do 200 i=1,m
      
        call dzero(S11,ncovs1*ncovs1)
        call dzero(S12,ncovs1*ncovs2)
        call dzero(S13,ncovs2*ncovs2)
        call dzero(S21,ncovs1)
        call dzero(S31,ncovs1)
        call dzero(S22,ncovs2)
        call dzero(S32,ncovs2)       
        call dzero(S2,ncovs1+ncovs2)
        call dzero(S3,ncovs1+ncovs2)
        call dzero(Sm,maxj)
        call dzero(Se,maxj)
        call dzero(Smx,maxj*ncovs1)
        call dzero(Sex,maxj*ncovs2)
        call dzero(Z2,ncovs1)
        call dzero(Z2d,ncovs2)

        iind0=iind1+1
        iind1=iind1+1
        iind0d=iind1d+1
        iind1d=iind1d+1

      
!     while loop to find start and end indices in ij for cluster i
 210    if((index1(iind1,icol).eq.i).and.(iind1.le.d1)) then
            iind1=iind1+1
            goto 210
        endif      
 211    if((index2(iind1d,icol).eq.i).and.(iind1d.le.d2)) then
            iind1d=iind1d+1
            goto 211
        endif     

!       Compute the Sm, Se vectors and Smx,Sex matrix
!       Sm=Sum(mu_ij*) for j=1..Ji
!       Smx=Sum(mu_ijkh*x_ijkh) for j=1..Ji
        iind1=iind1-1                
        do 220 ind=iind0,iind1
            j=index1(ind,jcol)
            smin=index1(ind,ismincol)
            smax=index1(ind,ismaxcol)
            r=index1(ind,ircol)
            time=times1(ind)
            bZ=ddot(ncovs1,Z(ind,1),d1,betahat,1)
            if(smax .GT. 0) then
                do 230 s=smin,smax
                    thisA=0.d0
                    call A(thisA,time,as,r,s,nr,ns)
                    mu=exp(bZ)*alphars(r,s)*thisA
                    Sm(j)=Sm(j)+mu                        
                    call daxpy(ncovs1,mu,Z(ind,1),d1,Smx(j,1),maxj)
 230            continue
            endif  
 220    continue 
        iind1d=iind1d-1                
        do 221 ind=iind0d,iind1d              
            j=index2(ind,jcol)
            smin=index2(ind,ismincol)
            smax=index2(ind,ismaxcol)
            r=index2(ind,ircol)
            time=times2(ind)
            bZ=ddot(ncovs2,Zd(ind,1),d2,betadhat,1)
            if((smax .GT. 0)) then
                do 240 s=smin,smax
                    thisA=0.d0
                    call A(thisA,time,asd,r,s,nr,nsd)
                    eta=exp(bZ)*alpharsd(r,s)*thisA
                    Se(j)=Se(j)+eta                        
                    call daxpy(ncovs2,eta,Zd(ind,1),d2,Sex(j,1),maxj)
 240            continue
            endif   
 221    continue 
        
 
! Compute phi
        do 300 j=1,Ji(i)
            phi(j)=(1.0d0+nu2*Sm(j))*(1.0d0+nu2d*Se(j))
     $      -theta*theta*Sm(j)*Se(j)   
 300    continue    
            
      
! Compute the various components of the matrix (pass 2)
        do 400 ind=iind0,iind1
            j=index1(ind,jcol)
            smin=index1(ind,ismincol)
            smax=index1(ind,ismaxcol)
            r=index1(ind,ircol)
            time=times1(ind)
            bZ=ddot(ncovs1,Z(ind,1),d1,betahat,1)
            call dcopy(ncovs1,Z(ind,1),d1,Z2(1),1)
            
            if(smax .GT. 0) then
                do 410 s=smin,smax
                    thisA=0.d0
                    call A(thisA,time,as,r,s,nr,ns)
                    mu=exp(bZ)*alphars(r,s)*thisA
      
                    call dger(ncovs1,ncovs1,mu,Z2,1,Z2,1,S11,ncovs1)
                    
                    w=1.0d0/(1.0d0+wij(i,j)*Sm(j))
                    call daxpy(ncovs1,w*mu,Z(ind,1),d1,S21(1),1)
                    
                    w= -theta*Se(j)/phi(j)
                    call daxpy(ncovs1,w*mu,Z(ind,1),d1,S31(1),1)
 410            continue
            endif
 400    continue
        do 401 ind=iind0d,iind1d
            j=index2(ind,jcol)
            smin=index2(ind,ismincol)
            smax=index2(ind,ismaxcol)
            r=index2(ind,ircol)
            time=times2(ind)
            bZ=ddot(ncovs2,Zd(ind,1),d2,betadhat,1)
            call dcopy(ncovs,Zd(ind,1),d2,Z2(1),1)
            if((smax .GT. 0)) then
                do 420 s=smin,smax
                    thisA=0.d0
                    call A(thisA,time,asd,r,s,nr,nsd)
                    eta=exp(bZ)*alpharsd(r,s)*thisA

                    call dger(ncovs2,ncovs2,eta,Z2,1,Z2,1,S13,ncovs2)
                    
                    w= -theta*Sm(j)/phi(j)
                    call daxpy(ncovs2,w*eta,Zd(ind,1),d2,S22(1),1)
                    
                    w=1.0d0/(1.0d0+zij(i,j)*Se(j))
                    call daxpy(ncovs2,w*eta,Zd(ind,1),d2,S32(1),1)
 420            continue
            endif
 401    continue
  
          

! Finish computation of the components via appropriate rank 1 transforms
        do 500 j=1,Ji(i)
         w=-wij(i,j)/(1.0d0+wij(i,j)*Sm(j))
       call dger(ncovs1,ncovs1,w,Smx(j,1),maxj,Smx(j,1),maxj,S11,ncovs1)
         w=-theta/phi(j)
       call dger(ncovs1,ncovs2,w,Smx(j,1),maxj,Sex(j,1),maxj,S12,ncovs1)
         w=-zij(i,j)/(1.0d0+zij(i,j)*Se(j))
       call dger(ncovs2,ncovs2,w,Sex(j,1),maxj,Sex(j,1),maxj,S13,ncovs2)
 500    continue 
        call dcopy(ncovs1,S21,1,S2,1)
        call dcopy(ncovs2,S22,1,S2(ncovs1+1),1)
        call dcopy(ncovs1,S31,1,S3,1)
        call dcopy(ncovs2,S32,1,S3(ncovs1+1),1)
 
 
      
! Concatenate the matrices S11,S12,S13 correctly into Si
! At the same time, apply the minus sign
        do 610 jj=1,ncovs1
            do 620 ii=1,ncovs1
                Sii(ii,jj)=-S11(ii,jj)
 620        continue
 610    continue        
        do 630 jj=1,ncovs2
            do 640 ii=1,ncovs1
                Sii(ii,jj+ncovs1)=-S12(ii,jj)
                Sii(jj+ncovs1,ii)=-S12(ii,jj)
 640        continue
 630    continue
        do 650 jj=1,ncovs2
            do 660 ii=1,ncovs2
                Sii(ii+ncovs1,jj+ncovs1)=-S13(ii,jj)
 660        continue
 650    continue
 

! Compute the ith sensitivity matrix with several rank 1 operations
        wj=wi(i)*sig2*(1.0d0+sig2d*si(i))
        call dsyr('u',ncovs1+ncovs2,wj,S2,1,Sii,ncovs1+ncovs2)
        wj=-wi(i)*sig2*sig2d*qi(i)
        call dsyr2('u',ncovs1+ncovs2,wj,S2,1,S3,1,Sii,ncovs1+ncovs2)
        wj=wi(i)*sig2d*(1.0d0+sig2*pi(i))
        call dsyr('u',ncovs1+ncovs2,wj,S3,1,Sii,ncovs1+ncovs2)
        

! Add this matrix to the total
        call daxpy((ncovs1+ncovs2)*(ncovs1+ncovs2),1.d0,Sii,1,Smat,1)
 
 200  continue    

! Make the matrix symmetrical again (copy U to L)
      do 700 j=1,ncovs1+ncovs2-1
         do 710 i=j+1,ncovs1+ncovs2
            Smat(i,j)=Smat(j,i)
 710     continue
 700  continue      
      
      end
!******

!cccccccccccccccccccccccccccccccccccccccc     

!****m* fortran/fmksens2full
! NAME
!  fmksens2full --- compute sensitivity matrix
! FUNCTION
!    Construct the full sensitivity matrix at the end of the procedure.
!    This computes sensitivity to both regression coefficient and
!    baseline hazard parameters
!    See also the R function makeSens2.
! INPUTS
!   Too many to list, but analogous to those in other routines in
!   this file. Inputs appended with d refer to process 2.
!   np, npd     number of parrameters for each of the processes
! OUTPUTS
!   Smat    matrix of size ncovs1+ncovs2 x ncovs2+ncovs2
!           sensitivity matrix
! SYNOPSIS
      subroutine fmksens2full(Smat,index1,index2,Z,Zd,alphars,alpharsd,
     $ as,asd,betahat,betadhat,times1,times2,pi,qi,ri,si,wi,wij,zij,
     $ sig2,sig2d,nu2,nu2d,theta,ncovs1,ncovs2,nr,ns,nsd,np,npd,
     $ d1,d2,m,Ji,maxj,Kcum,Kdcum)
! SOURCE
     
! output of frailty estimation
      double precision pi(m),qi(m),ri(m),si(m),wi(m),
     $  wij(m,maxj),zij(m,maxj)
! output of dispersion estimates
      double precision sig2,sig2d,nu2,nu2d,theta
! d1,d2: length of data matrix
! m, Ji, maxj: # of clusters, individuals, max individuals
      integer d1,d2,m,Ji(m),maxj,ncovs1,ncovs2,nr,ns,nsd,np,npd
! Sensitivity matrix (for output)
      double precision as(nr,ns), asd(nr,nsd)
      double precision Smat(np+npd,np+npd)
! matrix of alpharh values
      double precision alphars(nr,ns),alpharsd(nr,nsd)
! matrix of coefficient estimates
      double precision betahat(ncovs1),betadhat(ncovs2)
! matrix of indices
      integer index1(d1,6),index2(d2,6)
! covariates: matrix of covariates for each index entry
      double precision Z(d1,ncovs1)
      double precision Zd(d2,ncovs2)
      double precision times1(d1),times2(d2)
      integer Kcum(nr+1), Kdcum(nr+1)
! variables to hold the current values
      double precision bZd,bZ,mu,eta,w,wj,time,thisA
      integer i,ind,r,iind0,iind1,iind0d,iind1d,s,bst,bdst,sind
!      real*4 timex(2),ExecTime(3),LastTime(3),tottime,tottimex

! Constants for column indices
      parameter(icol=1, jcol=2, kcol=3, ircol=4, ismincol=5, ismaxcol=6)
      
      double precision S11(np,np)
      double precision S12(np,npd)
      double precision S13(npd,npd)
      double precision S21(np),S22(npd)
      double precision S31(np),S32(npd)      
      double precision S2(np+npd)
      double precision S3(np+npd)
      double precision Sii(np+npd,np+npd)
!      double precision Sinv(g+gd,g+gd)
      
      double precision Sm(maxj)
      double precision Se(maxj)
      double precision Smx(maxj,np)
      double precision Sex(maxj,npd)
      
      double precision phi(maxj),Z2(np),Z2d(npd)
      
! indices for the start and end of the i portion of the data matrix      
      iind0=0
      iind1=0
      iind0d=0
      iind1d=0
      
      bst=Kcum(nr+1)
      bdst=Kdcum(nr+1)
     
!      call etime(timex,tottimex)      
!      call etime(timex,tottime)      
!      LastTime(1)    = tottimex
!      LastTime(2)   = timex(1)
!      LastTime(3) = timex(2)
            
      call dzero(Smat,(np+npd)*(np+npd))
      
! Loop over the values of i
      do 200 i=1,m
      
!       if(i.eq.1) then
!          call etime(timex,tottimex)
!          ExecTime(1) = tottimex-LastTime(1)
!          ExecTime(2) = timex(1)-LastTime(2)
!          ExecTime(3) = timex(2)-LastTime(3)
!          LastTime(1)=tottimex
!          LastTime(2)=timex(1)
!          LastTime(3)=timex(2)
!          call realpr('Startup',-1,ExecTime,3)
!       endif
       
        call dzero(S11,np*np)
        call dzero(S12,np*npd)
        call dzero(S13,npd*npd)
        call dzero(S21,np)
        call dzero(S31,np)
        call dzero(S22,npd)
        call dzero(S32,npd)       
        call dzero(S2,np+npd)
        call dzero(S3,np+npd)
        call dzero(Sm,maxj)
        call dzero(Se,maxj)
        call dzero(Smx,maxj*np)
        call dzero(Sex,maxj*npd)
        call dzero(Z2,np)
        call dzero(Z2d,npd)

        iind0=iind1+1
        iind1=iind1+1
        iind0d=iind1d+1
        iind1d=iind1d+1
            
!     while loop to find start and end indices in ij for cluster i
 210    if((index1(iind1,icol).eq.i).and.(iind1.le.d1)) then
            iind1=iind1+1
            goto 210
        endif      
 211    if((index2(iind1d,icol).eq.i).and.(iind1d.le.d2)) then
            iind1d=iind1d+1
            goto 211
        endif     

!       Compute the Sm, Se vectors and Smx,Sex matrix
!       Sm=Sum(mu_ij*) for j=1..Ji
!       Smx=Sum(mu_ijkh*x_ijkh) for j=1..Ji
        iind1=iind1-1                
        do 220 ind=iind0,iind1
            j=index1(ind,jcol)
            smin=index1(ind,ismincol)
            smax=index1(ind,ismaxcol)
            r=index1(ind,ircol)
            time=times1(ind)
            bZ=ddot(ncovs1,Z(ind,1),d1,betahat,1)
            if(smax .GT. 0) then
                do 230 s=smin,smax
                    thisA=0.d0
                    call A(thisA,time,as,r,s,nr,ns)
                    mu=exp(bZ)*alphars(r,s)*thisA
                    Sm(j)=Sm(j)+mu 
                    Smx(j,Kcum(r)+s-1)=Smx(j,Kcum(r)+s-1)+mu        
                    call daxpy(ncovs1,mu,Z(ind,1),d1,Smx(j,bst),maxj)
 230            continue
            endif  
 220    continue 
        iind1d=iind1d-1                
        do 221 ind=iind0d,iind1d              
            j=index2(ind,jcol)
            smin=index2(ind,ismincol)
            smax=index2(ind,ismaxcol)
            r=index2(ind,ircol)
            time=times2(ind)
            bZ=ddot(ncovs2,Zd(ind,1),d2,betadhat,1)
            if((smax .GT. 0)) then
                do 240 s=smin,smax
                    thisA=0.d0
                    call A(thisA,time,asd,r,s,nr,nsd)
                    eta=exp(bZ)*alpharsd(r,s)*thisA
                    Se(j)=Se(j)+eta                        
                    Sex(j,Kdcum(r)+s-1)=Sex(j,Kdcum(r)+s-1)+eta        
                    call daxpy(ncovs2,eta,Zd(ind,1),d2,Sex(j,bdst),maxj)
 240            continue
            endif   
 221    continue 
        
 
! Compute phi
        do 300 j=1,Ji(i)
            phi(j)=(1.0d0+nu2*Sm(j))*(1.0d0+nu2d*Se(j))
     $      -theta*theta*Sm(j)*Se(j)   
 300    continue    
                   
! Compute the various components of the matrix (pass 2)
        do 400 ind=iind0,iind1
            j=index1(ind,jcol)
            smin=index1(ind,ismincol)
            smax=index1(ind,ismaxcol)
            r=index1(ind,ircol)
            time=times1(ind)
            bZ=ddot(ncovs1,Z(ind,1),d1,betahat,1)
!            call dcopy(ncovs,Z(ind,1),d1,Z2(bst),1)
            
            if(smax .GT. 0) then
                do 410 s=smin,smax
                    thisA=0.d0
                    call A(thisA,time,as,r,s,nr,ns)
                    mu=exp(bZ)*alphars(r,s)*thisA
                    sind=Kcum(r)+s-1
                    
                    S11(sind,sind)=S11(sind,sind)+mu
                    call daxpy(ncovs1,mu,Z(ind,1),d1,S11(sind,bst),np)
!                    call daxpy(ncovs,mu,Z(ind,1),d1,S11(bst,sind),1)
!                    call dger(np,np,mu,Z2,1,Z2,1,S11,np)
                    call dsyr('U',ncovs1,mu,Z(ind,1),d1,S11(bst,bst),np)
                    
                    w=1.0d0/(1.0d0+wij(i,j)*Sm(j))
                    S21(sind)=S21(sind)+w*mu
                    call daxpy(ncovs1,w*mu,Z(ind,1),d1,S21(bst),1)
                    
                    w= -theta*Se(j)/phi(j)
                    S31(sind)=S31(sind)+w*mu
                    call daxpy(ncovs1,w*mu,Z(ind,1),d1,S31(bst),1)
 410            continue
            endif
 400    continue
        do 401 ind=iind0d,iind1d
            j=index2(ind,jcol)
            smin=index2(ind,ismincol)
            smax=index2(ind,ismaxcol)
            r=index2(ind,ircol)
            time=times2(ind)
            bZ=ddot(ncovs2,Zd(ind,1),d2,betadhat,1)
!            call dcopy(ncovs,Zd(ind,1),d2,Z2d(bdst),1)
            if((smax .GT. 0)) then
                do 420 s=smin,smax
                    thisA=0.d0
                    call A(thisA,time,asd,r,s,nr,nsd)
                    eta=exp(bZ)*alpharsd(r,s)*thisA
                    sind=Kdcum(r)+s-1

                    S13(sind,sind)=S13(sind,sind)+eta
                  call daxpy(ncovs2,eta,Zd(ind,1),d2,S13(sind,bdst),npd)
!                  call daxpy(ncovs,eta,Zd(ind,1),d2,S13(bdst,sind),1)
                    call dger(npd,npd,eta,Z2d,1,Z2d,1,S13,npd)
              call dsyr('U',ncovs2,eta,Zd(ind,1),d2,S13(bdst,bdst),npd)

                    
                    w= -theta*Sm(j)/phi(j)
                    S22(sind)=S22(sind)+w*eta
                    call daxpy(ncovs2,w*eta,Zd(ind,1),d2,S22(bdst),1)
                    
                    w=1.0d0/(1.0d0+zij(i,j)*Se(j))
                    S32(sind)=S32(sind)+w*eta
                    call daxpy(ncovs2,w*eta,Zd(ind,1),d2,S32(bdst),1)
 420            continue
            endif
 401    continue
  
  
! Finish computation of the components via appropriate rank 1 transforms
        do 500 j=1,Ji(i)
             w=-wij(i,j)/(1.0d0+wij(i,j)*Sm(j))
!             call dger(np,np,w,Smx(j,1),maxj,Smx(j,1),maxj,S11,np)
             call dsyr('U',np,w,Smx(j,1),maxj,S11,np)
             w=-theta/phi(j)
             call dger(np,npd,w,Smx(j,1),maxj,Sex(j,1),maxj,S12,np)
             w=-zij(i,j)/(1.0d0+zij(i,j)*Se(j))
!             call dger(npd,npd,w,Sex(j,1),maxj,Sex(j,1),maxj,S13,npd)
             call dsyr('U',npd,w,Sex(j,1),maxj,S13,npd)

 500    continue 
        call dcopy(np,S21,1,S2,1)
        call dcopy(npd,S22,1,S2(np+1),1)
        call dcopy(np,S31,1,S3,1)
        call dcopy(npd,S32,1,S3(np+1),1)
       
! Concatenate the matrices S11,S12,S13 correctly into Si
! At the same time, apply the minus sign
        do 610 jj=1,np
            do 620 ii=1,jj
                Sii(ii,jj)=-S11(ii,jj)
 620        continue
 610    continue        
        do 630 jj=1,npd
            do 640 ii=1,np
                Sii(ii,jj+np)=-S12(ii,jj)
!                Sii(jj+np,ii)=-S12(ii,jj)
 640        continue
 630    continue
        do 650 jj=1,npd
            do 660 ii=1,jj
                Sii(ii+np,jj+np)=-S13(ii,jj)
 660        continue
 650    continue
 

! Compute the ith sensitivity matrix with several rank 1 operations
        wj=wi(i)*sig2*(1.0d0+sig2d*si(i))
        call dsyr('u',np+npd,wj,S2,1,Sii,np+npd)
        wj=-wi(i)*sig2*sig2d*qi(i)
        call dsyr2('u',np+npd,wj,S2,1,S3,1,Sii,np+npd)
        wj=wi(i)*sig2d*(1.0d0+sig2*pi(i))
        call dsyr('u',np+npd,wj,S3,1,Sii,np+npd)
               
! Add this matrix to the total
        call daxpy((np+npd)*(np+npd),1.d0,Sii,1,Smat,1)
 
 200  continue  
   
! Make the matrix symmetrical again (copy U to L)
      do 700 j=1,np+npd-1
         do 710 i=j+1,np+npd
            Smat(i,j)=Smat(j,i)
 710     continue
 700  continue      
 
      end
!******

!cccccccccccccccccccccccccccccccccccccccc     

!****m* fortran/fmkstderr
! NAME
!  fmkstderr --- compute standard errors
! FUNCTION
!    Compute regressin parameter standard errors using the full
!    sensitivity matrix.
!    See also the R function fmkstderr.
! INPUTS
!   Too many to list, but analogous to those in other routines in
!   this file. Inputs appended with d refer to process 2.
!   B       np+npd vector with 1 in the position of coefficients
! OUTPUTS
!   Smat    matrix of size ncovs1+ncovs2 x ncovs2+ncovs2
!           sensitivity matrix
!   B       standard errors of regression coefficients in positions
!           that were 1 in the input
! SYNOPSIS
      subroutine fmkstderr(Smat,B,index1,index2,Z,Zd,alphars,alpharsd,
     $ as,asd,betahat,betadhat,times1,times2,pi,qi,ri,si,wi,wij,zij,
     $ sig2,sig2d,nu2,nu2d,theta,ncovs1,ncovs2,nr,ns,nsd,np,npd,
     $ d1,d2,m,Ji,maxj,Kcum,Kdcum)
! SOURCE
      double precision pi(m),qi(m),ri(m),si(m),wi(m),
     $  wij(m,maxj),zij(m,maxj)
      double precision sig2,sig2d,nu2,nu2d,theta
      integer d1,d2,m,Ji(m),maxj,ncovs1,ncovs2,nr,ns,nsd,np,npd
      double precision Smat(np+npd,np+npd)
      double precision B(np+npd,ncovs1+ncovs2)
      double precision as(nr,ns), asd(nr,nsd)
      double precision alphars(nr,ns),alpharsd(nr,nsd)
      double precision betahat(ncovs1),betadhat(ncovs2)
      integer index1(d1,6),index2(d2,6)
      double precision Z(d1,ncovs1)
      double precision Zd(d2,ncovs2)
      double precision times1(d1),times2(d2)
      integer Kcum(nr+1), Kdcum(nr+1)
      double precision bZd,bZ,mu,eta,w,wj,time,thisA
      integer i,ind,r,iind0,iind1,iind0d,iind1d,s,bst,bdst,sind
      integer INFO

      call fmksens2full(Smat,index1,index2,Z,Zd,alphars,alpharsd,
     $ as,asd,betahat,betadhat,times1,times2,pi,qi,ri,si,wi,wij,zij,
     $ sig2,sig2d,nu2,nu2d,theta,ncovs1,ncovs2,nr,ns,nsd,np,npd,
     $ d1,d2,m,Ji,maxj,Kcum,Kdcum)
 
      call dscal((np+npd)*(np+npd),-1.d0,Smat,1)
     
      call dposv('u',np+npd,ncovs1+ncovs2,Smat,np+npd,B,
     $ np+npd,INFO)
      
      end
!******

!cccccccccccccccccccccccccccccccccccccccc     
     
!****m* fortran/dzero
! NAME
!  dzero --- zero out a double vector
! FUNCTION
!    Clears the values in a double vector and replaces them with 0.
! INPUTS
!   x       a double vector
!   len     length of x
! SYNOPSIS
      subroutine dzero(x,len)
! SOURCE
      integer len,i
      double precision x(len)
      do 100 i=1,len
        x(i)=0.d0
 100    continue
      end     
!******
