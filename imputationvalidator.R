

impvalidator=function(imputed,dataset,truevals,rules,frac=0,maxits=0,weight=NULL) {
  
  
  
  if (any(is.na(imputed))) {
    
    stop("(impvalidator) NA's in imputed dataset")
    
  }
  
  
  
  if (!is.null(weight)) {
    
    wts=weight
    
  } else {
    
    wts=rep(1,nrow(truevals))
    
  }
  
  
  
  totals=colSums(truevals*wts)
  
  totmiss2=totals-colSums(dataset*wts,na.rm = TRUE)
  
  
  
  truevars=numeric(0)
  
  truerows=numeric(0)
  
  trueweights=numeric(0)
  
  for (i in 1:ncol(truevals)) {
    
    tmp=which(is.na(dataset[,i]))
    
    truerows=c(truerows,tmp)
    
    truevars=c(truevars,rep(i,length(tmp)))
    
    trueweights=c(trueweights,wts[tmp])
    
  }
  
  validation=rbind(true=truevals[is.na(dataset)],variable=truevars,row=truerows,weights=trueweights,imputed[is.na(dataset)])
  
  
  
  vardif1=sapply(1:ncol(dataset),function(x){y=which(validation[2,]==x);return(sum(abs(validation[5,y]-validation[1,y])*validation[4,y]))})
  
  reldif=sum(na.omit(vardif1/totmiss2))
  
  vardif=c(vardif1,-1,vardif1/totmiss2)
  
  names(vardif)=c(colnames(dataset),"-",colnames(dataset))
  
  
  
  covdif1=log(abs(cov(imputed)/cov(truevals)))
  
  covdif2=which(!lower.tri(covdif1),arr.ind = TRUE)
  
  covdifnam=apply(covdif2,1,function(x){paste0(colnames(covdif1)[x[1]],"-",colnames(covdif1)[x[2]])})
  
  covdif=covdif1[!lower.tri(covdif1)]
  
  names(covdif)=covdifnam
  
  covdifmed=median(abs(covdif))
  
  
  
  totmiss=(totals-colSums(imputed*wts))/totmiss2
  
  
  
  meansn=0.1 #beware: unweighted means
  
  meansN=round(meansn*nrow(truevals))
  
  truemeans=colMeans(truevals)
  
  trues=apply(truevals,2,function(x){sd(x)/sqrt(meansN)})
  
  impsamples=replicate(1e4,{tmp=sample(nrow(imputed),meansN);colMeans(imputed[tmp,])})
  
  meansconf=sapply(1:ncol(truevals),function(x){mean(impsamples[x,]>=truemeans[x]-qnorm(0.975)*trues[x]&impsamples[x,]<=truemeans[x]+qnorm(0.975)*trues[x])})
  
  
  
  weightedks=function(a,b,wts=NULL) {
    
    library(Hmisc,quietly = TRUE)
    
    if (is.null(wts)) {
      
      w=rep(1,length(a))
      
    } else {
      
      w=wts
      
    }
    
    dwecdf=function(x,y1,y2,w) ecdf(Ecdf(y1,weights =  w,pl=FALSE)$x)(x)-ecdf(Ecdf(y2,weights = w,pl=FALSE)$x)(x)
    
    curve2=function (expr, from = NULL, to = NULL, n = 101, add = FALSE,
                     
                     type = "l", xname = "x", xlab = xname, ylab =
                       
                       NULL, log = NULL,
                     
                     xlim = NULL, ...)
      
    {
      
      
      
      sexpr <- substitute(expr)
      
      if (is.name(sexpr)) {
        
        expr <- call(as.character(sexpr), as.name(xname))
        
      }
      
      else {
        
        if (!((is.call(sexpr) || is.expression(sexpr)) && xname %in%
              
              all.vars(sexpr)))
          
          stop(gettextf("'expr' must be a function, or a call or an expression containing '%s'", xname), domain = NA)
        
        expr <- sexpr
        
      }
      
      if (dev.cur() == 1L && !identical(add, FALSE)) {
        
        warning("'add' will be ignored as there is no existing plot")
        
        add <- FALSE
        
      }
      
      addF <- identical(add, FALSE)
      
      if (is.null(ylab))
        
        ylab <- deparse(expr)
      
      if (is.null(from) || is.null(to)) {
        
        xl <- if (!is.null(xlim))
          
          xlim
        
        else if (!addF) {
          
          pu <- par("usr")[1L:2L]
          
          if (par("xaxs") == "r")
            
            pu <- extendrange(pu, f = -1/27)
          
          if (par("xlog"))
            
            10^pu
          
          else pu
          
        }
        
        else c(0, 1)
        
        if (is.null(from))
          
          from <- xl[1L]
        
        if (is.null(to))
          
          to <- xl[2L]
        
      }
      
      lg <- if (length(log))
        
        log
      
      else if (!addF && par("xlog"))
        
        "x"
      
      else ""
      
      if (length(lg) == 0)
        
        lg <- ""
      
      if (grepl("x", lg, fixed = TRUE)) {
        
        if (from <= 0 || to <= 0)
          
          stop("'from' and 'to' must be > 0 with log=\"x\"")
        
        x <- exp(seq.int(log(from), log(to), length.out = n))
        
      }
      
      else x <- seq.int(from, to, length.out = n)
      
      ll <- list(x = x)
      
      names(ll) <- xname
      
      y <- eval(expr, envir = ll, enclos = parent.frame())
      
      if (length(y) != length(x))
        
        stop("'expr' did not evaluate to an object of length 'n'")
      
      invisible(list(x = x, y = y))
      
    }
    
    diferr=curve2(dwecdf(x,a,b,w),from=min(a,b),to=max(a,b))
    
    return(max(abs(diferr$y)))
    
  }
  
  
  
  kstests=sapply(1:ncol(dataset),function(x){weightedks(truevals[,x],imputed[,x],wts=wts)})
  
  names(kstests)=colnames(dataset)
  
  
  
  errnames=colnames(rules)[!colnames(rules)%in%c("CONSTANT",colnames(dataset))]
  
  if (length(errnames)==0) {
    
    stop("(impvalidator) Rules have wrong format")
    
  }
  
  equalities=colnames(rules)[which(colnames(rules)%in%errnames&colSums(getAb(rules)!=0)>1)]
  
  solvedmaterr=function(crow,x1) {
    
    edmatr=substValue(rules,colnames(dataset),x1[crow,])
    
    solved=getAb(edmatr)[,c(errnames,"CONSTANT")]
    
    return(sapply(1:(ncol(solved)-1),function(x){tmp=solved[,"CONSTANT"]/solved[,x];tmp=na.omit(tmp);max(tmp[tmp!=Inf])}))
    
  }
  
  allnas=which(apply(dataset,1,function(x){any(is.na(x))}))
  
  finalerrs=t(sapply(1:length(allnas),solvedmaterr,x1=imputed[allnas,]))
  
  editerr=sum(sapply(1:length(errnames),function(x){ifelse(errnames[x]%in%equalities,sum(abs(finalerrs[,x])*wts[allnas]),sum(finalerrs[,x][finalerrs[,x]>0]*wts[allnas][finalerrs[,x]>0]))}))
  
  
  
  return(list(reldif=reldif,totmiss=totmiss,covdifmed=covdifmed,covdif=covdif,vardif=vardif,editerr=editerr,meansconf=meansconf,kstests=kstests,frac=frac,maxits=maxits))
  
}



validatecleaner=function(impvalidators) {
  
  reldif1=unlist(lapply(impvalidators,function(x){x[[1]]}))
  
  totmiss1=do.call(rbind,lapply(impvalidators,function(x){x[[2]]}))
  
  covdifmed1=unlist(lapply(impvalidators,function(x){x[[3]]}))
  
  covdif1=do.call(rbind,lapply(impvalidators,function(x){x[[4]]}))
  
  vardif1=do.call(rbind,lapply(impvalidators,function(x){x[[5]]}))
  
  editerr1=unlist(lapply(impvalidators,function(x){x[[6]]}))
  
  meansconf1=do.call(rbind,lapply(impvalidators,function(x){x[[7]]}))
  
  kstests1=do.call(rbind,lapply(impvalidators,function(x){x[[8]]}))
  
  frac1=unlist(lapply(impvalidators,function(x){x[[9]]}))
  
  maxits1=unlist(lapply(impvalidators,function(x){x[[10]]}))
  
  return(list(reldif=reldif1,totmiss=totmiss1,covdifmed=covdifmed1,covdif=covdif1,vardif=vardif1,editerr=editerr1,meansconf=meansconf1,kstests=kstests1,frac=frac1,maxits=maxits1))
  
}



rubinvalidator=function(datalist,datamissing,truevals,response,explanatory,frac=0,maxits=0,datinpool=0,logtrans=FALSE) {
  
  
  
  library(mice)
  
  
  
  poolcoefs=function(mira,wholepopulation=FALSE) {
    
    coefs=sapply(mira$analyses,function(a){summary(a)$coefficients[,1]})
    
    serrs=sapply(mira$analyses,function(a){summary(a)$coefficients[,2]})
    
    
    
    theta=rowMeans(coefs)
    
    if (wholepopulation) {
      
      Vw=rep(0,nrow(coefs))
      
    } else {
      
      Vw=rowMeans(serrs^2)
      
    }
    
    Vb=sqrt(rowSums((coefs-theta)^2)/(nrow(coefs)-1))
    
    Vt=Vw+Vb+Vb/nrow(coefs)
    
    df=(nrow(coefs)-1)*(1+nrow(coefs)*Vw/((nrow(coefs)+1)*Vb))^2
    
    return(cbind(coefs=theta,se=Vt,df=df))
    
  }
  
  
  
  datalist2=c(list(datamissing),datalist)
  
  Msets=lapply(1:length(datalist2),function(x){cbind(.imp=x-1,.id=1:nrow(datalist2[[x]]),datalist2[[x]])})
  
  Msets=do.call(rbind,Msets)
  
  Mmids=as.mids(Msets)
  
  if (logtrans) {
    
    Mmodel=with(Mmids,lm(as.formula(paste0("log(",colnames(datamissing)[response],"+1)~",paste(colnames(datamissing)[explanatory],collapse = "+")))))
    
  } else {
    
    Mmodel=with(Mmids,lm(as.formula(paste0(colnames(datamissing)[response],"~",paste(colnames(datamissing)[explanatory],collapse = "+")))))
    
  }
  
  #Mparams=matrix(as.vector(summary(pool(Mmodel))[,1:2]),nrow=1)
  
  #Mdf=summary(pool(Mmodel))[,4]
  
  Mparams=matrix(as.vector(poolcoefs(Mmodel)[,1:2]),nrow=1)
  
  Mdf=poolcoefs(Mmodel)[,3]
  
  if (class(Mparams[1])=="list") { #on some systems, unlist and as.vector work differently
    
    #Mparams=matrix(unlist(summary(pool(Mmodel))[,1:2]),nrow=1)
    
    Mparams=matrix(unlist(poolcoefs(Mmodel)[,1:2]),nrow=1)
    
  }
  
  
  
  if (logtrans) {
    
    truemod=lm(as.formula(paste0("log(",colnames(truevals)[response],"+1)~",paste(colnames(truevals)[explanatory],collapse = "+"))),data=truevals)
    
  } else {
    
    truemod=lm(as.formula(paste0(colnames(truevals)[response],"~",paste(colnames(truevals)[explanatory],collapse = "+"))),data=truevals)
    
  }
  
  truecoefs=summary(truemod)$coefficients
  
  
  
  if (nrow(na.omit(datamissing))>0) {
    
    if (logtrans) {
      
      ccmod=lm(as.formula(paste0("log(",colnames(truevals)[response],"+1)~",paste(colnames(truevals)[explanatory],collapse  = "+"))),data=na.omit(datamissing))
      
    } else {
      
      ccmod=lm(as.formula(paste0(colnames(truevals)[response],"~",paste(colnames(truevals)[explanatory],collapse  = "+"))),data=na.omit(datamissing))
      
    }
    
    cccoefs=summary(ccmod)$coefficients
    
  } else {
    
    cccoefs=truecoefs
    
    cccoefs[!is.na(truecoefs)]=NA
    
  }
  
  
  
  Mintervals=apply(Mparams,1,function(x){ints=rbind(x[1:(ncol(Mparams)/2)]+qnorm(0.975)*x[(ncol(Mparams)/2+1):ncol(Mparams)],x[1:(ncol(Mparams)/2)]-qnorm(0.975)*x[(ncol(Mparams)/2+1):ncol(Mparams)]);truecoefs[,1]<=ints[1,]&truecoefs[,1]>=ints[2,]})
  
  #Mintervals=apply(Mparams,1,function(x){ints=rbind(x[1:(ncol(Mparams)/2)]+qt(0.975,Mdf)*x[(ncol(Mparams)/2+1):ncol(Mparams)],x[1:(ncol(Mparams)/2)]-qT(0.975,Mdf)*x[(ncol(Mparams)/2+1):ncol(Mparams)]);truecoefs[,1]<=ints[1,]&truecoefs[,1]>=ints[2,]})
  
  Mparams=rbind(as.vector(truecoefs[,1:2]),as.vector(cccoefs[,1:2]),c(as.numeric(Mintervals),rep(NA,ncol(Mparams)/2)),Mparams)
  
  colnames(Mparams)=c(paste("b",rownames(summary(pool(Mmodel)))),paste("se",rownames(summary(pool(Mmodel)))))
  
  Mparams=cbind(frac=c(NA,NA,NA,frac),maxits=c(NA,NA,NA,maxits),datinpool=c(NA,NA,NA,datinpool),Mparams)
  
  return(Mparams)
  
}



rubincleaner=function(rubinlist) {
  
  metadat=do.call(rbind,lapply(rubinlist,function(x){x[4,c(1:3)]}))
  
  uniquemetadat=unique(metadat)
  
  
  
  rubinlist2=lapply(rubinlist,function(x){x[,-c(1:3)]})
  
  
  
  originalcoefs=do.call(rbind,lapply(rubinlist2,function(x){x[1,]}))
  
  originalequal=all(apply(originalcoefs,2,function(x){!any(is.na(x))&&all(x==x[1])}))
  
  imputecoefs=lapply(rubinlist2,function(x){x[c(4:nrow(x)),]})
  
  
  
  Mintervals=list(0)
  
  correctcoefs=do.call(rbind,lapply(rubinlist2,function(x){x[3,]}))
  
  for (i in 1:nrow(uniquemetadat)) {
    
    correctcoefs2=correctcoefs[apply(metadat,1,function(x){all(x==uniquemetadat[i,])}),,drop=FALSE]
    
    correctcoefs2[is.na(correctcoefs2[,ncol(correctcoefs2)]),ncol(correctcoefs2)]=1
    
    totalpools=sum(correctcoefs2[,ncol(correctcoefs2)])
    
    Mintervals[[i]]=colSums(correctcoefs2*correctcoefs2[,ncol(correctcoefs2)])/totalpools
    
    Mintervals[[i]][length(Mintervals[[length(Mintervals)]])]=totalpools
    
    Mintervals[[i]]=rbind(Mintervals[[length(Mintervals)]],do.call(rbind,imputecoefs[apply(metadat,1,function(x){all(x==uniquemetadat[i,])})]))
    
  }
  
  
  
  cccoefs=do.call(rbind,lapply(rubinlist2,function(x){x[2,]}))
  
  for (i in 1:nrow(uniquemetadat)) {
    
    cccoefs2=cccoefs[apply(metadat,1,function(x){all(x==uniquemetadat[i,])}),,drop=FALSE]
    
    ccequal=all(sapply(1:ncol(cccoefs2),function(x){all(!is.na(cccoefs2[,x]))&&all(cccoefs2[,x]==cccoefs2[1,x])}))
    
    if (ccequal) {
      
      Mintervals[[i]]=rbind(cccoefs2[1,],Mintervals[[i]])
      
    } else {
      
      Mintervals[[i]]=rbind(NA,Mintervals[[i]])
      
    }
    
  }
  
  
  
  for (i in 1:nrow(uniquemetadat)) {
    
    originalcoefs2=originalcoefs[apply(metadat,1,function(x){all(x==uniquemetadat[i,])}),,drop=FALSE]
    
    originalequal=all(sapply(1:ncol(originalcoefs2),function(x){all(!is.na(originalcoefs2[,x]))&&all(originalcoefs2[,x]==originalcoefs2[1,x])}))
    
    if (originalequal) {
      
      Mintervals[[i]]=rbind(originalcoefs2[1,],Mintervals[[i]])
      
    } else {
      
      Mintervals[[i]]=rbind(NA,Mintervals[[i]])
      
    }
    
  }
  
  
  
  names(Mintervals)=NULL
  
  
  
  for (i in 1:nrow(uniquemetadat)) {
    
    Mintervals[[i]]=cbind(frac=uniquemetadat[i,1],maxits=uniquemetadat[i,2],datinpool=uniquemetadat[i,3],Mintervals[[i]])
    
  }
  
  
  
  return(Mintervals)
  
}
