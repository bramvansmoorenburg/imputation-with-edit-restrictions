

sampleimp=function(data,fun,edits,totals,weight=NULL,init.values=NULL,maxit=5) {
  
  
  
  if (is.null(weight)) {
    
    wts=rep(1,nrow(data))
    
  } else {
    
    wts=weight
    
  }
  
  
  
  if (is.editmatrix(edits)) {
    
    rules=edits
    
  } else {
    
    rules=editmatrix(edits)
    
  }
  
  errnames=colnames(rules)[!colnames(rules)%in%c("CONSTANT",colnames(data))]
  
  if (length(errnames)==0) {
    
    stop("(sampleimp) Edits have no error terms")
    
  }
  
  if (!any(is.na(data))) {
    
    warning("(sampleimp) There are no missing values in the dataset")
    
  }
  
  
  
  equalities=colnames(rules)[which(colnames(rules)%in%errnames&colSums(getAb(rules)!=0)>1)]
  
  
  
  totmiss=totals-colSums(data*wts,na.rm=TRUE)
  
  
  
  solvedmaterr=function(crow,x1) {
    
    edmatr=substValue(rules,colnames(data),x1[crow,])
    
    solved=getAb(edmatr)[,c(errnames,"CONSTANT")]
    
    return(sapply(1:(ncol(solved)-1),function(x){tmp=solved[,"CONSTANT"]/solved[,x];tmp=na.omit(tmp);max(tmp[tmp!=Inf])}))
    
  }
  
  
  
  chooseimp=function(crow,old=oldimprows,new=improws,better=impbetter) {
    
    if (any(is.na(old[crow,]))) {
      
      stop("(sampleimp) NA's in initial values")
      
    } else if (any(is.na(new[crow,]))) {
      
      stop("(sampleimp) Imputation function introduced NA's in dataset")
      
    } else if (all(better[crow,])) {
      
      return(TRUE)
      
    } else {
      
      return(FALSE)
      
    }
    
  }
  
  
  
  if (is.null(init.values)) {
    
    if (length(formals(fun))!=1) {
      
      stop("(sampleimp) Imputation function requires too many arguments")
      
    }
    
    inits=fun(data)
    
  } else {
    
    inits=init.values
    
  }
  
  
  
  narows=which(apply(data,1,function(x){any(is.na(x))}))
  
  navars=which(apply(data,2,function(x){any(is.na(x))}))
  
  oldimprows=inits[narows,]
  
  oldimperrs=t(sapply(1:nrow(oldimprows),solvedmaterr,x1=oldimprows))
  
  temperr=((colSums(oldimprows*wts[narows]*is.na(data[narows,]))-totmiss)/totmiss)[navars] #note, this will give errors if nonnegative in freeimp is set to false and totmiss in one variable is exactly zero while there are still missing values
  
  startp=ifelse(length(formals(fun))==1,2,1)
  
  
  
  for (i in startp:maxit) {
    
    if (length(formals(fun))==1) {
      
      impdata=fun(data)
      
    } else if (length(formals(fun))==2) {
      
      impdata=fun(data,inits)
      
    } else {
      
      stop("(sampleimp) Imputation function has too many or too few arguments")
      
    }
    
    improws=impdata[narows,]
    
    imperrs=t(sapply(1:nrow(improws),solvedmaterr,x1=improws))
    
    impbetter=matrix(FALSE,nrow = nrow(imperrs),ncol = ncol(imperrs))
    
    for (j in 1:ncol(imperrs)) {
      
      if (errnames[j]%in%equalities) {
        
        impbetter[,j]=abs(imperrs[,j])<=abs(oldimperrs[,j])
        
      } else {
        
        impbetter[,j]=imperrs[,j]<=oldimperrs[,j]|oldimperrs[,j]<=0&imperrs[,j]<=0
        
      }
      
    }
    
    impchoice=sapply(1:nrow(improws),chooseimp)
    
    impchange=t(t((improws-oldimprows)*wts[narows])/totmiss)[,navars,drop=FALSE]
    
    for (k in sample(which(impchoice))) {
      
      if (sum(na.omit(impchange[k,])*sign(na.omit(temperr)))<=0) {
        
        temperr=temperr+impchange[k,]
        
      } else {
        
        impchoice[k]=FALSE
        
      }
      
    }
    
    oldimprows=improws*impchoice+oldimprows*!impchoice
    
    oldimperrs=imperrs*impchoice+oldimperrs*!impchoice
    
  }
  
  impdata=data
  
  impdata[narows,]=oldimprows
  
  
  
  return(impdata)
  
}
