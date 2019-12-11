

getinits=function(data,rules,totals,Weights=NULL) {
  
  
  
  if (!is.null(Weights)) {
    
    wts=Weights
    
  } else {
    
    wts=rep(1,nrow(data))
    
  }
  
  
  
  errnames=colnames(rules)[!colnames(rules)%in%c("CONSTANT",colnames(data))]
  
  if (length(errnames)==0) {
    
    stop("edits in incorrect format")
    
  }
  
  equalities=colnames(rules)[which(colnames(rules)%in%errnames&colSums(getAb(rules)!=0)>1)]
  
  if (length(equalities)>0) {
    
    equalvars=sapply(equalities,function(x){abrules=getAb(rules);sapply(1:ncol(data),function(y){any(abrules[,colnames(data)][,y]!=0&abrules[,x]!=0)})})
    
    equalvars=colnames(data)[apply(equalvars,1,any)]
    
  }
  
  
  
  solvedmat=function(crow,x1,cname1,cnames1,tolfac=1) {
    
    edmat=sapply(helplist$pattern,function(b){ifelse(all(c(cnames1,cname1,errnames)%in%b)&all(b%in%c(cnames1,cname1,errnames)),TRUE,FALSE)})
    
    edmat=which(edmat)
    
    
    
    edmatr=helplist$matrices[[edmat]]
    
    if (length(equalities)>0) {
      
      edmatr=substValue(edmatr,equalities,rep((max(inits[crow,equalvars])+100)*tolfac,length(equalities)))
      
      edmatr=substValue(edmatr,errnames[!errnames%in%equalities],rep(0,length(errnames[!errnames%in%equalities])))
      
    } else {
      
      edmatr=substValue(edmatr,errnames,rep(0,length(errnames)))
      
    }
    
    
    
    solved=substValue(edmatr,cnames1,x1[crow,])
    
    solved=getAb(solved)
    
    if (length(solved)==0||all(solved[,cname1]==0)) {
      
      return(c(-Inf,Inf))
      
    }
    
    const=solved[,"CONSTANT"]*wts[crow]/solved[,cname1]
    
    if (any(const==Inf|const==-Inf)) {
      
      solved=solved[-which(const==Inf|const==-Inf),,drop=FALSE]
      
      const=const[-which(const==Inf|const==-Inf)]
      
    }
    
    if (all(solved[,cname1]<=0)) {
      
      return(c(max(const),Inf))
      
    } else if (all(solved[,cname1]>=0)) {
      
      return(c(-Inf,min(const)))
      
    } else {
      
      out=c(max(const[solved[,cname1]<0]),min(const[solved[,cname1]>0]))
      
      if (out[2]<out[1]) {
        
        tmp=out[2]
        
        out[2]=out[1]
        
        out[1]=tmp
        
      }
      
      return(out)
      
    }
    
  }
  
  
  
  totmiss=totals-colSums(data*wts,na.rm=TRUE)
  
  
  
  narows=which(apply(data,1,function(x){any(is.na(x))}))
  
  navars=which(apply(data,2,function(x){any(is.na(x))}))
  
  
  
  colsample=sample(navars)
  
  varcombs=unique(apply(data,1,function(y){which(is.na(y))}))
  
  varcombs=varcombs[lengths(varcombs)>0]
  
  varcombstemp=varcombs
  
  for (i in colsample) {
    
    if (length(varcombstemp)==0) {
      
      break
      
    }
    
    varcombstemp=unique(lapply(varcombstemp,function(x){if (any(x==i)) {return(x[-which(x==i)])} else {return(x)}}))
    
    varcombstemp=varcombstemp[lengths(varcombstemp)>0]
    
    varcombs=unique(c(varcombs,varcombstemp))
    
  }
  
  
  
  helplist=list(matrices=list(rules),pattern=list(colnames(rules)[-length(colnames(rules))]))
  
  for (i in 1:length(varcombs)) {
    
    varcombstemp=varcombs[[i]]
    
    elim=helplist$pattern[[1]]
    
    elim=elim[which(!elim%in%c(colnames(data)[-varcombstemp],errnames))]
    
    M=helplist$matrices[[1]]
    
    for (a in 1:length(elim)) {
      
      M=eliminate(M,elim[a])
      
    }
    
    helplist$matrices[[i+1]]=M
    
    helplist$pattern[[i+1]]=c(colnames(data)[-varcombstemp],errnames)
    
  }
  
  
  
  inits=data-data
  
  inits[is.na(inits)]=0
  
  inits2=data
  
  for (i in colnames(data)[colsample]) {
    
    nas=which(is.na(inits2[,i]))
    
    bounds=t(vapply(nas,function(x){solvedmat(crow=x,x1=inits2[,!is.na(inits2[x,]),drop=FALSE],cname1=i,cnames1=colnames(data)[!is.na(inits2[x,])],tolfac=0)},FUN.VALUE= numeric(2)))
    
    bounds=bounds/wts[nas]
    
    inits2[nas,i]=sapply(1:length(nas),function(x){runif(1,bounds[x,1],bounds[x,2])})
    
  }
  
  inits=inits2
  
}
