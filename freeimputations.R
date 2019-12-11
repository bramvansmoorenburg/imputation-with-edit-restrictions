

freeimp=function(data,rules,totals,Weight=NULL,maxerr=1e-5,nonnegative=FALSE) {
  
  
  
  if (!is.null(Weight)) {
    
    wts=Weight
    
  } else {
    
    wts=rep(1,nrow(data))
    
  }
  
  
  
  errnames=colnames(rules)[!colnames(rules)%in%c("CONSTANT",colnames(data))]
  
  if (length(errnames)==0) {
    
    stop("(freeimp) edits have no error terms")
    
  }
  
  
  
  totmiss=totals-colSums(data*wts,na.rm=TRUE)
  
  
  
  solvedmat=function(crow,x1,cname1,cnames1) {
    
    edmat=sapply(helplist$pattern,function(b){ifelse(all(c(cnames1,cname1,errnames)%in%b)&all(b%in%c(cnames1,cname1,errnames)),TRUE,FALSE)})
    
    edmat=which(edmat)
    
    edmatr=helplist$matrices[[edmat]]
    
    
    
    edmatr=substValue(edmatr,errnames,rep(0,length(errnames)))
    
    
    
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
  
  
  
  nahist=sum(is.na(data))
  
  data2=data
  
  repeat {
    
    navars=which(apply(data2,2,function(x){any(is.na(x))}))
    
    for (i in colnames(data)[navars]) {
      
      nas=which(is.na(data2[,i]))
      
      helplist=list(matrices=list(rules),pattern=list(colnames(rules)[-length(colnames(rules))]))
      
      varcombs=lapply(nas,function(x){which(is.na(data2[x,]))})
      
      varcombs=unique(lapply(varcombs,function(x){x[-which(x==which(colnames(data)==i))]}))
      
      for (a in 1:length(varcombs)) {
        
        varcombstemp=as.numeric(varcombs[[a]])
        
        elim=helplist$pattern[[1]]
        
        elim=elim[which(!elim%in%c(colnames(data)[-varcombstemp],errnames))]
        
        M=helplist$matrices[[1]]
        
        for (b in 1:length(elim)) {
          
          M=try(eliminate(M,elim[b]))
          
          if (class(M)=="try-error") {stop("Function 'eliminate' did not work. Be sure to use the function eliminate from package editrules.")}
          
        }
        
        helplist$matrices[[a+1]]=M
        
        helplist$pattern[[a+1]]=c(colnames(data)[-varcombstemp],errnames)
        
      }
      
      bounds=t(vapply(nas,function(x){solvedmat(crow=x,x1=data2[,!is.na(data2[x,]),drop=FALSE],cname1=i,cnames1=colnames(data)[!is.na(data2[x,])])},FUN.VALUE = numeric(2)))
      
      boundsequal=abs(bounds[,1]-bounds[,2])<maxerr
      
      data2[nas[boundsequal],i]=bounds[boundsequal,1]/wts[nas[boundsequal]]
      
    }
    
    totmiss2=totals-colSums(data2*wts,na.rm=TRUE)
    
    if (nonnegative&&!is.na(any(na.omit(totmiss2)==0))&&any(na.omit(totmiss2)==0)) {
      
      for (i in which(totmiss2==0)) {
        
        data2[is.na(data2[,i]),i]=0
        
      }
      
    }
    
    if (any(apply(data2,2,function(s){sum(is.na(s))})==1&!is.na(totals))) {
      
      for (i in which(apply(data2,2,function(s){sum(is.na(s))})==1&!is.na(totals))) {
        
        data2[is.na(data2[,i]),i]=totmiss2[i]/wts[is.na(data2[,i])]
        
      }
      
    }
    
    if (sum(is.na(data2))==nahist) {
      
      break
      
    } else {
      
      nahist=sum(is.na(data2))
      
    }
    
  }
  
  
  
  return(data2)
  
}
