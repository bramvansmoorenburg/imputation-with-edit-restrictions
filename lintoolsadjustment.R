

lintadj=function(imputed,data,totals,edits,Weights=NULL,maxits=5,maxits2=1e3) {
  
  
  
  if (any(is.na(imputed))) {
    
    stop("(lintadj) NA's in imputed dataset")
    
  } else if (any(abs(imputed)==Inf)) {
    
    stop("(lintadj) infinite values in imputed dataset")
    
  }
  
  
  
  if (!is.null(Weights)) {
    
    wts=Weights
    
  } else {
    
    wts=rep(1,nrow(imputed))
    
  }
  
  meanwts=wts/mean(wts)
  
  
  
  errnames=colnames(edits)[!colnames(edits)%in%c("CONSTANT",colnames(data))]
  
  if (length(errnames)!=0) {
    
    rules=reduce(substValue(edits,errnames,rep(0,length(errnames))))
    
  } else {
    
    rules=edits
    
  }
  
  rules=rules[,c(colnames(data),"CONSTANT")]
  
  
  
  totmiss=totals-colSums(data*wts,na.rm=TRUE)
  
  
  
  out=imputed
  
  narows=which(apply(data,1,function(x){any(is.na(x))}))
  
  navars=which(apply(data,2,function(x){any(is.na(x))})&!is.na(totals))
  
  totchange=1
  
  it=1
  
  while (it<=maxits&&totchange>0) {
    
    for (j in navars) {
      
      if (!is.na(totals[j])) {
        
        nas=which(is.na(data[,j]))
        
        if (all(out[nas,j]>=0)||all(out[nas,j]<=0)) {
          
          sumout=ifelse(sum(out[nas,j]*wts[nas])!=0,sum(out[nas,j]*wts[nas]),1)
          
          out[nas,j]=out[nas,j]*totmiss[j]/sumout
          
        } else {
          
          tmpchange=(totmiss[j]-sum(out[nas,j]*wts[nas]))
          
          sumout=ifelse(sum(out[nas,j][out[nas,j]<0]*wts[nas][out[nas,j]<0])!=0,sum(out[nas,j][out[nas,j]<0]*wts[nas][out[nas,j]<0]),1)
          
          tmpratio=-sum(out[nas,j][out[nas,j]<0]*wts[nas][out[nas,j]<0])/sum(abs(out[nas,j])*wts[nas])
          
          out[nas,j][out[nas,j]<0]=out[nas,j][out[nas,j]<0]+out[nas,j][out[nas,j]<0]*tmpchange*tmpratio/sumout
          
          out[nas,j][out[nas,j]>0]=out[nas,j][out[nas,j]>0]+out[nas,j][out[nas,j]>0]*tmpchange*(1-tmpratio)/sumout
          
        }
        
      }
      
    }
    
    tempchange=numeric(0)
    
    for (i in narows) {
      
      editwts=abs(1/as.numeric(out[i,]))
      
      editwts[editwts==Inf]=2
      
      tmprules=substValue(rules,colnames(data)[!is.na(data[i,])],data[i,!is.na(data[i,])])
      
      tempadj=project(as.numeric(out[i,]),getA(tmprules),getb(tmprules),0,w=editwts,eps=.Machine$double.eps,maxiter = max(10,round(maxits2*meanwts[i])))
      
      tempchange=c(tempchange,tempadj$objective)
      
      out[i,]=tempadj$x
      
    }
    
    totchange=sum(tempchange)
    
    it=it+1
    
  }
  
  
  
  return(out)
  
}
