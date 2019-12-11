

missmech=function(data) {
  
  missind=apply(is.na(data),1,any)
  
  varind=apply(data,2,function(x){any(!is.na(x)&missind)})
  
  varmods=lapply(1:ncol(data),function(x){if (varind[x]) {return(glm(missind~data[,x],family = binomial))} else {return(NA)}})
  
  varsummary=lapply(1:ncol(data),function(x){if (varind[x]) {return(summary(varmods[[x]])$coefficients)} else {return(NA)}})
  
  varout=sapply(1:ncol(data),function(x){if (varind[x]) {return(varsummary[[x]][-1,1][which.max(abs(varsummary[[x]][-1,1]))])} else {return(NA)}})
  
  varout=rbind(varout,sapply(1:ncol(data),function(x){if (varind[x]) {return(min(varsummary[[x]][-1,4]))} else {return(NA)}}))
  
  rownames(varout)=c("coefficient","p.value")
  
  colnames(varout)=colnames(data)
  
  return(varout)
  
}
