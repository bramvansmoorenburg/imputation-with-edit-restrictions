
starttime=Sys.time()



library(lintools)

library(editrules)

library(mice)

library(Hmisc)

library(parallel)

library(mvtnorm)



response=2
explanatory=c(5,6,8)

truevals=read.csv("csv datRcom.csv")
original=truevals[,-c(1:3)]
orignames=colnames(original)
original=original[,-3]
mu=colMeans(original)
Sigma=cov(original)
Sigma2=Sigma/10
diag(Sigma2)=diag(Sigma2)*10
wts=rep(1,nrow(truevals))
truevals=truevals[,-c(1:3)]
totals=colSums(truevals*wts)



parallelimp=6 #detectCores()-1

rubindatasets=20

rubindatasets2=10

rubindatasets3=5

holebatches=5

datbatches=5

maxhours=Inf

frac1=0.05

frac2=0.1

frac3=0.2

maxits1=5

maxits2=10

maxits3=20

vars=1:ncol(truevals)

n=nrow(truevals)

startseed=round(as.numeric(starttime))

logtransform=TRUE

givefeedback=TRUE



source("freeimputations.R")

source("sampleimp.R")

source("lintoolsadjustment.R")

source("imputationvalidator.R")



rules=editmatrix(as.character(unlist(read.table("test edits (errors).txt"))))

errnames=colnames(rules)[!colnames(rules)%in%c("CONSTANT",orignames)]

rulesraw=reduce(substValue(rules,errnames,rep(0,length(errnames))))



basic=function(x,freeimp=freeimp,lintadj=lintadj,data=data,rules=rules,totals=totals,wts=wts,impmeth=impmeth) {
  
  library(lintools)
  
  library(editrules)
  
  library(mice)
  
  data2=freeimp(data,rules,totals,wts)
  
  impdata=complete(mice(data2,1,impmeth,threshold=0.999999999999999,print=FALSE))
  
  out=lintadj(impdata,data2,totals,rules,wts)
  
  return(out)
  
}

advanced=function(x,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=data,rules=rules,totals=totals,wts=wts,impmeth=impmeth,maxits=20) {
  
  library(lintools)
  
  library(editrules)
  
  library(mice)
  
  data2=freeimp(data,rules,totals,wts)
  
  impdata=sampleimp(data2,function(y){complete(suppressWarnings(mice(y,1,impmeth,threshold=0.999999999999999,print=FALSE)))},rules,totals,weight=wts,maxit=maxits)
  
  out=lintadj(impdata,data2,totals,rules,wts)
  
  return(out)
  
}

impmeth="pmm"



holecycle=ceiling(rubindatasets/parallelimp)*holebatches

it=0

outlist=outlistcont=outlistsim1=outlistsim2=outlistsimcont1=outlistsimcont2=NULL

imputedsim111=imputedsim112=imputedsim113=imputedsim121=imputedsim122=imputedsim123=imputedsim131=imputedsim132=imputedsim133=imputedsim211=imputedsim212=imputedsim213=imputedsim221=imputedsim222=imputedsim223=imputedsim231=imputedsim232=imputedsim233=list()

imputedsimcont11=imputedsimcont12=imputedsimcont13=imputedsimcont21=imputedsimcont22=imputedsimcont23=list()

imputed11=imputed12=imputed13=imputed21=imputed22=imputed23=imputed31=imputed32=imputed33=list()

imputedcont1=imputedcont2=imputedcont3=list()

rubinlist=rubinlistcont=rubinlistmice=rubinlistsim1=rubinlistsim2=rubinlistsimcont1=rubinlistsimcont2=rubinlistmicesim1=rubinlistmicesim2=list()



seeddat=seedhole=startseed

clu=makeCluster(parallelimp)

repeat {
  
  it=it+parallelimp
  
  if (givefeedback) {
    
    cat("\n","Iteration: ",it,", ",sep="")
    
  }
  
  if (((it/parallelimp)-1)%%(holecycle*datbatches)==0) { #making new datasets
    
    set.seed(seeddat)
    
    seeddat=.Random.seed[3]
    
    set.seed(seeddat)
    
    truevalssim1=truevalssim2=matrix(0,nrow=0,ncol = length(mu))
    
    colnames(truevalssim1)=names(mu)
    
    colnames(truevalssim2)=names(mu)
    
    truevalssim1=as.data.frame(truevalssim1)
    
    truevalssim2=as.data.frame(truevalssim2)
    
    while (nrow(truevalssim1)<n) { #making dataset 1
      
      tempvals=as.data.frame(rmvnorm(n,mu,Sigma))
      
      colnames(tempvals)=names(mu)
      
      tempvals$PURESALE=tempvals$PURTOT-tempvals$PUROTHAL
      
      violations=violatedEdits(rulesraw,tempvals)
      
      tempvals=tempvals[which(!apply(violations,1,any)),]
      
      truevalssim1=rbind(truevalssim1,tempvals)
      
    }
    
    while (nrow(truevalssim2)<n) { #making dataset 2
      
      tempvals=as.data.frame(rmvnorm(n,mu,Sigma2))
      
      colnames(tempvals)=names(mu)
      
      tempvals$PURESALE=tempvals$PURTOT-tempvals$PUROTHAL
      
      violations=violatedEdits(rulesraw,tempvals)
      
      tempvals=tempvals[which(!apply(violations,1,any)),]
      
      truevalssim2=rbind(truevalssim2,tempvals)
      
    }
    
    truevalssim1=truevalssim1[1:n,colnames(truevals)]
    
    truevalssim2=truevalssim2[1:n,colnames(truevals)]
    
    wtssim1=rep(1,nrow(truevalssim1))
    
    wtssim2=rep(1,nrow(truevalssim2))
    
    totalssim1=colSums(truevalssim1*wtssim1)
    
    totalssim2=colSums(truevalssim2*wtssim2)
    
    if (givefeedback) {
      
      cat("made new datasets, ")
      
    }
    
  }
  
  if (((it/parallelimp)-1)%%holecycle==0) {
    
    set.seed(seedhole)
    
    seedhole=.Random.seed[4]
    
    set.seed(seedhole)
    
    holessim11=holessim12=holessim13=as.matrix(truevalssim1[,vars])
    
    holessim21=holessim22=holessim23=as.matrix(truevalssim2[,vars])
    
    nassim11=sample(prod(dim(holessim11)),frac1*prod(dim(holessim11)))
    
    nassim12=sample(prod(dim(holessim12)),frac2*prod(dim(holessim12)))
    
    nassim13=sample(prod(dim(holessim13)),frac3*prod(dim(holessim13)))
    
    nassim21=sample(prod(dim(holessim21)),frac1*prod(dim(holessim21)))
    
    nassim22=sample(prod(dim(holessim22)),frac2*prod(dim(holessim22)))
    
    nassim23=sample(prod(dim(holessim23)),frac3*prod(dim(holessim23)))
    
    holessim11[nassim11]=NA
    
    holessim12[nassim12]=NA
    
    holessim13[nassim13]=NA
    
    holessim21[nassim21]=NA
    
    holessim22[nassim22]=NA
    
    holessim23[nassim23]=NA
    
    holessim11=as.data.frame(holessim11)
    
    holessim12=as.data.frame(holessim12)
    
    holessim13=as.data.frame(holessim13)
    
    holessim21=as.data.frame(holessim21)
    
    holessim22=as.data.frame(holessim22)
    
    holessim23=as.data.frame(holessim23)
    
    datasim11=datasim12=datasim13=as.data.frame(truevalssim1)
    
    datasim21=datasim22=datasim23=as.data.frame(truevalssim2)
    
    datasim11[,vars]=holessim11
    
    datasim12[,vars]=holessim12
    
    datasim13[,vars]=holessim13
    
    datasim21[,vars]=holessim21
    
    datasim22[,vars]=holessim22
    
    datasim23[,vars]=holessim23
    
    datafreesim11=freeimp(datasim11,rules,totalssim1,wtssim1)
    
    datafreesim12=freeimp(datasim12,rules,totalssim1,wtssim1)
    
    datafreesim13=freeimp(datasim13,rules,totalssim1,wtssim1)
    
    datafreesim21=freeimp(datasim21,rules,totalssim2,wtssim2)
    
    datafreesim22=freeimp(datasim22,rules,totalssim2,wtssim2)
    
    datafreesim23=freeimp(datasim23,rules,totalssim2,wtssim2)
    
    holes1=holes2=holes3=as.matrix(truevals[,vars])
    
    nas1=sample(prod(dim(holes1)),frac1*prod(dim(holes1)))
    
    nas2=sample(prod(dim(holes2)),frac2*prod(dim(holes2)))
    
    nas3=sample(prod(dim(holes3)),frac3*prod(dim(holes3)))
    
    holes1[nas1]=NA
    
    holes2[nas2]=NA
    
    holes3[nas3]=NA
    
    holes1=as.data.frame(holes1)
    
    holes2=as.data.frame(holes2)
    
    holes3=as.data.frame(holes3)
    
    data1=data2=data3=as.data.frame(truevals)
    
    data1[,vars]=holes1
    
    data2[,vars]=holes2
    
    data3[,vars]=holes3
    
    datafree1=freeimp(data1,rules,totals,wts)
    
    datafree2=freeimp(data2,rules,totals,wts)
    
    datafree3=freeimp(data3,rules,totals,wts)
    
    imputedsim111=imputedsim112=imputedsim113=imputedsim121=imputedsim122=imputedsim123=imputedsim131=imputedsim132=imputedsim133=imputedsim211=imputedsim212=imputedsim213=imputedsim221=imputedsim222=imputedsim223=imputedsim231=imputedsim232=imputedsim233=list()
    
    imputedsimcont11=imputedsimcont12=imputedsimcont13=imputedsimcont21=imputedsimcont22=imputedsimcont23=list()
    
    imputed11=imputed12=imputed13=imputed21=imputed22=imputed23=imputed31=imputed32=imputed33=list()
    
    imputedcont1=imputedcont2=imputedcont3=list()
    
    if (givefeedback) {
      
      cat("made new holes, ")
      
    }
    
  }
  
  
  
  #imputed simulations algorithm
  
  imputedsim111=c(imputedsim111,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim11,rules=rules,totals=totalssim1,wts=wtssim1,impmeth=impmeth,maxits=maxits1))
  
  imputedsim112=c(imputedsim112,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim11,rules=rules,totals=totalssim1,wts=wtssim1,impmeth=impmeth,maxits=maxits2))
  
  imputedsim113=c(imputedsim113,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim11,rules=rules,totals=totalssim1,wts=wtssim1,impmeth=impmeth,maxits=maxits3))
  
  imputedsim121=c(imputedsim121,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim12,rules=rules,totals=totalssim1,wts=wtssim1,impmeth=impmeth,maxits=maxits1))
  
  imputedsim122=c(imputedsim122,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim12,rules=rules,totals=totalssim1,wts=wtssim1,impmeth=impmeth,maxits=maxits2))
  
  imputedsim123=c(imputedsim123,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim12,rules=rules,totals=totalssim1,wts=wtssim1,impmeth=impmeth,maxits=maxits3))
  
  imputedsim131=c(imputedsim131,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim13,rules=rules,totals=totalssim1,wts=wtssim1,impmeth=impmeth,maxits=maxits1))
  
  imputedsim132=c(imputedsim132,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim13,rules=rules,totals=totalssim1,wts=wtssim1,impmeth=impmeth,maxits=maxits2))
  
  imputedsim133=c(imputedsim133,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim13,rules=rules,totals=totalssim1,wts=wtssim1,impmeth=impmeth,maxits=maxits3))
  
  imputedsim211=c(imputedsim211,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim21,rules=rules,totals=totalssim2,wts=wtssim2,impmeth=impmeth,maxits=maxits1))
  
  imputedsim212=c(imputedsim212,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim21,rules=rules,totals=totalssim2,wts=wtssim2,impmeth=impmeth,maxits=maxits2))
  
  imputedsim213=c(imputedsim213,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim21,rules=rules,totals=totalssim2,wts=wtssim2,impmeth=impmeth,maxits=maxits3))
  
  imputedsim221=c(imputedsim221,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim22,rules=rules,totals=totalssim2,wts=wtssim2,impmeth=impmeth,maxits=maxits1))
  
  imputedsim222=c(imputedsim222,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim22,rules=rules,totals=totalssim2,wts=wtssim2,impmeth=impmeth,maxits=maxits2))
  
  imputedsim223=c(imputedsim223,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim22,rules=rules,totals=totalssim2,wts=wtssim2,impmeth=impmeth,maxits=maxits3))
  
  imputedsim231=c(imputedsim231,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim23,rules=rules,totals=totalssim2,wts=wtssim2,impmeth=impmeth,maxits=maxits1))
  
  imputedsim232=c(imputedsim232,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim23,rules=rules,totals=totalssim2,wts=wtssim2,impmeth=impmeth,maxits=maxits2))
  
  imputedsim233=c(imputedsim233,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=datasim23,rules=rules,totals=totalssim2,wts=wtssim2,impmeth=impmeth,maxits=maxits3))
  
  
  
  #imputed simulations control
  
  imputedsimcont11=c(imputedsimcont11,parLapply(clu,1:parallelimp,basic,freeimp=freeimp,lintadj=lintadj,data=datasim11,rules=rules,totals=totalssim1,wts=wtssim1,impmeth=impmeth))
  
  imputedsimcont12=c(imputedsimcont12,parLapply(clu,1:parallelimp,basic,freeimp=freeimp,lintadj=lintadj,data=datasim12,rules=rules,totals=totalssim1,wts=wtssim1,impmeth=impmeth))
  
  imputedsimcont13=c(imputedsimcont13,parLapply(clu,1:parallelimp,basic,freeimp=freeimp,lintadj=lintadj,data=datasim13,rules=rules,totals=totalssim1,wts=wtssim1,impmeth=impmeth))
  
  imputedsimcont21=c(imputedsimcont21,parLapply(clu,1:parallelimp,basic,freeimp=freeimp,lintadj=lintadj,data=datasim21,rules=rules,totals=totalssim2,wts=wtssim2,impmeth=impmeth))
  
  imputedsimcont22=c(imputedsimcont22,parLapply(clu,1:parallelimp,basic,freeimp=freeimp,lintadj=lintadj,data=datasim22,rules=rules,totals=totalssim2,wts=wtssim2,impmeth=impmeth))
  
  imputedsimcont23=c(imputedsimcont23,parLapply(clu,1:parallelimp,basic,freeimp=freeimp,lintadj=lintadj,data=datasim23,rules=rules,totals=totalssim2,wts=wtssim2,impmeth=impmeth))
  
  
  
  #imputed real algorithm
  
  imputed11=c(imputed11,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=data1,rules=rules,totals=totals,wts=wts,impmeth=impmeth,maxits=maxits1))
  
  imputed12=c(imputed12,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=data1,rules=rules,totals=totals,wts=wts,impmeth=impmeth,maxits=maxits2))
  
  imputed13=c(imputed13,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=data1,rules=rules,totals=totals,wts=wts,impmeth=impmeth,maxits=maxits3))
  
  imputed21=c(imputed21,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=data2,rules=rules,totals=totals,wts=wts,impmeth=impmeth,maxits=maxits1))
  
  imputed22=c(imputed22,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=data2,rules=rules,totals=totals,wts=wts,impmeth=impmeth,maxits=maxits2))
  
  imputed23=c(imputed23,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=data2,rules=rules,totals=totals,wts=wts,impmeth=impmeth,maxits=maxits3))
  
  imputed31=c(imputed31,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=data3,rules=rules,totals=totals,wts=wts,impmeth=impmeth,maxits=maxits1))
  
  imputed32=c(imputed32,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=data3,rules=rules,totals=totals,wts=wts,impmeth=impmeth,maxits=maxits2))
  
  imputed33=c(imputed33,parLapply(clu,1:parallelimp,advanced,freeimp=freeimp,sampleimp=sampleimp,lintadj=lintadj,data=data3,rules=rules,totals=totals,wts=wts,impmeth=impmeth,maxits=maxits3))
  
  
  
  #imputed real control
  
  imputedcont1=c(imputedcont1,parLapply(clu,1:parallelimp,basic,freeimp=freeimp,lintadj=lintadj,data=data1,rules=rules,totals=totals,wts=wts,impmeth=impmeth))
  
  imputedcont2=c(imputedcont2,parLapply(clu,1:parallelimp,basic,freeimp=freeimp,lintadj=lintadj,data=data2,rules=rules,totals=totals,wts=wts,impmeth=impmeth))
  
  imputedcont3=c(imputedcont3,parLapply(clu,1:parallelimp,basic,freeimp=freeimp,lintadj=lintadj,data=data3,rules=rules,totals=totals,wts=wts,impmeth=impmeth))
  
  
  
  if (givefeedback) {
    
    cat("validating, ")
    
  }
  
  #validated simulations algorithm
  
  validatedsim111=parLapply(clu,imputedsim111[(length(imputedsim111)-parallelimp+1):length(imputedsim111)],impvalidator,dataset=datasim11,truevals=truevalssim1,rules=rules,weight=wtssim1,frac=frac1,maxits=maxits1)
  
  validatedsim112=parLapply(clu,imputedsim112[(length(imputedsim112)-parallelimp+1):length(imputedsim112)],impvalidator,dataset=datasim11,truevals=truevalssim1,rules=rules,weight=wtssim1,frac=frac1,maxits=maxits2)
  
  validatedsim113=parLapply(clu,imputedsim113[(length(imputedsim113)-parallelimp+1):length(imputedsim113)],impvalidator,dataset=datasim11,truevals=truevalssim1,rules=rules,weight=wtssim1,frac=frac1,maxits=maxits3)
  
  validatedsim121=parLapply(clu,imputedsim121[(length(imputedsim121)-parallelimp+1):length(imputedsim121)],impvalidator,dataset=datasim12,truevals=truevalssim1,rules=rules,weight=wtssim1,frac=frac2,maxits=maxits1)
  
  validatedsim122=parLapply(clu,imputedsim122[(length(imputedsim122)-parallelimp+1):length(imputedsim122)],impvalidator,dataset=datasim12,truevals=truevalssim1,rules=rules,weight=wtssim1,frac=frac2,maxits=maxits2)
  
  validatedsim123=parLapply(clu,imputedsim123[(length(imputedsim123)-parallelimp+1):length(imputedsim123)],impvalidator,dataset=datasim12,truevals=truevalssim1,rules=rules,weight=wtssim1,frac=frac2,maxits=maxits3)
  
  validatedsim131=parLapply(clu,imputedsim131[(length(imputedsim131)-parallelimp+1):length(imputedsim131)],impvalidator,dataset=datasim13,truevals=truevalssim1,rules=rules,weight=wtssim1,frac=frac3,maxits=maxits1)
  
  validatedsim132=parLapply(clu,imputedsim132[(length(imputedsim132)-parallelimp+1):length(imputedsim132)],impvalidator,dataset=datasim13,truevals=truevalssim1,rules=rules,weight=wtssim1,frac=frac3,maxits=maxits2)
  
  validatedsim133=parLapply(clu,imputedsim133[(length(imputedsim133)-parallelimp+1):length(imputedsim133)],impvalidator,dataset=datasim13,truevals=truevalssim1,rules=rules,weight=wtssim1,frac=frac3,maxits=maxits3)
  
  validatedsim211=parLapply(clu,imputedsim211[(length(imputedsim111)-parallelimp+1):length(imputedsim211)],impvalidator,dataset=datasim21,truevals=truevalssim2,rules=rules,weight=wtssim2,frac=frac1,maxits=maxits1)
  
  validatedsim212=parLapply(clu,imputedsim212[(length(imputedsim212)-parallelimp+1):length(imputedsim212)],impvalidator,dataset=datasim21,truevals=truevalssim2,rules=rules,weight=wtssim2,frac=frac1,maxits=maxits2)
  
  validatedsim213=parLapply(clu,imputedsim213[(length(imputedsim213)-parallelimp+1):length(imputedsim213)],impvalidator,dataset=datasim21,truevals=truevalssim2,rules=rules,weight=wtssim2,frac=frac1,maxits=maxits3)
  
  validatedsim221=parLapply(clu,imputedsim221[(length(imputedsim221)-parallelimp+1):length(imputedsim221)],impvalidator,dataset=datasim22,truevals=truevalssim2,rules=rules,weight=wtssim2,frac=frac2,maxits=maxits1)
  
  validatedsim222=parLapply(clu,imputedsim222[(length(imputedsim222)-parallelimp+1):length(imputedsim222)],impvalidator,dataset=datasim22,truevals=truevalssim2,rules=rules,weight=wtssim2,frac=frac2,maxits=maxits2)
  
  validatedsim223=parLapply(clu,imputedsim223[(length(imputedsim223)-parallelimp+1):length(imputedsim223)],impvalidator,dataset=datasim22,truevals=truevalssim2,rules=rules,weight=wtssim2,frac=frac2,maxits=maxits3)
  
  validatedsim231=parLapply(clu,imputedsim231[(length(imputedsim231)-parallelimp+1):length(imputedsim231)],impvalidator,dataset=datasim23,truevals=truevalssim2,rules=rules,weight=wtssim2,frac=frac3,maxits=maxits1)
  
  validatedsim232=parLapply(clu,imputedsim232[(length(imputedsim232)-parallelimp+1):length(imputedsim232)],impvalidator,dataset=datasim23,truevals=truevalssim2,rules=rules,weight=wtssim2,frac=frac3,maxits=maxits2)
  
  validatedsim233=parLapply(clu,imputedsim233[(length(imputedsim233)-parallelimp+1):length(imputedsim233)],impvalidator,dataset=datasim23,truevals=truevalssim2,rules=rules,weight=wtssim2,frac=frac3,maxits=maxits3)
  
  
  
  #validated simulations control
  
  validatedsimcont11=parLapply(clu,imputedsimcont11[(length(imputedsimcont11)-parallelimp+1):length(imputedsimcont11)],impvalidator,dataset=datasim11,truevals=truevalssim1,rules=rules,weight=wtssim1,frac=frac1,maxits=1)
  
  validatedsimcont12=parLapply(clu,imputedsimcont12[(length(imputedsimcont12)-parallelimp+1):length(imputedsimcont12)],impvalidator,dataset=datasim12,truevals=truevalssim1,rules=rules,weight=wtssim1,frac=frac2,maxits=1)
  
  validatedsimcont13=parLapply(clu,imputedsimcont13[(length(imputedsimcont13)-parallelimp+1):length(imputedsimcont13)],impvalidator,dataset=datasim13,truevals=truevalssim1,rules=rules,weight=wtssim1,frac=frac3,maxits=1)
  
  validatedsimcont21=parLapply(clu,imputedsimcont21[(length(imputedsimcont21)-parallelimp+1):length(imputedsimcont21)],impvalidator,dataset=datasim21,truevals=truevalssim2,rules=rules,weight=wtssim2,frac=frac1,maxits=1)
  
  validatedsimcont22=parLapply(clu,imputedsimcont22[(length(imputedsimcont22)-parallelimp+1):length(imputedsimcont22)],impvalidator,dataset=datasim22,truevals=truevalssim2,rules=rules,weight=wtssim2,frac=frac2,maxits=1)
  
  validatedsimcont23=parLapply(clu,imputedsimcont23[(length(imputedsimcont23)-parallelimp+1):length(imputedsimcont23)],impvalidator,dataset=datasim23,truevals=truevalssim2,rules=rules,weight=wtssim2,frac=frac3,maxits=1)
  
  
  
  #Validated real algorithm
  
  validated11=parLapply(clu,imputed11[(length(imputed11)-parallelimp+1):length(imputed11)],impvalidator,dataset=data1,truevals=truevals,rules=rules,weight=wts,frac=frac1,maxits=maxits1)
  
  validated12=parLapply(clu,imputed12[(length(imputed12)-parallelimp+1):length(imputed12)],impvalidator,dataset=data1,truevals=truevals,rules=rules,weight=wts,frac=frac1,maxits=maxits2)
  
  validated13=parLapply(clu,imputed13[(length(imputed13)-parallelimp+1):length(imputed13)],impvalidator,dataset=data1,truevals=truevals,rules=rules,weight=wts,frac=frac1,maxits=maxits3)
  
  validated21=parLapply(clu,imputed21[(length(imputed21)-parallelimp+1):length(imputed11)],impvalidator,dataset=data2,truevals=truevals,rules=rules,weight=wts,frac=frac2,maxits=maxits1)
  
  validated22=parLapply(clu,imputed22[(length(imputed22)-parallelimp+1):length(imputed22)],impvalidator,dataset=data2,truevals=truevals,rules=rules,weight=wts,frac=frac2,maxits=maxits2)
  
  validated23=parLapply(clu,imputed23[(length(imputed23)-parallelimp+1):length(imputed23)],impvalidator,dataset=data2,truevals=truevals,rules=rules,weight=wts,frac=frac2,maxits=maxits3)
  
  validated31=parLapply(clu,imputed31[(length(imputed31)-parallelimp+1):length(imputed31)],impvalidator,dataset=data3,truevals=truevals,rules=rules,weight=wts,frac=frac3,maxits=maxits1)
  
  validated32=parLapply(clu,imputed32[(length(imputed32)-parallelimp+1):length(imputed32)],impvalidator,dataset=data3,truevals=truevals,rules=rules,weight=wts,frac=frac3,maxits=maxits2)
  
  validated33=parLapply(clu,imputed33[(length(imputed33)-parallelimp+1):length(imputed33)],impvalidator,dataset=data3,truevals=truevals,rules=rules,weight=wts,frac=frac3,maxits=maxits3)
  
  
  
  #validated real control
  
  validatedcont1=parLapply(clu,imputedcont1[(length(imputedcont1)-parallelimp+1):length(imputedcont1)],impvalidator,dataset=data1,truevals=truevals,rules=rules,weight=wts,frac=frac1,maxits=1)
  
  validatedcont2=parLapply(clu,imputedcont2[(length(imputedcont2)-parallelimp+1):length(imputedcont2)],impvalidator,dataset=data2,truevals=truevals,rules=rules,weight=wts,frac=frac2,maxits=1)
  
  validatedcont3=parLapply(clu,imputedcont3[(length(imputedcont1)-parallelimp+1):length(imputedcont3)],impvalidator,dataset=data3,truevals=truevals,rules=rules,weight=wts,frac=frac3,maxits=1)
  
  
  
  if (length(outlistsim1)!=0) {
    
    outlistsim1=list(outlistsim1)
    
  }
  
  if (length(outlistsimcont1)!=0) {
    
    outlistsimcont1=list(outlistsimcont1)
    
  }
  
  if (length(outlistsim2)!=0) {
    
    outlistsim2=list(outlistsim2)
    
  }
  
  if (length(outlistsimcont2)!=0) {
    
    outlistsimcont2=list(outlistsimcont2)
    
  }
  
  if (length(outlist)!=0) {
    
    outlist=list(outlist)
    
  }
  
  if (length(outlistcont)!=0) {
    
    outlistcont=list(outlistcont)
    
  }
  
  outlistsim1=validatecleaner(c(outlistsim1,validatedsim111,validatedsim112,validatedsim113,validatedsim121,validatedsim122,validatedsim123,validatedsim131,validatedsim132,validatedsim133))
  
  outlistsim2=validatecleaner(c(outlistsim2,validatedsim211,validatedsim212,validatedsim213,validatedsim221,validatedsim222,validatedsim223,validatedsim231,validatedsim232,validatedsim233))
  
  outlistsimcont1=validatecleaner(c(outlistsimcont1,validatedsimcont11,validatedsimcont12,validatedsimcont13))
  
  outlistsimcont2=validatecleaner(c(outlistsimcont2,validatedsimcont21,validatedsimcont22,validatedsimcont23))
  
  outlist=validatecleaner(c(outlist,validated11,validated12,validated13,validated21,validated22,validated23,validated31,validated32,validated33))
  
  outlistcont=validatecleaner(c(outlistcont,validatedcont1,validatedcont2,validatedcont3))
  
  
  
  if (length(imputedsim111)>=rubindatasets) {
    
    if (givefeedback) {
      
      cat("rubin-validating")
      
    }
    
    for (xa in 1:(length(imputedsim111)%/%rubindatasets)) {
      
      
      
      #rubin-validating simulations algorithm
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim111[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim11,truevalssim1,response,explanatory,frac=frac1,maxits=maxits1,datinpool=rubindatasets)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim112[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim11,truevalssim1,response,explanatory,frac=frac1,maxits=maxits2,datinpool=rubindatasets)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim113[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim11,truevalssim1,response,explanatory,frac=frac1,maxits=maxits3,datinpool=rubindatasets)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim121[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim12,truevalssim1,response,explanatory,frac=frac2,maxits=maxits1,datinpool=rubindatasets)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim122[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim12,truevalssim1,response,explanatory,frac=frac2,maxits=maxits2,datinpool=rubindatasets)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim123[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim12,truevalssim1,response,explanatory,frac=frac2,maxits=maxits3,datinpool=rubindatasets)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim131[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim13,truevalssim1,response,explanatory,frac=frac3,maxits=maxits1,datinpool=rubindatasets)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim132[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim13,truevalssim1,response,explanatory,frac=frac3,maxits=maxits2,datinpool=rubindatasets)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim133[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim13,truevalssim1,response,explanatory,frac=frac3,maxits=maxits3,datinpool=rubindatasets)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim211[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim21,truevalssim2,response,explanatory,frac=frac1,maxits=maxits1,datinpool=rubindatasets)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim212[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim21,truevalssim2,response,explanatory,frac=frac1,maxits=maxits2,datinpool=rubindatasets)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim213[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim21,truevalssim2,response,explanatory,frac=frac1,maxits=maxits3,datinpool=rubindatasets)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim221[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim22,truevalssim2,response,explanatory,frac=frac2,maxits=maxits1,datinpool=rubindatasets)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim222[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim22,truevalssim2,response,explanatory,frac=frac2,maxits=maxits2,datinpool=rubindatasets)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim223[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim22,truevalssim2,response,explanatory,frac=frac2,maxits=maxits3,datinpool=rubindatasets)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim231[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim23,truevalssim2,response,explanatory,frac=frac3,maxits=maxits1,datinpool=rubindatasets)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim232[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim23,truevalssim2,response,explanatory,frac=frac3,maxits=maxits2,datinpool=rubindatasets)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim233[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim23,truevalssim2,response,explanatory,frac=frac3,maxits=maxits3,datinpool=rubindatasets)
      
      
      
      #rubin-validating simulations control
      
      rubinlistsimcont1[[length(rubinlistsimcont1)+1]]=rubinvalidator(imputedsimcont11[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim11,truevalssim1,response,explanatory,frac=frac1,maxits=1,datinpool=rubindatasets)
      
      rubinlistsimcont1[[length(rubinlistsimcont1)+1]]=rubinvalidator(imputedsimcont12[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim12,truevalssim1,response,explanatory,frac=frac2,maxits=1,datinpool=rubindatasets)
      
      rubinlistsimcont1[[length(rubinlistsimcont1)+1]]=rubinvalidator(imputedsimcont13[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim13,truevalssim1,response,explanatory,frac=frac3,maxits=1,datinpool=rubindatasets)
      
      rubinlistsimcont2[[length(rubinlistsimcont2)+1]]=rubinvalidator(imputedsimcont21[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim21,truevalssim2,response,explanatory,frac=frac1,maxits=1,datinpool=rubindatasets)
      
      rubinlistsimcont2[[length(rubinlistsimcont2)+1]]=rubinvalidator(imputedsimcont22[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim22,truevalssim2,response,explanatory,frac=frac2,maxits=1,datinpool=rubindatasets)
      
      rubinlistsimcont2[[length(rubinlistsimcont2)+1]]=rubinvalidator(imputedsimcont23[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafreesim23,truevalssim2,response,explanatory,frac=frac3,maxits=1,datinpool=rubindatasets)
      
      
      
      #rubin-validating real algorithm
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed11[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafree1,truevals,response,explanatory,frac=frac1,maxits=maxits1,datinpool=rubindatasets,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed12[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafree1,truevals,response,explanatory,frac=frac1,maxits=maxits2,datinpool=rubindatasets,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed13[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafree1,truevals,response,explanatory,frac=frac1,maxits=maxits3,datinpool=rubindatasets,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed21[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafree2,truevals,response,explanatory,frac=frac2,maxits=maxits1,datinpool=rubindatasets,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed22[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafree2,truevals,response,explanatory,frac=frac2,maxits=maxits2,datinpool=rubindatasets,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed23[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafree2,truevals,response,explanatory,frac=frac2,maxits=maxits3,datinpool=rubindatasets,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed31[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafree3,truevals,response,explanatory,frac=frac3,maxits=maxits1,datinpool=rubindatasets,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed32[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafree3,truevals,response,explanatory,frac=frac3,maxits=maxits2,datinpool=rubindatasets,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed33[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafree3,truevals,response,explanatory,frac=frac3,maxits=maxits3,datinpool=rubindatasets,logtrans=logtransform)
      
      
      
      #rubin-validating real control
      
      rubinlistcont[[length(rubinlistcont)+1]]=rubinvalidator(imputedcont1[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafree1,truevals,response,explanatory,frac=frac1,maxits=1,datinpool=rubindatasets,logtrans=logtransform)
      
      rubinlistcont[[length(rubinlistcont)+1]]=rubinvalidator(imputedcont2[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafree2,truevals,response,explanatory,frac=frac2,maxits=1,datinpool=rubindatasets,logtrans=logtransform)
      
      rubinlistcont[[length(rubinlistcont)+1]]=rubinvalidator(imputedcont3[(1+(xa-1)*rubindatasets):(xa*rubindatasets)],datafree3,truevals,response,explanatory,frac=frac3,maxits=1,datinpool=rubindatasets,logtrans=logtransform)
      
      
      
      #rubin-validating simulated mice
      
      micedats=parLapply(clu,1:rubindatasets,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim11,impmeth=impmeth)
      
      rubinlistmicesim1[[length(rubinlistmicesim1)+1]]=rubinvalidator(micedats,datasim11,truevalssim1,response,explanatory,frac=frac1,maxits=0,datinpool=rubindatasets)
      
      micedats=parLapply(clu,1:rubindatasets,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim12,impmeth=impmeth)
      
      rubinlistmicesim1[[length(rubinlistmicesim1)+1]]=rubinvalidator(micedats,datasim12,truevalssim1,response,explanatory,frac=frac2,maxits=0,datinpool=rubindatasets)
      
      micedats=parLapply(clu,1:rubindatasets,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim13,impmeth=impmeth)
      
      rubinlistmicesim1[[length(rubinlistmicesim1)+1]]=rubinvalidator(micedats,datasim13,truevalssim1,response,explanatory,frac=frac3,maxits=0,datinpool=rubindatasets)
      
      micedats=parLapply(clu,1:rubindatasets,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim21,impmeth=impmeth)
      
      rubinlistmicesim2[[length(rubinlistmicesim2)+1]]=rubinvalidator(micedats,datasim21,truevalssim2,response,explanatory,frac=frac1,maxits=0,datinpool=rubindatasets)
      
      micedats=parLapply(clu,1:rubindatasets,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim22,impmeth=impmeth)
      
      rubinlistmicesim2[[length(rubinlistmicesim2)+1]]=rubinvalidator(micedats,datasim21,truevalssim2,response,explanatory,frac=frac2,maxits=0,datinpool=rubindatasets)
      
      micedats=parLapply(clu,1:rubindatasets,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim23,impmeth=impmeth)
      
      rubinlistmicesim2[[length(rubinlistmicesim2)+1]]=rubinvalidator(micedats,datasim23,truevalssim2,response,explanatory,frac=frac3,maxits=0,datinpool=rubindatasets)
      
      
      
      #rubin-validating real mice
      
      micedats=parLapply(clu,1:rubindatasets,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=data1,impmeth=impmeth)
      
      rubinlistmice[[length(rubinlistmice)+1]]=rubinvalidator(micedats,data1,truevals,response,explanatory,frac=frac1,maxits=0,datinpool=rubindatasets,logtrans=logtransform)
      
      micedats=parLapply(clu,1:rubindatasets,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=data2,impmeth=impmeth)
      
      rubinlistmice[[length(rubinlistmice)+1]]=rubinvalidator(micedats,data2,truevals,response,explanatory,frac=frac2,maxits=0,datinpool=rubindatasets,logtrans=logtransform)
      
      micedats=parLapply(clu,1:rubindatasets,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=data3,impmeth=impmeth)
      
      rubinlistmice[[length(rubinlistmice)+1]]=rubinvalidator(micedats,data3,truevals,response,explanatory,frac=frac3,maxits=0,datinpool=rubindatasets,logtrans=logtransform)
      
    }
    
    
    
    #rubin-validating with rubindatasets2
    
    for (xb in 1:(length(imputedsim111)%/%rubindatasets2)) {
      
      
      
      #rubin-validating simulations algorithm
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim111[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim11,truevalssim1,response,explanatory,frac=frac1,maxits=maxits1,datinpool=rubindatasets2)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim112[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim11,truevalssim1,response,explanatory,frac=frac1,maxits=maxits2,datinpool=rubindatasets2)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim113[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim11,truevalssim1,response,explanatory,frac=frac1,maxits=maxits3,datinpool=rubindatasets2)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim121[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim12,truevalssim1,response,explanatory,frac=frac2,maxits=maxits1,datinpool=rubindatasets2)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim122[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim12,truevalssim1,response,explanatory,frac=frac2,maxits=maxits2,datinpool=rubindatasets2)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim123[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim12,truevalssim1,response,explanatory,frac=frac2,maxits=maxits3,datinpool=rubindatasets2)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim131[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim13,truevalssim1,response,explanatory,frac=frac3,maxits=maxits1,datinpool=rubindatasets2)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim132[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim13,truevalssim1,response,explanatory,frac=frac3,maxits=maxits2,datinpool=rubindatasets2)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim133[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim13,truevalssim1,response,explanatory,frac=frac3,maxits=maxits3,datinpool=rubindatasets2)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim211[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim21,truevalssim2,response,explanatory,frac=frac1,maxits=maxits1,datinpool=rubindatasets2)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim212[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim21,truevalssim2,response,explanatory,frac=frac1,maxits=maxits2,datinpool=rubindatasets2)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim213[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim21,truevalssim2,response,explanatory,frac=frac1,maxits=maxits3,datinpool=rubindatasets2)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim221[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim22,truevalssim2,response,explanatory,frac=frac2,maxits=maxits1,datinpool=rubindatasets2)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim222[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim22,truevalssim2,response,explanatory,frac=frac2,maxits=maxits2,datinpool=rubindatasets2)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim223[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim22,truevalssim2,response,explanatory,frac=frac2,maxits=maxits3,datinpool=rubindatasets2)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim231[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim23,truevalssim2,response,explanatory,frac=frac3,maxits=maxits1,datinpool=rubindatasets2)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim232[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim23,truevalssim2,response,explanatory,frac=frac3,maxits=maxits2,datinpool=rubindatasets2)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim233[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim23,truevalssim2,response,explanatory,frac=frac3,maxits=maxits3,datinpool=rubindatasets2)
      
      
      
      #rubin-validating simulations control
      
      rubinlistsimcont1[[length(rubinlistsimcont1)+1]]=rubinvalidator(imputedsimcont11[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim11,truevalssim1,response,explanatory,frac=frac1,maxits=1,datinpool=rubindatasets2)
      
      rubinlistsimcont1[[length(rubinlistsimcont1)+1]]=rubinvalidator(imputedsimcont12[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim12,truevalssim1,response,explanatory,frac=frac2,maxits=1,datinpool=rubindatasets2)
      
      rubinlistsimcont1[[length(rubinlistsimcont1)+1]]=rubinvalidator(imputedsimcont13[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim13,truevalssim1,response,explanatory,frac=frac3,maxits=1,datinpool=rubindatasets2)
      
      rubinlistsimcont2[[length(rubinlistsimcont2)+1]]=rubinvalidator(imputedsimcont21[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim21,truevalssim2,response,explanatory,frac=frac1,maxits=1,datinpool=rubindatasets2)
      
      rubinlistsimcont2[[length(rubinlistsimcont2)+1]]=rubinvalidator(imputedsimcont22[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim22,truevalssim2,response,explanatory,frac=frac2,maxits=1,datinpool=rubindatasets2)
      
      rubinlistsimcont2[[length(rubinlistsimcont2)+1]]=rubinvalidator(imputedsimcont23[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafreesim23,truevalssim2,response,explanatory,frac=frac3,maxits=1,datinpool=rubindatasets2)
      
      
      
      #rubin-validating real algorithm
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed11[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafree1,truevals,response,explanatory,frac=frac1,maxits=maxits1,datinpool=rubindatasets2,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed12[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafree1,truevals,response,explanatory,frac=frac1,maxits=maxits2,datinpool=rubindatasets2,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed13[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafree1,truevals,response,explanatory,frac=frac1,maxits=maxits3,datinpool=rubindatasets2,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed21[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafree2,truevals,response,explanatory,frac=frac2,maxits=maxits1,datinpool=rubindatasets2,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed22[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafree2,truevals,response,explanatory,frac=frac2,maxits=maxits2,datinpool=rubindatasets2,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed23[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafree2,truevals,response,explanatory,frac=frac2,maxits=maxits3,datinpool=rubindatasets2,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed31[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafree3,truevals,response,explanatory,frac=frac3,maxits=maxits1,datinpool=rubindatasets2,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed32[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafree3,truevals,response,explanatory,frac=frac3,maxits=maxits2,datinpool=rubindatasets2,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed33[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafree3,truevals,response,explanatory,frac=frac3,maxits=maxits3,datinpool=rubindatasets2,logtrans=logtransform)
      
      
      
      #rubin-validating real control
      
      rubinlistcont[[length(rubinlistcont)+1]]=rubinvalidator(imputedcont1[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafree1,truevals,response,explanatory,frac=frac1,maxits=1,datinpool=rubindatasets2,logtrans=logtransform)
      
      rubinlistcont[[length(rubinlistcont)+1]]=rubinvalidator(imputedcont2[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafree2,truevals,response,explanatory,frac=frac2,maxits=1,datinpool=rubindatasets2,logtrans=logtransform)
      
      rubinlistcont[[length(rubinlistcont)+1]]=rubinvalidator(imputedcont3[(1+(xb-1)*rubindatasets2):(xb*rubindatasets2)],datafree3,truevals,response,explanatory,frac=frac3,maxits=1,datinpool=rubindatasets2,logtrans=logtransform)
      
      
      
      #rubin-validating simulated mice
      
      micedats=parLapply(clu,1:rubindatasets2,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim11,impmeth=impmeth)
      
      rubinlistmicesim1[[length(rubinlistmicesim1)+1]]=rubinvalidator(micedats,datasim11,truevalssim1,response,explanatory,frac=frac1,maxits=0,datinpool=rubindatasets2)
      
      micedats=parLapply(clu,1:rubindatasets2,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim12,impmeth=impmeth)
      
      rubinlistmicesim1[[length(rubinlistmicesim1)+1]]=rubinvalidator(micedats,datasim12,truevalssim1,response,explanatory,frac=frac2,maxits=0,datinpool=rubindatasets2)
      
      micedats=parLapply(clu,1:rubindatasets2,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim13,impmeth=impmeth)
      
      rubinlistmicesim1[[length(rubinlistmicesim1)+1]]=rubinvalidator(micedats,datasim13,truevalssim1,response,explanatory,frac=frac3,maxits=0,datinpool=rubindatasets2)
      
      micedats=parLapply(clu,1:rubindatasets2,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim21,impmeth=impmeth)
      
      rubinlistmicesim2[[length(rubinlistmicesim2)+1]]=rubinvalidator(micedats,datasim21,truevalssim2,response,explanatory,frac=frac1,maxits=0,datinpool=rubindatasets2)
      
      micedats=parLapply(clu,1:rubindatasets2,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim22,impmeth=impmeth)
      
      rubinlistmicesim2[[length(rubinlistmicesim2)+1]]=rubinvalidator(micedats,datasim21,truevalssim2,response,explanatory,frac=frac2,maxits=0,datinpool=rubindatasets2)
      
      micedats=parLapply(clu,1:rubindatasets2,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim23,impmeth=impmeth)
      
      rubinlistmicesim2[[length(rubinlistmicesim2)+1]]=rubinvalidator(micedats,datasim23,truevalssim2,response,explanatory,frac=frac3,maxits=0,datinpool=rubindatasets2)
      
      
      
      #rubin-validating real mice
      
      micedats=parLapply(clu,1:rubindatasets2,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=data1,impmeth=impmeth)
      
      rubinlistmice[[length(rubinlistmice)+1]]=rubinvalidator(micedats,data1,truevals,response,explanatory,frac=frac1,maxits=0,datinpool=rubindatasets2,logtrans=logtransform)
      
      micedats=parLapply(clu,1:rubindatasets2,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=data2,impmeth=impmeth)
      
      rubinlistmice[[length(rubinlistmice)+1]]=rubinvalidator(micedats,data2,truevals,response,explanatory,frac=frac2,maxits=0,datinpool=rubindatasets2,logtrans=logtransform)
      
      micedats=parLapply(clu,1:rubindatasets2,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=data3,impmeth=impmeth)
      
      rubinlistmice[[length(rubinlistmice)+1]]=rubinvalidator(micedats,data3,truevals,response,explanatory,frac=frac3,maxits=0,datinpool=rubindatasets2,logtrans=logtransform)
      
    }
    
    
    
    #rubin-validating with rubindatasets3
    
    for (xc in 1:(length(imputedsim111)%/%rubindatasets3)) {
      
      
      
      #rubin-validating simulations algorithm
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim111[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim11,truevalssim1,response,explanatory,frac=frac1,maxits=maxits1,datinpool=rubindatasets3)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim112[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim11,truevalssim1,response,explanatory,frac=frac1,maxits=maxits2,datinpool=rubindatasets3)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim113[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim11,truevalssim1,response,explanatory,frac=frac1,maxits=maxits3,datinpool=rubindatasets3)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim121[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim12,truevalssim1,response,explanatory,frac=frac2,maxits=maxits1,datinpool=rubindatasets3)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim122[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim12,truevalssim1,response,explanatory,frac=frac2,maxits=maxits2,datinpool=rubindatasets3)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim123[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim12,truevalssim1,response,explanatory,frac=frac2,maxits=maxits3,datinpool=rubindatasets3)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim131[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim13,truevalssim1,response,explanatory,frac=frac3,maxits=maxits1,datinpool=rubindatasets3)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim132[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim13,truevalssim1,response,explanatory,frac=frac3,maxits=maxits2,datinpool=rubindatasets3)
      
      rubinlistsim1[[length(rubinlistsim1)+1]]=rubinvalidator(imputedsim133[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim13,truevalssim1,response,explanatory,frac=frac3,maxits=maxits3,datinpool=rubindatasets3)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim211[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim21,truevalssim2,response,explanatory,frac=frac1,maxits=maxits1,datinpool=rubindatasets3)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim212[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim21,truevalssim2,response,explanatory,frac=frac1,maxits=maxits2,datinpool=rubindatasets3)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim213[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim21,truevalssim2,response,explanatory,frac=frac1,maxits=maxits3,datinpool=rubindatasets3)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim221[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim22,truevalssim2,response,explanatory,frac=frac2,maxits=maxits1,datinpool=rubindatasets3)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim222[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim22,truevalssim2,response,explanatory,frac=frac2,maxits=maxits2,datinpool=rubindatasets3)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim223[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim22,truevalssim2,response,explanatory,frac=frac2,maxits=maxits3,datinpool=rubindatasets3)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim231[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim23,truevalssim2,response,explanatory,frac=frac3,maxits=maxits1,datinpool=rubindatasets3)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim232[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim23,truevalssim2,response,explanatory,frac=frac3,maxits=maxits2,datinpool=rubindatasets3)
      
      rubinlistsim2[[length(rubinlistsim2)+1]]=rubinvalidator(imputedsim233[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim23,truevalssim2,response,explanatory,frac=frac3,maxits=maxits3,datinpool=rubindatasets3)
      
      
      
      #rubin-validating simulations control
      
      rubinlistsimcont1[[length(rubinlistsimcont1)+1]]=rubinvalidator(imputedsimcont11[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim11,truevalssim1,response,explanatory,frac=frac1,maxits=1,datinpool=rubindatasets3)
      
      rubinlistsimcont1[[length(rubinlistsimcont1)+1]]=rubinvalidator(imputedsimcont12[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim12,truevalssim1,response,explanatory,frac=frac2,maxits=1,datinpool=rubindatasets3)
      
      rubinlistsimcont1[[length(rubinlistsimcont1)+1]]=rubinvalidator(imputedsimcont13[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim13,truevalssim1,response,explanatory,frac=frac3,maxits=1,datinpool=rubindatasets3)
      
      rubinlistsimcont2[[length(rubinlistsimcont2)+1]]=rubinvalidator(imputedsimcont21[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim21,truevalssim2,response,explanatory,frac=frac1,maxits=1,datinpool=rubindatasets3)
      
      rubinlistsimcont2[[length(rubinlistsimcont2)+1]]=rubinvalidator(imputedsimcont22[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim22,truevalssim2,response,explanatory,frac=frac2,maxits=1,datinpool=rubindatasets3)
      
      rubinlistsimcont2[[length(rubinlistsimcont2)+1]]=rubinvalidator(imputedsimcont23[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafreesim23,truevalssim2,response,explanatory,frac=frac3,maxits=1,datinpool=rubindatasets3)
      
      
      
      #rubin-validating real algorithm
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed11[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafree1,truevals,response,explanatory,frac=frac1,maxits=maxits1,datinpool=rubindatasets3,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed12[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafree1,truevals,response,explanatory,frac=frac1,maxits=maxits2,datinpool=rubindatasets3,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed13[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafree1,truevals,response,explanatory,frac=frac1,maxits=maxits3,datinpool=rubindatasets3,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed21[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafree2,truevals,response,explanatory,frac=frac2,maxits=maxits1,datinpool=rubindatasets3,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed22[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafree2,truevals,response,explanatory,frac=frac2,maxits=maxits2,datinpool=rubindatasets3,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed23[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafree2,truevals,response,explanatory,frac=frac2,maxits=maxits3,datinpool=rubindatasets3,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed31[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafree3,truevals,response,explanatory,frac=frac3,maxits=maxits1,datinpool=rubindatasets3,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed32[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafree3,truevals,response,explanatory,frac=frac3,maxits=maxits2,datinpool=rubindatasets3,logtrans=logtransform)
      
      rubinlist[[length(rubinlist)+1]]=rubinvalidator(imputed33[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafree3,truevals,response,explanatory,frac=frac3,maxits=maxits3,datinpool=rubindatasets3,logtrans=logtransform)
      
      
      
      #rubin-validating real control
      
      rubinlistcont[[length(rubinlistcont)+1]]=rubinvalidator(imputedcont1[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafree1,truevals,response,explanatory,frac=frac1,maxits=1,datinpool=rubindatasets3,logtrans=logtransform)
      
      rubinlistcont[[length(rubinlistcont)+1]]=rubinvalidator(imputedcont2[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafree2,truevals,response,explanatory,frac=frac2,maxits=1,datinpool=rubindatasets3,logtrans=logtransform)
      
      rubinlistcont[[length(rubinlistcont)+1]]=rubinvalidator(imputedcont3[(1+(xc-1)*rubindatasets3):(xc*rubindatasets3)],datafree3,truevals,response,explanatory,frac=frac3,maxits=1,datinpool=rubindatasets3,logtrans=logtransform)
      
      
      
      #rubin-validating simulated mice
      
      micedats=parLapply(clu,1:rubindatasets3,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim11,impmeth=impmeth)
      
      rubinlistmicesim1[[length(rubinlistmicesim1)+1]]=rubinvalidator(micedats,datasim11,truevalssim1,response,explanatory,frac=frac1,maxits=0,datinpool=rubindatasets3)
      
      micedats=parLapply(clu,1:rubindatasets3,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim12,impmeth=impmeth)
      
      rubinlistmicesim1[[length(rubinlistmicesim1)+1]]=rubinvalidator(micedats,datasim12,truevalssim1,response,explanatory,frac=frac2,maxits=0,datinpool=rubindatasets3)
      
      micedats=parLapply(clu,1:rubindatasets3,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim13,impmeth=impmeth)
      
      rubinlistmicesim1[[length(rubinlistmicesim1)+1]]=rubinvalidator(micedats,datasim13,truevalssim1,response,explanatory,frac=frac3,maxits=0,datinpool=rubindatasets3)
      
      micedats=parLapply(clu,1:rubindatasets3,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim21,impmeth=impmeth)
      
      rubinlistmicesim2[[length(rubinlistmicesim2)+1]]=rubinvalidator(micedats,datasim21,truevalssim2,response,explanatory,frac=frac1,maxits=0,datinpool=rubindatasets3)
      
      micedats=parLapply(clu,1:rubindatasets3,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim22,impmeth=impmeth)
      
      rubinlistmicesim2[[length(rubinlistmicesim2)+1]]=rubinvalidator(micedats,datasim21,truevalssim2,response,explanatory,frac=frac2,maxits=0,datinpool=rubindatasets3)
      
      micedats=parLapply(clu,1:rubindatasets3,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=datasim23,impmeth=impmeth)
      
      rubinlistmicesim2[[length(rubinlistmicesim2)+1]]=rubinvalidator(micedats,datasim23,truevalssim2,response,explanatory,frac=frac3,maxits=0,datinpool=rubindatasets3)
      
      
      
      #rubin-validating real mice
      
      micedats=parLapply(clu,1:rubindatasets3,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=data1,impmeth=impmeth)
      
      rubinlistmice[[length(rubinlistmice)+1]]=rubinvalidator(micedats,data1,truevals,response,explanatory,frac=frac1,maxits=0,datinpool=rubindatasets3,logtrans=logtransform)
      
      micedats=parLapply(clu,1:rubindatasets3,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=data2,impmeth=impmeth)
      
      rubinlistmice[[length(rubinlistmice)+1]]=rubinvalidator(micedats,data2,truevals,response,explanatory,frac=frac2,maxits=0,datinpool=rubindatasets3,logtrans=logtransform)
      
      micedats=parLapply(clu,1:rubindatasets3,function(x,data,impmeth){library(mice);complete(suppressWarnings(mice(data,1,impmeth,threshold=0.999999999999999,print=FALSE)))},data=data3,impmeth=impmeth)
      
      rubinlistmice[[length(rubinlistmice)+1]]=rubinvalidator(micedats,data3,truevals,response,explanatory,frac=frac3,maxits=0,datinpool=rubindatasets3,logtrans=logtransform)
      
    }
    
    
    
    rubinlistsim1=rubincleaner(rubinlistsim1)
    
    rubinlistsim2=rubincleaner(rubinlistsim2)
    
    rubinlistsimcont1=rubincleaner(rubinlistsimcont1)
    
    rubinlistsimcont2=rubincleaner(rubinlistsimcont2)
    
    rubinlistmicesim1=rubincleaner(rubinlistmicesim1)
    
    rubinlistmicesim2=rubincleaner(rubinlistmicesim2)
    
    rubinlist=rubincleaner(rubinlist)
    
    rubinlistcont=rubincleaner(rubinlistcont)
    
    rubinlistmice=rubincleaner(rubinlistmice)
    
    
    
    imputedsim111=imputedsim111[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim112=imputedsim112[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim113=imputedsim113[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim121=imputedsim121[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim122=imputedsim122[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim123=imputedsim123[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim131=imputedsim131[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim132=imputedsim132[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim133=imputedsim133[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim211=imputedsim211[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim212=imputedsim212[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim213=imputedsim213[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim221=imputedsim221[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim222=imputedsim222[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim223=imputedsim223[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim231=imputedsim231[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim232=imputedsim232[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsim233=imputedsim233[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    
    
    imputedsimcont11=imputedsimcont11[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsimcont12=imputedsimcont12[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsimcont13=imputedsimcont13[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsimcont21=imputedsimcont21[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsimcont22=imputedsimcont22[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedsimcont23=imputedsimcont23[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    
    
    imputed11=imputed11[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputed12=imputed12[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputed13=imputed13[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputed21=imputed21[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputed22=imputed22[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputed23=imputed23[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputed31=imputed31[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputed32=imputed32[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputed33=imputed33[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    
    
    imputedcont1=imputedcont1[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedcont2=imputedcont2[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
    imputedcont3=imputedcont3[-c(1:max(xa*rubindatasets,xb*rubindatasets2,xc*rubindatasets3))]
    
  }
  
  
  
  if ((it/parallelimp)%%(holecycle*datbatches)==0|difftime(Sys.time(),starttime,units = "hour")>=maxhours) {
    
    datasetnr=ceiling((it/parallelimp)/(holecycle*datbatches))
    
    dirnamemain=paste0("./",round(as.numeric(starttime)),"it_",datasetnr)
    
    if ((it/parallelimp)/(holecycle*datbatches)!=ceiling((it/parallelimp)/(holecycle*datbatches))) {
      
      dirnamemain=paste0(dirnamemain,"#")
      
    }
    
    dirnamesim1=paste0(dirnamemain,"/algorithmsim1")
    
    dirnamesim2=paste0(dirnamemain,"/algorithmsim2")
    
    dirnamesimcont1=paste0(dirnamemain,"/algorithmsimcont1")
    
    dirnamesimcont2=paste0(dirnamemain,"/algorithmsimcont2")
    
    dirname=paste0(dirnamemain,"/algorithm")
    
    dirnamecont=paste0(dirnamemain,"/algorithmcont")
    
    dir.create(dirnamemain)
    
    dir.create(dirnamesim1)
    
    dir.create(dirnamesim2)
    
    dir.create(dirnamesimcont1)
    
    dir.create(dirnamesimcont2)
    
    dir.create(dirname)
    
    dir.create(dirnamecont)
    
    
    
    NAset=(rep(1:datbatches,each=holecycle*parallelimp*9)+(datasetnr-1)*datbatches)[1:length(outlistsim1$covdifmed)]
    
    NAsetcont=(rep(1:datbatches,each=holecycle*parallelimp*3)+(datasetnr-1)*datbatches)[1:length(outlistsimcont1$covdifmed)]
    
    
    
    valid=cbind(dataset=datasetnr*2-1,NAset=NAset,frac=outlistsim1$frac,maxits=outlistsim1$maxits,reldif=outlistsim1$reldif,covdifmed=abs(outlistsim1$covdifmed),editerr=outlistsim1$editerr)
    
    write.csv(valid,paste0(dirnamesim1,"/validation.csv"))
    
    write.csv(outlistsim1$covdif,paste0(dirnamesim1,"/covdif.csv"))
    
    write.csv(outlistsim1$vardif,paste0(dirnamesim1,"/vardif.csv"))
    
    write.csv(outlistsim1$kstests,paste0(dirnamesim1,"/kstests.csv"))
    
    write.csv(outlistsim1$meansconf,paste0(dirnamesim1,"/meansconf.csv"))
    
    write.csv(outlistsim1$totmiss,paste0(dirnamesim1,"/totmiss.csv"))
    
    rubinlistsim1=lapply(rubinlistsim1,function(x){rbind(x,-1)})
    
    rubinlistsim1=try(do.call(rbind,rubinlistsim1))
    
    if (!class(rubinlistsim1)%in%c("NULL","try-error")) {write.csv(rubinlistsim1,paste0(dirnamesim1,"/rubin_parameters.csv"))}
    
    rubinlistmicesim1=lapply(rubinlistmicesim1,function(x){rbind(x,-1)})
    
    if (!class(rubinlistsim1)%in%c("NULL","try-error")) {rubinlistmicesim1=do.call(rbind,rubinlistmicesim1);write.csv(rubinlistmicesim1,paste0(dirnamesim1,"/rubin_parameters_mice.csv"))}
    
    
    
    valid=cbind(dataset=datasetnr*2,NAset=NAset,frac=outlistsim2$frac,maxits=outlistsim2$maxits,reldif=outlistsim2$reldif,covdifmed=abs(outlistsim2$covdifmed),editerr=outlistsim2$editerr)
    
    write.csv(valid,paste0(dirnamesim2,"/validation.csv"))
    
    write.csv(outlistsim2$covdif,paste0(dirnamesim2,"/covdif.csv"))
    
    write.csv(outlistsim2$vardif,paste0(dirnamesim2,"/vardif.csv"))
    
    write.csv(outlistsim2$kstests,paste0(dirnamesim2,"/kstests.csv"))
    
    write.csv(outlistsim2$meansconf,paste0(dirnamesim2,"/meansconf.csv"))
    
    write.csv(outlistsim2$totmiss,paste0(dirnamesim2,"/totmiss.csv"))
    
    rubinlistsim2=lapply(rubinlistsim2,function(x){rbind(x,-1)})
    
    if (!class(rubinlistsim1)%in%c("NULL","try-error")) {rubinlistsim2=do.call(rbind,rubinlistsim2);write.csv(rubinlistsim2,paste0(dirnamesim2,"/rubin_parameters.csv"))}
    
    rubinlistmicesim2=lapply(rubinlistmicesim2,function(x){rbind(x,-1)})
    
    if (!class(rubinlistsim1)%in%c("NULL","try-error")) {rubinlistmicesim2=do.call(rbind,rubinlistmicesim2);write.csv(rubinlistmicesim2,paste0(dirnamesim2,"/rubin_parameters_mice.csv"))}
    
    
    
    valid=cbind(dataset=datasetnr*2-1,NAset=NAsetcont,frac=outlistsimcont1$frac,maxits=outlistsimcont1$maxits,reldif=outlistsimcont1$reldif,covdifmed=abs(outlistsimcont1$covdifmed),editerr=outlistsimcont1$editerr)
    
    write.csv(valid,paste0(dirnamesimcont1,"/validation.csv"))
    
    write.csv(outlistsimcont1$covdif,paste0(dirnamesimcont1,"/covdif.csv"))
    
    write.csv(outlistsimcont1$vardif,paste0(dirnamesimcont1,"/vardif.csv"))
    
    write.csv(outlistsimcont1$kstests,paste0(dirnamesimcont1,"/kstests.csv"))
    
    write.csv(outlistsimcont1$meansconf,paste0(dirnamesimcont1,"/meansconf.csv"))
    
    write.csv(outlistsimcont1$totmiss,paste0(dirnamesimcont1,"/totmiss.csv"))
    
    rubinlistsimcont1=lapply(rubinlistsimcont1,function(x){rbind(x,-1)})
    
    if (!class(rubinlistsim1)%in%c("NULL","try-error")) {rubinlistsimcont1=do.call(rbind,rubinlistsimcont1);write.csv(rubinlistsimcont1,paste0(dirnamesimcont1,"/rubin_parameters.csv"))}
    
    
    
    valid=cbind(dataset=datasetnr*2,NAset=NAsetcont,frac=outlistsimcont2$frac,maxits=outlistsimcont2$maxits,reldif=outlistsimcont2$reldif,covdifmed=abs(outlistsimcont2$covdifmed),editerr=outlistsimcont2$editerr)
    
    write.csv(valid,paste0(dirnamesimcont2,"/validation.csv"))
    
    write.csv(outlistsimcont2$covdif,paste0(dirnamesimcont2,"/covdif.csv"))
    
    write.csv(outlistsimcont2$vardif,paste0(dirnamesimcont2,"/vardif.csv"))
    
    write.csv(outlistsimcont2$kstests,paste0(dirnamesimcont2,"/kstests.csv"))
    
    write.csv(outlistsimcont2$meansconf,paste0(dirnamesimcont2,"/meansconf.csv"))
    
    write.csv(outlistsimcont2$totmiss,paste0(dirnamesimcont2,"/totmiss.csv"))
    
    rubinlistsimcont2=lapply(rubinlistsimcont2,function(x){rbind(x,-1)})
    
    if (!class(rubinlistsim1)%in%c("NULL","try-error")) {rubinlistsimcont2=do.call(rbind,rubinlistsimcont2);write.csv(rubinlistsimcont2,paste0(dirnamesimcont2,"/rubin_parameters.csv"))}
    
    
    
    valid=cbind(dataset=0,NAset=NAset,frac=outlist$frac,maxits=outlist$maxits,reldif=outlist$reldif,covdifmed=abs(outlist$covdifmed),editerr=outlist$editerr)
    
    write.csv(valid,paste0(dirname,"/validation.csv"))
    
    write.csv(outlist$covdif,paste0(dirname,"/covdif.csv"))
    
    write.csv(outlist$vardif,paste0(dirname,"/vardif.csv"))
    
    write.csv(outlist$kstests,paste0(dirname,"/kstests.csv"))
    
    write.csv(outlist$meansconf,paste0(dirname,"/meansconf.csv"))
    
    write.csv(outlist$totmiss,paste0(dirname,"/totmiss.csv"))
    
    rubinlist=lapply(rubinlist,function(x){rbind(x,-1)})
    
    if (!class(rubinlistsim1)%in%c("NULL","try-error")) {rubinlist=do.call(rbind,rubinlist);write.csv(rubinlist,paste0(dirname,"/rubin_parameters.csv"))}
    
    rubinlistmice=lapply(rubinlistmice,function(x){rbind(x,-1)})
    
    if (!class(rubinlistsim1)%in%c("NULL","try-error")) {rubinlistmice=do.call(rbind,rubinlistmice);write.csv(rubinlistmice,paste0(dirname,"/rubin_parameters_mice.csv"))}
    
    
    
    valid=cbind(dataset=0,NAset=NAsetcont,frac=outlistcont$frac,maxits=outlistcont$maxits,reldif=outlistcont$reldif,covdifmed=abs(outlistcont$covdifmed),editerr=outlistcont$editerr)
    
    write.csv(valid,paste0(dirnamecont,"/validation.csv"))
    
    write.csv(outlistcont$covdif,paste0(dirnamecont,"/covdif.csv"))
    
    write.csv(outlistcont$vardif,paste0(dirnamecont,"/vardif.csv"))
    
    write.csv(outlistcont$kstests,paste0(dirnamecont,"/kstests.csv"))
    
    write.csv(outlistcont$meansconf,paste0(dirnamecont,"/meansconf.csv"))
    
    write.csv(outlistcont$totmiss,paste0(dirnamecont,"/totmiss.csv"))
    
    rubinlistcont=lapply(rubinlistcont,function(x){rbind(x,-1)})
    
    if (!class(rubinlistsim1)%in%c("NULL","try-error")) {rubinlistcont=do.call(rbind,rubinlistcont);write.csv(rubinlistcont,paste0(dirnamecont,"/rubin_parameters.csv"))}
    
    
    
    outlist=outlistcont=outlistsim1=outlistsim2=outlistsimcont1=outlistsimcont2=NULL
    
    imputedsim111=imputedsim112=imputedsim113=imputedsim121=imputedsim122=imputedsim123=imputedsim131=imputedsim132=imputedsim133=imputedsim211=imputedsim212=imputedsim213=imputedsim221=imputedsim222=imputedsim223=imputedsim231=imputedsim232=imputedsim233=list()
    
    imputedsimcont11=imputedsimcont12=imputedsimcont13=imputedsimcont21=imputedsimcont22=imputedsimcont23=list()
    
    imputed11=imputed12=imputed13=imputed21=imputed22=imputed23=imputed31=imputed32=imputed33=list()
    
    imputedcont1=imputedcont2=imputedcont3=list()
    
    rubinlist=rubinlistcont=rubinlistmice=rubinlistsim1=rubinlistsim2=rubinlistsimcont1=rubinlistsimcont2=rubinlistmicesim1=rubinlistmicesim2=list()
    
  }
  
  
  
  if (difftime(Sys.time(),starttime,units = "hour")>=maxhours) {
    
    break
    
  }
}

stopCluster(clu)
gc()
