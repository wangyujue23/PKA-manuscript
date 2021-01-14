EMU_Model <- list()
require(DEoptim)
require(nloptr)
require(MASS)
require(pracma)
require(numDeriv)


Test.Labeling <- read.csv(file="PKA data.csv")

######################################
group<-2
glcr.fcirc <-Test.Labeling[1,1+group]
lac.fcirc <-Test.Labeling[2,1+group]
glc.fcirc <-Test.Labeling[3,1+group]
E.lac.glcr <-Test.Labeling[4,1+group]
E.glc.glcr <-Test.Labeling[5,1+group]
E.glcr.lac <-Test.Labeling[6,1+group]
E.glc.lac <-Test.Labeling[7,1+group]
E.glcr.glc <-Test.Labeling[8,1+group]
E.lac.glc <-Test.Labeling[9,1+group]
E.total <- Test.Labeling[4:9,1+group]


Fluxify <- function(f) {
  
  v <- rep(0,12)
  v[2] <- f[1]
  v[4] <- f[3]
  v[7] <- f[4]
  v[11] <- f[5]
  v[12] <- f[6]
  v[1] <- glcr.fcirc - v[4] - v[7]
  v[9] <- glc.fcirc - v[7] - v[12]
  
  if(group%%2==1){
    v[3] <- f[2]
    v[6] <- lac.fcirc - v[3] - v[9]
  }
  else{
    v[6] <- f[2]
    v[3] <- lac.fcirc - v[6] - v[9]
  }
  
  v[8] <- glcr.fcirc - v[2] - v[3]
  v[10] <- glc.fcirc - v[8] - v[11]
  v[5] <- lac.fcirc - v[4] - v[10]

  return(v)
}


#########

MID <- function(f) {#Basic function
  v <- Fluxify(f)
  result<-rep(0,6)
  result[1]<- (v[3]+v[8]*v[9]/glc.fcirc)/lac.fcirc
  result[2]<- (v[8]+v[3]*v[10]/lac.fcirc)/glc.fcirc
  result[3]<- (v[4]+v[10]*v[7]/glc.fcirc)/glcr.fcirc
  result[4]<- (v[10]+v[4]*v[8]/glcr.fcirc)/glc.fcirc
  result[5]<- (v[7]+v[9]*v[4]/lac.fcirc)/glcr.fcirc
  result[6]<- (v[9]+v[7]*v[3]/glcr.fcirc)/lac.fcirc
  return(result)
}

###################
CI<-function(fluxset,threshold,fluxnumber){
  #upper boundary
  fluxlength<-length(fluxset)
  FluxPath <- matrix(0,nrow=1000,ncol=fluxlength+2)
  # colnames(FluxPath) <- c("1-Glycerol","2-Lactate","3-FFA","4-Protein1","5-Protein2","6-Protein3","7-Glycolysis","8-PDH","9-Suc_DH","10-Glc_tracer","11-Glcr_lac","12-Glcr_DHAP","Res","threshold")
  FluxPath[1,1:fluxlength] <- fluxset
  FluxPath[1,fluxlength+1] <- Res(fluxset)
  FluxPath[1,fluxlength+2] <- threshold
  CurrentFlux <- fluxset
  goodcount<-0
  highlimit<-0
  lowlimit<-0
  changingset<-CurrentFlux>0
  changingset[fluxnumber]<-FALSE
  
  for (i in 1:500) {
    J <- grad(Res,CurrentFlux)
    H <- hessian(Res,CurrentFlux)
    if(CurrentFlux[fluxnumber]<0.01){
      h <- 2*CurrentFlux[fluxnumber]
    } else {
      h <- 0.02
    }
    b <- -J[changingset]-H[changingset,fluxnumber]*h
    A <- H[changingset,changingset]
    u <- ginv(A,tol=1e-5) %*% b
    update <- rep(0,fluxlength)
    update[fluxnumber]<-h
    update[changingset]<-u[1:sum(changingset)]
    CurrentFlux <- CurrentFlux + update
    while(!prod(CurrentFlux>0)){
      changingset<-changingset*CurrentFlux>0
      CurrentFlux <- CurrentFlux - update
      if (sum(changingset)==0){
        break
        } else{
      b <- -J[changingset]-H[changingset,fluxnumber]*h
      A <- H[changingset,changingset]
      u <- ginv(A,tol=1e-5) %*% b
      update <- rep(0,fluxlength)
      update[fluxnumber]<-h
      update[changingset]<-u[1:sum(changingset)]
      CurrentFlux <- CurrentFlux + update
      }
    }
    FluxPath[1+i,1:fluxlength] <- CurrentFlux
    FluxPath[1+i,fluxlength+1] <- Res(CurrentFlux)
    if(Res(CurrentFlux)>=threshold){
      goodcount <- goodcount + 1
    }
    if(Res(CurrentFlux)<threshold){
      goodcount <- 0
      highlimit <- CurrentFlux[fluxnumber]
    }
    flush.console()
    if(i%%1==0) print(paste(FluxPath[1+i,]))
    if(sum(changingset)==0) break
    if(goodcount == 5) break
  }
  highPath<-FluxPath[1:(i+2),]
  
  #lower boundary
  FluxPath <- matrix(0,nrow=400,ncol=fluxlength+2)
  # colnames(FluxPath) <- c("1-Glycerol","2-Lactate","3-FFA","4-Protein1","5-Protein2","6-Protein3","7-Glycolysis","8-PDH","9-Suc_DH","10-Glc_tracer","11-Glcr_lac","12-Glcr_DHAP","Res","threshold")
  FluxPath[1,1:fluxlength] <- fluxset
  FluxPath[1,fluxlength+1] <- Res(fluxset)
  FluxPath[1,fluxlength+2] <- threshold
  CurrentFlux <- fluxset
  goodcount<-0
  changingset<-CurrentFlux>0
  changingset[fluxnumber]<-FALSE
  
  for (i in 1:500) {
    J <- grad(Res,CurrentFlux)
    H <- hessian(Res,CurrentFlux)
    if(CurrentFlux[fluxnumber]<0.02){
      break
    }
    else{
      h <- -0.02
    }
    b <- -J[changingset]-H[changingset,fluxnumber]*h
    A <- H[changingset,changingset]
    u <- ginv(A,tol=1e-5) %*% b
    update <- rep(0,fluxlength)
    update[fluxnumber]<-h
    update[changingset]<-u[1:sum(changingset)]
    CurrentFlux <- CurrentFlux + update
    while(!prod(CurrentFlux>0)){
      changingset<-changingset*CurrentFlux>0
      CurrentFlux <- CurrentFlux - update
      if (sum(changingset)==0){
        break
      } else{
        b <- -J[changingset]-H[changingset,fluxnumber]*h
        A <- H[changingset,changingset]
        u <- ginv(A,tol=1e-5) %*% b
        update <- rep(0,fluxlength)
        update[fluxnumber]<-h
        update[changingset]<-u[1:sum(changingset)]
        CurrentFlux <- CurrentFlux + update
      }
    }
    FluxPath[1+i,1:fluxlength] <- CurrentFlux
    FluxPath[1+i,fluxlength+1] <- Res(CurrentFlux)
    if(Res(CurrentFlux)>=threshold){
      goodcount <- goodcount + 1
    }
    if(Res(CurrentFlux)<threshold){
      goodcount <- 0
      lowlimit <- CurrentFlux[fluxnumber]
    }
    flush.console()
    if(i%%1==0) print(paste(FluxPath[1+i,]))
    if(sum(changingset)==0) break
    if(goodcount == 5) break
  }
  lowPath<-FluxPath[1:(i+1),]
  resultPath <- rbind(highPath,lowPath,c(highlimit-fluxset[fluxnumber],lowlimit-fluxset[fluxnumber],rep(0,fluxlength)))
  return(resultPath)
}

####################
Res <- function(f) {
  Label <- MID(f)
  Residual <- sum((E.total - Label)^2)
  return(Residual)
}

########
uplimit <- c(glcr.fcirc,
             lac.fcirc,
             min(glcr.fcirc,lac.fcirc),
             min(glc.fcirc,glcr.fcirc),
             glc.fcirc,
             glc.fcirc)

fluxsDE <- DEoptim(fn = Res,lower=rep(0.0001,6),upper=uplimit,control=list(itermax=2000))
write.csv(c(fluxsDE$optim$bestval,fluxsDE$optim$bestmem),file=paste(group,"_flux.csv",sep=""),row.names=TRUE)
bestflux<-fluxsDE$optim$bestmem
Fluxify((bestflux))


#6.24.2020 Fcirc from c13 glucose
#                       glcr_out    v3/v6       lac->glcr   glc->glcr   glc         glc_out   
#0.088239 bestmemit:    0.000100    0.962816    0.850332    0.000100    0.000100    0.000100     GFP               
#0.034556 bestmemit:    0.131477    0.000100    1.373058    0.000100    0.000100    0.229854     HF
#0.080254 bestmemit:    0.000100    0.000100    0.485897    0.000100    0.000100    0.440069     PKA
#0.037724 bestmemit:    0.000100    0.000100    1.326680    0.000100    0.000100    6.704199     PKAHF


threshold<-Res(bestflux)+0.0001*qchisq(c(0.8,0.9,0.95),1)[3]#Chi square threshold

for(fluxnumber in c(1:6)){
  write.csv(CI(bestflux,threshold,fluxnumber),file=paste("group",group,"_v",fluxnumber,".csv",sep=""))
}


flx <- read.csv(file="1_flux.csv")
flux1<-flx[2:7,2]
flx <- read.csv(file="2_flux.csv")
flux2<-flx[2:7,2]
flx <- read.csv(file="3_flux.csv")
flux3<-flx[2:7,2]
flx <- read.csv(file="4_flux.csv")
flux4<-flx[2:7,2]

fluxlist<- list(flux1,flux2,flux3,flux4)
######
fluxnames <-c(1:12,"Res")
Sum <- matrix(0,nrow=4,ncol=13,byrow=TRUE)
colnames(Sum) <- fluxnames
rownames(Sum) <- c("GFP","HF","PKA","PKAHF")

for (group in 1:4){
  glcr.fcirc <-Test.Labeling[1,1+group]
  lac.fcirc <-Test.Labeling[2,1+group]
  glc.fcirc <-Test.Labeling[3,1+group]
  E.total <- Test.Labeling[4:9,1+group]
  
  Sum[group,] <- c(Fluxify(fluxlist[[group]]),Res(fluxlist[[group]]))
}


write.csv(Sum,file=paste("Result",".csv",sep=""),row.names=TRUE)

###########################


