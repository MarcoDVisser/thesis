 setwd("/home/marcov/Dropbox/Backup/Work/Thesis/")

## SETUP

## USE C++?
CPP <- TRUE

## USE PARALELL?
paral <- FALSE

################################################################################
## Pacala's 1997 models on dispersal and Janzen-Connell.
## adapted with NDD dispersal
################################################################################
## D community death rate
## S local conspecific death rate
## E fraction that disperses outward
## N number of iterations
## R initial species richness
## P number of patches 

NDDe <- function(x,E,a=0.15) SSlogis(x,E,a,-E/10)
###############################################################################
## R versions
################################################################################

runModel <- function(S,E,D=.05,N=1e4,R=50,P=1000,vizualize=FALSE){

    ## initialize model
    adults <- sample(as.character(1:R),P,replace=TRUE)
    ## even
   # adults <- rep(1:R,each=P/R)
    adults <- factor(adults,levels=1:R)
    

    richness <- vector("numeric",N)

    for(j in 1:N){
        deaths <- sample(1:P,round(P*D))

        oldadults <- adults
        
    for(i in deaths){
        ## for each patch with death calculate dispersal
         ## from natal and outside
        spN <- as.numeric(as.character(oldadults[i]))
        seeds <- immigrants <- E*(table(oldadults[-i])/(P-1))
        locals <- (1-E)*S  #* seq(10000,100,length.out=R)[spN]
        seeds[spN] <- (immigrants[spN]*S) + locals
        pseeds <- seeds/sum(seeds)
        patchwinner <- sample(1:R,1,prob=pseeds)
        adults[i] <- factor(patchwinner,levels=1:R)
    }
    SR <- table(adults)
    Div <- (SR/sum(SR))[SR!=0]
    
    richness[j] <- sum(table(adults)>0)
    cat("\t \r",signif(j/N,4)*100,"%", "R =",richness[j],
        " D=", -sum(Div*log(Div)), "R//t =", (richness[j]-R)/j)

        if(vizualize){
            plot(table(adults),xaxt="n",yaxt="n")
        }
}
    cat("\n")

    return(list(richness,data.frame(R=richness[N],S,E)))
    
}

################################################################################
## now the NDDe version
runModelNDDe <- function(S,E,D=.05,N=1e4,R=50,P=1000,vizualize=FALSE){

    ## initialize model
    adults <- sample(as.character(1:R),P,replace=TRUE)
    ## even
#    adults <- rep(1:R,each=P/R)

    adults <- factor(adults,levels=1:R)
    

    richness <- vector("numeric",N)

    for(j in 1:N){
        deaths <- sample(1:P,round(P*D))
        Ri <- table(adults)
        oldadults <- adults
    for(i in deaths){
        ## for each patch with death calculate dispersal
        ## from natal and outside
        spN <- as.numeric(as.character(oldadults[i]))
        Esp <- NDDe(Ri/sum(Ri),E=E)
        seeds <- immigrants <- Esp*(table(oldadults[-i])/(P-1)) #*seq(10000,100,length.out=R)
        locals <- (1-Esp[spN])*(S)  #* seq(10000,100,length.out=R)[spN]
        seeds[spN] <- (immigrants[spN]*(S)) + locals
        pseeds <- seeds/sum(seeds)
        patchwinner <- sample(1:R,1,prob=pseeds)
        adults[i] <- factor(patchwinner,levels=1:R)
    }
    SR <- table(adults)
    Div <- (SR/sum(SR))[SR!=0]
    
    richness[j] <- sum(table(adults)>0)
    cat("\t \r",signif(j/N,4)*100,"%", "R =",richness[j],
        " D=", -sum(Div*log(Div)), "R//t =", (richness[j]-R)/j)

        if(vizualize){
            plot(table(adults),xaxt="",yaxt="")
        }
}
    cat("\n")

    return(list(richness,data.frame(R=richness[N],S,E)))
    
}


################################################################################
## Cpp version
################################################################################

if(CPP){
## load packages 
require(Rcpp)
require(RcppArmadillo)
require(inline)
sourceCpp("./R/PacalaModel.cpp")
}


################################################################################

dims <- c(20,20)

evec <- seq(0.01,.9,length.out=dims[1])
svec <- seq(0,1,length.out=dims[2])

AllScens <- expand.grid(evec,svec)
colnames(AllScens) <- c("e","s")

Allres <- vector("list",nrow(AllScens))
Allresndd <- vector("list",nrow(AllScens))
Nini <- 2.5e5
Rini <- 50
Dini <- 0.3
Pini <- 1000

if(paral){
JOBS <- list(seq(1,400,40),seq(40,400,40))

children1 <- children2 <- list()

require(parallel)
ncores <- detectCores()/2


    
    for(i in 1:ncores){
        X <- JOBS[[1]][i]:JOBS[[2]][i]
    children1[[i]] <- mcparallel(runModel(S=AllScens$s[X],E=AllScens$e[X],N=Nini,
                                     D=Dini,R=Rini,P=Pini))
    }
    
Allres <- mccollect(children1)

    for(i in 1:ncores){
       X <- JOBS[[1]][i]:JOBS[[2]][i]
    children2[[i]] <- mcparallel(runModelNDDe(S=AllScens$s[X],E=AllScens$e[X],N=Nini,
                                     D=Dini,R=Rini,P=Pini))
}

Allresndd <- mccollect(children2)

    saveRDS(list(Allres,Allresndd),file="~/NDDres.rds")
    
}

if(!paral&!CPP){
for(scen in 1:nrow(AllScens)){

    cat("S=",AllScens$s[scen],"E=",AllScens$e[scen],"\n")
    set.seed(1)
    Allres[[scen]] <- runModel(S=AllScens$s[scen],E=AllScens$e[scen],N=Nini,
                               D=Dini,R=Rini,P=Pini)
    set.seed(1)
    Allresndd[[scen]] <- runModelNDDe(S=AllScens$s[scen],E=AllScens$e[scen],N=Nini,
                                  D=Dini,R=Rini,Pini)
}
}

if(CPP){
for(scen in 1:nrow(AllScens)){
    cat("S=",AllScens$s[scen],"E=",AllScens$e[scen],"\n")
    set.seed(1)
    Allres[[scen]] <- runModelcpp(S=AllScens$s[scen],E=AllScens$e[scen],N=Nini,
                                 D=round(Dini*Pini),R=0:(Rini-1),P=Pini)

    Allres[[scen]]$R <- sum(Allres[[scen]][[3]]>1)

    Allresndd[[scen]] <- runModelNDDecpp(S=AllScens$s[scen],E=AllScens$e[scen],
                                         N=Nini,
                                         D=round(Dini*Pini),R=0:(Rini-1),P=Pini)
    Allresndd[[scen]]$R <- sum(Allresndd[[scen]][[3]]>1)
    cat("Normal =", Allres[[scen]]$R,"NDDe =", Allresndd[[scen]]$R,"\t")
    cat(timestamp(),"\n")
}
}
 
AllScens$R <- sapply(Allres, function(X) X$R)
AllScens$Rndd <-  sapply(Allresndd, function(X) X$R)

saveRDS(list(Allres,Allresndd),file="./objects/PacalaModelRun.rds")
saveRDS(AllScens,file="./objects/PacalaModelRunSummary.rds")
saveRDS(dims,file="./objects/finaldimsPacalaModelRun.rds")
saveRDS(list(svec,evec),file="./objects/finalESModel.rds")

normat <- matrix(AllScens$R,ncol=dims[1],nrow=dims[2],byrow=TRUE)
rownames(normat) <- svec
colnames(normat) <- evec

nddmat <- matrix(AllScens$Rndd,ncol=dims[1],nrow=dims[2],byrow=TRUE)
rownames(nddmat) <- svec
colnames(nddmat) <- evec

plot(svec,normat[,1],type="l",lwd=2)
lines(svec,normat[,2],lwd=2,col="red")
lines(svec,normat[,3],lwd=2,col="blue")
lines(svec,normat[,4],lwd=2,col="green")



