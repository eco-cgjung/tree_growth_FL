rm(list=ls())
lapply(c("ggplot2","mvnfast"), require, character.only = T)

# setwd("C:/Users/Alicia/Documents/GitHub/FL_Carbon/Longleaf Remeasurment")
envdata<-read.csv("LongleafEnvData1.csv", header=T, sep=",") 
load("longleafAgeTotals.rdata")
load("longleafAgeTotalsEnd.rdata")
load("longleafDIATotalsEnd.rdata")
load("longleafDIATotals.rdata")
load("longleafPlotStart.rdata")
load("longleafPlotEnd.rdata")

# setwd("C:/Users/Alicia/Documents/GitHub/tree_growth_FL")
source("generate_pars.R")

mylist<-list()

Diameter.all<-list()



# Plot 170 ----------------------------------------------
#set prior ranges
pmin <- c()
pmax <- c()
#   a1             b1              b2               b3           b4               b5            b6                b7                b8    
pmin[1] <- 0.8;pmin[2] <-  0.00;pmin[3] <-  0.00;pmin[4] <- 1;pmin[5] <-  0.01;pmin[6] <- -0.72;pmin[7] <- 0.00;pmin[8] <- -0.01;pmin[9] <- 0.00001
pmax[1] <- 2.8;pmax[2] <-  0.04;pmax[3] <-  0.05;pmax[4] <- 5;pmax[5] <-  1.00;pmax[6] <- -0.12;pmax[7] <- 0.4;pmax[8] <- -0.003;pmax[9] <- 0.00015

# growth<-(10^(1.872009 - 0.016069*CN_scale + 0.019731*CN_scale*aridity - 1.912994*aridity))*(0.5718692*age[i,j]^(-0.4281308))
# if (Diameter[i,j]>0) {M<- rbinom(1,1,(1.636e-01 +(-8.008e-03 *Diameter[i,j]*2.54)+(9.980e-05*(Diameter[i,j]*2.54^2)))/4)}

par.name <- c("a1","b1","b2","b3","b4","b5", "b6", "b7", "b8")
p_op <- c(1.872009, 0.016069, 0.019731, 1.912994, 0.5718692, -0.4281308, 1.636e-01, -8.008e-03, 9.980e-05)
no.simu <- 100000
d <- 6
J_last <- 100000
updated <- 0

#
p_rec <- matrix(NA,length(pmin), no.simu)
p_upgraded <- matrix(NA, length(pmin), no.simu)
J_keep <- rep(NA,no.simu)
J <- c()
DJ <- c()
J_new <- c()
sd <- 2.38/sqrt(length(pmin))
simu <- 0

# system.time({
for (d1 in 1:no.simu) {
  simu <- simu+1
  if (simu <= 2000) { #two steps; 1-default sampling; 2-sampling from covariances of each pars
    pnew <- generate.pars(p_op, pmin, pmax, d)
  } else {
    pnew <- generate.pars.cov(p_op, pmin, pmax, covars)
  }
  
  p_rec[,simu] <- pnew #save pnew
  
  #assign pars
  for (b in 1:length(par.name)) {
    assign(par.name[b], pnew[b])
  }
  
  
  
  a<-c(1, 3, 4, 6, 8, 9, 11:21, 23)
  
  
  for (s in a){
    
    # set stand age and density
    plot_density<-diameter.totals[[s]]
    observed.a<-envdata[s,6]
    CN <- envdata[s,16]
    aridity <- envdata[s,17]
    ages<-age.totals[[s]]
    predict.tasb<-matrix(nrow = 10, ncol = 1,0)
    predict.d<-matrix(nrow = 10, ncol = 1,0)
    
    
    # for (o in 1:10){
    age<-matrix(0, nrow=observed.a, ncol=length(plot_density)) # initialize the age matrix
    Diameter<-matrix(0, nrow=observed.a, length(plot_density)) # initialize the diameter matrix
    TASB<-matrix(0, nrow=observed.a, length(plot_density)) # initialize the total above-stump biomass matrix
    
    # initialize the diameter for the first year
    Diameter[1,]<-plot_density
    age[1,]<-ages
    
    for (j in 1:length(plot_density)){ # specify tree per hectare
      
      for (i in 2:observed.a){ # specify how long to run the simulation (years)
        
        age[i,j]<-age[i-1,j]+1
        
        
        # growth<-(10^(a1 - b1*envdata[s,16] + b2*envdata[s,16]*envdata[s,17] - b3*envdata[s,17]))*(b4*age[i,j]^(b5))
        # growth<-(10^(a1 - b1*CN + b2*CN*aridity - b3*aridity))*(b4*age[i,j]^(b5))
        growth<-(10^(a1 - b1*CN + b2*CN*aridity - b3*aridity))*(b4*age[i,j]^(b5))
        # growth <- 0.1
        #growth<-(10^(0.9333281 - 0.009259*CN + 0.011340*CN*aridity - 0.977477*aridity))*(0.5718692*age[i,j]^(-0.4281308))
        
        # define the mortality rate here
        # initialize as a numeric with only 1 possible value
        M <- numeric(length = 1)
        
        # Mortality based on diameter class
        # mort <- (b6 +(b7 *Diameter[i,j]*2.54)+(b8*(Diameter[i,j]*2.54^2)))/4;mort.binom <- rbinom(1,1,mort)
        # print(paste(mort, mort.binom, i))
        if (Diameter[i,j]>0) {
          M<- rbinom(1,1,(b6 +(b7 *Diameter[i,j]*2.54)+(b8*(Diameter[i,j]*2.54^2)))/4)
        }
        
        # else {M=0}
        # print(M)
        # if (Diameter[i,j]>0) { M<- rbinom(1,1,(2.109e-02 +(-2.662e-03 *Diameter[i-1,j])+(8.540e-05*(Diameter[i-1,j]^2))))}
        # else {M=0}
        # print(M)
        # if (Diameter[i,j]<=3.94){ M<- rbinom(1,1,(.162/8))}
        # else if (Diameter[i,j]>3.94 & Diameter[i,j]<=7.97){M<-rbinom(1,1,(.026/8))}
        # else if (Diameter[i,j]>7.97 & Diameter[i,j]<=11.81){M<-rbinom(1,1,(.006/8))}
        # else if (Diameter[i,j]>11.81 & Diameter[i,j]<=15.75){M<-rbinom(1,1,(.013/8))}
        # else if (Diameter[i,j]>15.75 & Diameter[i,j]<=19.69){M<-rbinom(1,1,(.024/8))}
        # else if (Diameter[i,j]>19.69 & Diameter[i,j]<=23.62){M<-rbinom(1,1,(.047/8))}
        # else if (Diameter[i,j]>23.62 & Diameter[i,j]<=27.56){M<-rbinom(1,1,(.060/8))}
        # else if (Diameter[i,j]>27.56){M<-rbinom(1,1,(0.129/8))}
        
        
        # Calculate the diameter for jth tree for the ith observed year
        Diameter[i,j]<-Diameter[i-1,j] + growth - M*(Diameter[i-1,j]+growth)
        
        #use diameter and age to calculate total aboveground biomass of the jth tree in the ith year
        TASB [i,j]<-(0.0725*((Diameter[i,j]*2.54)^2.5074))+(0.0016*((Diameter[i,j]*2.54)^3.0786))+(0.0214*((Diameter[i,j]*2.54)^2.0051))
        
        # If the tree dies, plant a new tree (age = 0)
        if (M==1){
          age[i,j]<-0
        }
      }
    }
    
    #   # save average modeled diameter
    #   predict.d[o,1]<-mean(Diameter[observed.a,])
    #   predict.tasb[o,1]<-sum(TASB[observed.a, ])*.5*(1/4047)
    # }
    
    # observed.d<-mean(diameter.totals[[s]])
    # observed.tasb<-sum(plot_data_end[[s]]$TASB)*.5*(1/4047)
    # modeled.d<-mean(predict.d)
    # modeled.tasb<-mean(predict.tasb)
    # sd.tasb<-sd(predict.tasb)
    # sd.dbh<-sd(predict.d)
    # hist(Diameter[observed.a,], main=paste("End simulation", s), xlab = "Diameter (in)")
    
    # set up dataframe to store simulated data
    # df<-cbind(observed.d, modeled.d, observed.a, length(plot_density), envdata[s,12], envdata[s,17], observed.tasb, modeled.tasb, 
    #           sd.tasb, sd.dbh)
    # 
    # df<-data.frame(df)
    # 
    # colnames(df)<-c("Observed_Diameter","Modeled_Diameter","Age", "Tree_Density", "Temperature",  "Aridity", "Observed_Biomass", 
    #                 "Modeled_Biomass", "sd.tasb", "sd.dbh")
    
    # data will be saved as list 1
    mylist[[s]] <- df
    Diameter.all[[s]]<-Diameter[observed.a,]
  }
  
  #calculate quantiles
  # for (z in a) {
  # simu.quan <- quantile(Diameter.all[[z]], prob = seq(0,1,0.03),na.rm = T)
  # obs.quan <- quantile(diameter.totals.end[[z]], prob = seq(0,1,0.03))
  J_new <- mapply(function(x, y){
    calc <- (sort(x) - sort(y))^2/(2*(0.6*y)^2)
    return(calc)
    }, Diameter.all, diameter.totals.end)
  # J[[z]] <- (sort(Diameter.all[[z]])) 
  #calculate J
  # J_new[[z]] <- sum((simu.quan - obs.quan)^2/(2*(0.6*obs.quan)^2)) #changed 0.6 to 2 from Sasha's recommendation

   # }
  J_new1 <-unlist(J_new)
  delta.J <- sum(J_new1) - J_last
  
  # print(sum(J_new1))
  if (min(1, exp(-delta.J)) > runif(1)) {
    p_op <- pnew
    J_last <- sum(J_new1)
    updated <- updated + 1
    J_keep[updated] <- sum(J_new1) 
    p_upgraded[,updated] <- p_op
    # plot(prob_denom.keep)
    if (updated %in% c(100*1:100)) {
      par(mfrow=c(1,9))
      par(mar=c(2,3,1,1))
      for (par.no in 1:9) {
        hist(p_upgraded[par.no,(updated/2):updated],xlim = c(pmin[par.no],pmax[par.no]), main = par.name[par.no], xlab =NA)
      }
      #   hist(obs.d, breaks.his, main = "Observed diameter", xlab="Diameter (in)")
      #   hist(simul.tree.d.last, breaks.his, main = "Simulated diameter", xlab="Diameter (in)")
      #   plot(sort(simul.tree.d.last), sort(obs.d), xlab = "observed DBH", ylab = "simulated DBH", xlim = c(4.5,16), ylim = c(4.5,16))
      #   abline(0,1)
    }
  }
  
  if (simu == 2000) {
    covars <- cov(t(p_rec[,1:simu]))
  }
  if (simu > 2000) {
    covars <- sd*cov(t(p_rec[,1:simu]))
  }
  print(paste("simu =", simu, "updated =", updated))
}







# save(simu.diameters, file="C:/Users/Alicia/Documents/GitHub/FL_Carbon/Longleaf Remeasurment/longleafDIATotals.simu2.4.rdata")

par(mfrow=c(1,1))

final_list <- do.call(rbind.data.frame, mylist)
row.names(final_list) <- c(a)
# final_list$Plot<- c(1,12,14:18)
# final_list$Plot<-as.factor(final_list$Plot)
sdev<-as.vector(final_list$sd.tasb)
sdev.dbh<-as.vector(final_list$sd.dbh)

model.1<-lm(data = final_list, log10(Observed_Biomass/Modeled_Biomass)~Tree_Density +Temperature + Aridity)
summary(model.1)
# plot(Temperature~Aridity, data= final_list[[1]])
# 
# res.age<-residuals(model.1)
# plot(res.age~df2$Age)
# plot(res.age~df2$Aridity)

model.2<-lm(data = final_list_slash, log10(Observed_Diameter/Modeled_Diameter)~Tree_Density)
summary(model.2)

model.3<-lm(data = final_list, Modeled_Biomass~Observed_Biomass)
summary(model.3)

plot(data = final_list, Observed_Diameter~Modeled_Diameter,  xlim = c(1,12), ylim = c(1,12), xlab="Modeled DBH (in)", ylab="Observed DBH (in)",
     main = "Before parameter correction", col.axis="#027368", col="#75BFBF", pch=16, type="p") + abline(0,1, col="#048ABF")
text(Observed_Diameter~Modeled_Diameter, labels=rownames(final_list),data=final_list, cex=0.9, font=2, pos=4)
arrows(final_list$Modeled_Diameter-sdev.dbh, final_list$Observed_Diameter, final_list$Modeled_Diameter+sdev.dbh, final_list$Observed_Diameter, length=0.05, angle=90, code=3)

plot(data = final_list, Observed_Biomass~Modeled_Biomass, xlim = c(0,6), ylim = c(0,6),col = "#75BFBF", xlab="Modeled", ylab="Observed", main ="AGB (kgC/m^2)",
     col.axis="#027368", pch=16, type="p") + abline(0,1, col="#048ABF")
text(Observed_Biomass~Modeled_Biomass, labels=rownames(final_list),data=final_list, cex=0.9, font=2, pos=4)
arrows(final_list$Modeled_Biomass-sdev, final_list$Observed_Biomass, final_list$Modeled_Biomass+sdev, final_list$Observed_Biomass, length=0.05, angle=90, code=3)
