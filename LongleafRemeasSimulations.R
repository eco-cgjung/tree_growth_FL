# Run longleaf pine simulations ----------------------------------------------

rm(list=ls())
source("generate_pars.R")
# setwd("C:/Users/Alicia/Documents/GitHub/FL_Carbon")
envdata<-read.csv("LongleafRemeasEnvData.csv", header=T, sep=",") 

# setwd("C:/Users/Alicia/Desktop/FL")
load("longleaf_remeas_start.Rdata")
load("longleaf_remeas_end.Rdata")

# envdata$PLT_CN<-sort(envdata$PLT_CN, decreasing = FALSE)
# 
# envdata$PLT_CN<-as.character(envdata$PLT_CN)


par(mfrow=c(2,5))

mylist<-list()


# Plot 170 ----------------------------------------------
diameter.totals<-list()
age.totals<-list()
diameter.totals.end<-list()
age.totals.end<-list()
plot_data_end<-list()
plot_data_start<- list()
Diameter.all <- list()

envdata$SOILGRIDS_C_AVG<-((envdata$SOC0_5_NAD/10)*(1/3))+((envdata$SOC5_15NAD/10)*(2/3))

envdata$SOILGRIDS_N_AVG<-((envdata$N0_5_NAD/100)*(1/3))+((envdata$N5_15_NAD/100)*(2/3))                 


envdata$SOILGRIDS_CN<-envdata$SOILGRIDS_C_AVG/envdata$SOILGRIDS_N_AVG 

envdata$SOILGRIDS_CN_SCALE<-(envdata$SOILGRIDS_CN*6.212)+24.634

a<-c(1, 3:23)


#set prior ranges
pmin <- c()
pmax <- c()
#   a             b1              b2               b3           b4               b5             
pmin[1] <- 0.8;pmin[2] <-  0.00;pmin[3] <-  0.00;pmin[4] <- 1;pmin[5] <-  0.01;pmin[6] <- -0.72
pmax[1] <- 2.8;pmax[2] <-  0.04;pmax[3] <-  0.05;pmax[4] <- 5;pmax[5] <-  1.00;pmax[6] <- -0.12

# growth<-(10^(1.872009 - 0.016069*CN_scale + 0.019731*CN_scale*aridity - 1.912994*aridity))*(0.5718692*age[i,j]^(-0.4281308))

par.name <- c("a1","b1","b2","b3","b4","b5")
p_op <- c(1.872009, 0.016069, 0.019731, 1.912994, 0.5718692, -0.4281308)
no.simu <- 100000
d <- 6
J_last <- 20000
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

for (d in 1:no.simu) {
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
  
  # system.time({
for (s in a){
  
  aridity<-envdata[s, 6]*0.0001
  
  temp<-envdata[s, 11]
  
  CN_scale<-envdata[s,15]
  
  plot_data_start[[s]]<-plots.start[[s]] %>%
    filter(STATUSCD=="1") %>%
    # mutate(TPA=ifelse(is.na(TPA_UNADJ), TPAGROW_UNADJ, TPA_UNADJ)) %>%
    mutate(TPA_total=sum(round(TPA_UNADJ))) %>%
    mutate(age=round((10^(-3.273492 + 0.02809908*CN_scale - 0.03450265*CN_scale*aridity + 3.34516*aridity))*(DIA^1.748652)))
  # mutate(age=round((10^(-3.018915 + 1.2741*aridity))*(DIA^2.307065)))  
  
  
  for (h in 1:length(plot_data_start$DIA)){
    diameter.totals[[s]]<-rep((plot_data_start[[s]]$DIA), round(plot_data_start[[s]]$TPA_UNADJ))
  age.totals[[s]]<-rep((plot_data_start[[s]]$age), round(plot_data_start[[s]]$TPA_UNADJ))
  }
  
  # hist(diameter.totals[[s]], main = paste("Start plot", s), xlab = "Diameter (in)")
  
  plot_data_end[[s]]<-plots.end[[s]] %>%
    filter(STATUSCD=="1") %>%
    # mutate(TPA=ifelse(is.na(TPA_UNADJ), TPAGROW_UNADJ, TPA_UNADJ)) %>%
    mutate(TPA_total=sum(round(TPA_UNADJ))) %>%
    mutate(age=round((10^(-3.273492 + 0.02809908*CN_scale - 0.03450265*CN_scale*aridity + 3.34516*aridity))*(DIA^1.748652))) %>%
    mutate(TASB=((0.0725*((DIA*2.54)^2.5074))+(0.0016*((DIA*2.54)^3.0786))+(0.0214*((DIA*2.54)^2.0051)))*(round(TPA_UNADJ)))
    # mutate(TASB=(0.041281*((DIA*2.54)^2.722214))*(round(TPA_UNADJ)))
  
  
  for (h in 1:length(plot_data_end$DIA)){
    diameter.totals.end[[s]]<-rep((plot_data_end[[s]]$DIA), round(plot_data_end[[s]]$TPA_UNADJ))
  age.totals.end[[s]]<-rep((plot_data_end[[s]]$age), round(plot_data_end[[s]]$TPA_UNADJ))
  }
  
  # hist(diameter.totals.end[[s]], main =paste("End plot", s), xlab = "Diameter (in)")
  
  
  
  
  
  # # set parameters for growth equation
  # aridity<-envdata[a, 7]*0.0001
  # 
  # temp<-envdata[a, 6]
  # 
  # temp2<-envdata[a, 6]^2
  # 
  # # C_AVG<-((envdata[1,11]/10)*(1/3))+((envdata[1,12]/10)*(2/3))
  # # 
  # # N_AVG<-((envdata[1,13]/100)*(1/3))+((envdata[1,14]/100)*(2/3))
  # # 
  # # CN<-C_AVG/N_AVG
  # # 
  # # CN_SCALE<-(CN*4.941)+29.777
  
  
  # set stand age and density
  plot_density<-diameter.totals[[s]]
  observed.a<-envdata[s,5]
  ages<-age.totals[[s]]
  predict.tasb<-matrix(nrow = 10, ncol = 1,0)
  predict.d<-matrix(nrow = 10, ncol = 1,0)
  
  for (o in 1:10){
    age<-matrix(0, nrow=observed.a, ncol=length(plot_density)) # initialize the age matrix
    Diameter<-matrix(0, nrow=observed.a, length(plot_density)) # initialize the diameter matrix
    TASB<-matrix(0, nrow=observed.a, length(plot_density)) # initialize the total above-stump biomass matrix
    
    # initialize the diameter for the first year
    Diameter[1,]<-plot_density
    age[1,]<-ages
    
    for (j in 1:length(plot_density)){ # specify tree per hectare
      
      for (i in 2:observed.a){ # specify how long to run the simulation (years)
        
        age[i,j]<-age[i-1,j]+1
        

        growth<-(10^(a1 - b1*CN_scale + b2*CN_scale*aridity - b3*aridity))*(b4*age[i,j]^(b5))
        # growth<-(10^(1.872009 - 0.016069*CN_scale + 0.019731*CN_scale*aridity - 1.912994*aridity))*(0.5718692*age[i,j]^(-0.4281308))
        #growth<-(10^(0.9333281 - 0.009259*CN + 0.011340*CN*aridity - 0.977477*aridity))*(0.5718692*age[i,j]^(-0.4281308))
        
        # define the mortality rate here
        # initialize as a numeric with only 1 possible value
        M <- numeric(length = 1)
        
        # Mortality based on diameter class
        if (Diameter[i,j]<=3.94){ M<- rbinom(1,1,(.162/8))}
        else if (Diameter[i,j]>3.94 & Diameter[i,j]<=7.97){M<-rbinom(1,1,(.026/8))}
        else if (Diameter[i,j]>7.97 & Diameter[i,j]<=11.81){M<-rbinom(1,1,(.006/8))}
        else if (Diameter[i,j]>11.81 & Diameter[i,j]<=15.75){M<-rbinom(1,1,(.013/8))}
        else if (Diameter[i,j]>15.75 & Diameter[i,j]<=19.69){M<-rbinom(1,1,(.024/8))}
        else if (Diameter[i,j]>19.69 & Diameter[i,j]<=23.62){M<-rbinom(1,1,(.047/8))}
        else if (Diameter[i,j]>23.62 & Diameter[i,j]<=27.56){M<-rbinom(1,1,(.060/8))}
        else if (Diameter[i,j]>27.56){M<-rbinom(1,1,(0.129/8))}
        
        # Calculate the diameter for jth tree for the ith observed year
        Diameter[i,j]<-Diameter[i-1,j] + growth - M*(Diameter[i-1,j]+growth)
        
        #use diameter and age to calculate total aboveground biomass of the jth tree in the ith year
        TASB [i,j]<-(0.0725*((Diameter[i,j]*2.54)^2.5074))+(0.0016*((Diameter[i,j]*2.54)^3.0786))+(0.0214*((Diameter[i,j]*2.54)^2.0051))
        
        # If the tree dies, plant a new tree (age = 0)
        if (M==1){age[i,j]<-0}
      }
    }
    # save average modeled diameter
    predict.d[o,1]<-mean(Diameter[observed.a,])
    predict.tasb[o,1]<-sum(TASB[observed.a, ])*.5*(1/4047)
  }
  
  observed.d<-mean(diameter.totals[[s]])
  observed.tasb<-sum(plot_data_end[[s]]$TASB)*.5*(1/4047)
  modeled.d<-mean(predict.d)
  modeled.tasb<-mean(predict.tasb)
  sd.tasb<-sd(predict.tasb)
  sd.dbh<-sd(predict.d)
  # hist(Diameter[observed.a,], main=paste("End simulation", s), xlab = "Diameter (in)")
  
  # set up dataframe to store simulated data
  df<-cbind(observed.d, modeled.d, observed.a, length(plot_density), temp, aridity, observed.tasb, modeled.tasb, 
            sd.tasb, sd.dbh)
  
  df<-data.frame(df)
  
  colnames(df)<-c("Observed_Diameter","Modeled_Diameter","Age", "Tree_Density", "Temperature",  "Aridity", "Observed_Biomass", 
                  "Modeled_Biomass", "sd.tasb", "sd.dbh")
  
  # data will be saved as list 1
  mylist[[s]] <- df

  Diameter.all[[s]] <- Diameter[observed.a,]  
}
  # })
  # diameter.totals.end
  # 
  # Diameter.all
  # 
  
  #1
  # diameter.totals.end[[1]]
  # Diameter.all[[1]]
  # 
  # diameter.totals.end[[4]]
  # Diameter.all[[4]]
  # 
  # breaks.his <- seq(min(c(diameter.totals.end[[4]],Diameter.all[[4]])),max(c(diameter.totals.end[[4]],Diameter.all[[4]])), length.out = 10) 
  # obs.count <- hist(diameter.totals.end[[4]], breaks.his,plot = F)$counts
  # simul.count <- hist(Diameter.all[[4]], breaks.his, plot = F)$counts
  # 
  for (z in a) {
    breaks.his <- seq(min(c(diameter.totals.end[[z]],Diameter.all[[z]])),max(c(diameter.totals.end[[z]],Diameter.all[[z]])), length.out = 10) 
    obs.count <- hist(diameter.totals.end[[z]], breaks.his,plot = F)$counts
    simul.count <- hist(Diameter.all[[z]], breaks.his, plot = F)$counts
    J <- (simul.count - obs.count)^2
    DJ <- 2*(2*obs.count)^2
    DJ[DJ == 0] <- 1
    J_new[[z]] <- J/DJ
  }
  
  J_new1 <- unlist(J_new)
  delta.J <- sum(J_new1) - J_last
  
  if (min(1, exp(-delta.J)) > runif(1)) {
    p_op <- pnew
    J_last <- sum(J_new1)
    updated <- updated + 1
    J_keep[updated] <- sum(J_new1) 
    p_upgraded[,updated] <- p_op
    # plot(prob_denom.keep)
    # if (updated %in% c(100*1:100)) {
    #   par(mfrow=c(2,7))
    #   par(mar=c(5,5,1,1))
    #   for (i in 1:7) {
    #     hist(p_upgraded[i,(updated/2):updated],xlim = c(pmin[i],pmax[i]), main = par.name[i], xlab =NA)
    #   }
    #   hist(obs.d, breaks.his, main = "Observed diameter", xlab="Diameter (in)")
    #   hist(simul.tree.d.last, breaks.his, main = "Simulated diameter", xlab="Diameter (in)")
    #   plot(sort(simul.tree.d.last), sort(obs.d), xlab = "observed DBH", ylab = "simulated DBH", xlim = c(4.5,16), ylim = c(4.5,16))
    #   abline(0,1)
    # }
  }
  
  if (simu == 2000) {
    covars <- cov(t(p_rec[,1:simu]))
  }
  if (simu > 2000) {
    covars <- sd*cov(t(p_rec[,1:simu]))
  }
  print(paste("simu =", simu, "updated =", updated))
}
  ##
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

