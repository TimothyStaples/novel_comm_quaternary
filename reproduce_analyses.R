# ######################################################################### ####
# Emergence patterns of locally novel plant communities driven by past      ####
# climate change and modern anthropogenic impacts                           ####
# Author:    Timothy L Staples                                              ####
# Collaborators: John Pandolfi                                              #### 
#                Wolfgang Kiessling                                         ####
# ######################################################################### ####
# Script purpose ####

# This script imports the processed novel community detection framework data and
# runs all main text and most supplementary analyses, producing all tables and figures.
# The full analysis pathway requires large data files that are not provided by this
# repository, but can be recreated using the "data_acquisition_processing.R" and 
# "full_analysis.R" scripts.

# Global attributes & working directories ####

rm(list=ls())

setwd("/Users/uqtstapl/Dropbox/Tim/Post-doc/Research projects/ecolLettNeotoma")

# Packages & functions ####

# source functions from 'functions' sub-folder
sapply(paste0("./functions/", list.files("./functions")), source)

package.loader(c("sp", "maptools", "rworldmap", "lmtest",
                 "vegan", "DHARMa", "nlme",
                 "lme4", "fossil", "multcomp", "vegclust",
                 "betareg", "mgcv", "gamm4", "shape", "merTools",
                 "MuMIn", "rgdal", "rgeos", "raster"))

# A little custom function to add dates to output files
date.wrap <- function(string, ext){
  paste0(string, " ", Sys.Date(), ext)
}

colour.mat <- read.csv("./raw.datafiles/colour.mat.csv",
                       stringsAsFactors = FALSE)

# IMPORT NOVELTY ####

plant.with.genus <- read.csv(unz("./raw.datafiles/processed_genus_records.zip",
                                 "processed_genus_records.csv"))

plant.with.genus$site = plant.with.genus$site.id
colnames(plant.with.genus)[colnames(plant.with.genus)=="site.1"] = "site"
site.df <- read.csv("./outputs/siteDf.csv")
colnames(site.df)[colnames(site.df) == "site.id"] = "site"

all.novel <- readRDS("./outputs/all neotoma novelty (sub-sampled).rds")

# IMPORT CLIMATE DATA ####
#           Convert temp records onto same scale ####

mod.env.data <- read.csv("raw.datafiles/Marcott2013_5x5deg.csv")

mod.env.data1 <- read.csv("./raw.datafiles/temp12k_allmethods_percentiles.csv")

cor.factor <- mean(mod.env.data$temp[mod.env.data1$age <= 11500 &
                                       mod.env.data1$age >= 6500])

old.env.data <- read.csv("raw.datafiles/Shakun2012_retreat_temp.csv")
old.env.data$temp = old.env.data$temp + cor.factor
old.env.data$age = old.env.data$age * 1000

plot(mod.env.data$temp ~ mod.env.data$age, type='l', xlim=c(0,20000), ylim=c(-4,0.6))
lines(mod.env.data1$global_median ~ mod.env.data1$ages, col="Red")
lines(old.env.data$temp ~ old.env.data$age, col="blue")

pdf("./plots/shakun marcott curve compare.pdf", height=3, width=5)
par(mar=c(3,2,1,1), ps=8)
plot(x=NULL, y=NULL, type="n", xlim=c(-100,22000), ylim=c(-4,1))
lines(mod.env.data$temp ~ mod.env.data$age, lwd=1)
lines(mod.env.data1$global_median ~ mod.env.data1$ages, col="Red")
lines(old.env.data$temp ~ old.env.data$age, lwd=1, col="blue")
dev.off()

# Interpolate temp for each 200 year bin from reconstructions
mod.env.sub <- mod.env.data1[,c("ages","global_median")]
colnames(mod.env.sub) <- c("age", "temp")
comb.env.data <- rbind(old.env.data[old.env.data$age>max(mod.env.data1$ages),1:2], 
                       mod.env.sub)

temp.gam <- gam(temp ~ s(age, bs="cr", k=200), data=comb.env.data)
summary(temp.gam)

pred.df <- data.frame(age = seq(-100, 22000, 50))
temp.df <- cbind(pred.df,
                 as.data.frame(predict(temp.gam, newdata=pred.df, se.fit=TRUE)))

# pretty good fit
pdf("./plots/temp curve compare.pdf", height=3, width=5)
par(mar=c(3,2,1,1), ps=8)
plot(comb.env.data$temp ~ comb.env.data$age, pch=16, cex=1, type="n", xlim=c(-100,38000), ylim=c(-8,1))
#rect(xleft=-100, xright=100, ybottom=-10, ytop=10, col="grey90")
lines(temp.df$fit ~ pred.df$age, col="red", lwd=2)
#abline(v=seq(0,20000,1000))
lines(y=snyder.env.data$temp, x=snyder.env.data$age*1000, lwd=2, col="black")
dev.off()

recent.env.dat$age <- 1950 - recent.env.dat$Year
recent.env.dat$temp <- recent.env.dat$NTREND2015spl20
str(recent.env.dat)
recent.bin <- bin.env.data(env.data=recent.env.dat,
                           bin.width=200,
                           lims=c(-100,1300),
                           env.var = "temp")

existing.bin <- bin.env.data(env.data=temp.df,
                             bin.width=200,
                             lims=c(-100,1300),
                             env.var = "fit")

pdf("./plots/modern temp curve compare.pdf", height=3, width=5)
par(mar=c(3,2,1,1), ps=8)
plot(existing.bin$env ~ existing.bin$bin, pch=16, cex=1, type="n", xlim=c(-100,1000), ylim=c(-1,1.5))
segments(x0=existing.bin$bin +25, 
         x1=existing.bin$bin +25,
         y0 = existing.bin$env + 1.96 * existing.bin$env.se,
         y1 = existing.bin$env - 1.96 * existing.bin$env.se)
points(y=existing.bin$env, x=existing.bin$bin +25, pch=21, bg="red")

segments(x0=recent.bin$bin -25, 
         x1=recent.bin$bin -25,
         y0 = recent.bin$env + 1.96 * recent.bin$env.se,
         y1 = recent.bin$env - 1.96 * recent.bin$env.se)
points(y=recent.bin$env, x=recent.bin$bin-25, pch=21, bg="grey")
dev.off()

# modern comparison
temp.data <-  bin.env.data(env.data=temp.df,
                           bin.width=200,
                           lims=c(-100,23000),
                           env.var = "fit")
colnames(temp.data) <- gsub("env", "temp", colnames(temp.data))

write.csv(temp.data, "./outputs/interpolated temp data.csv")

#           local climate conditions ####

# This is a arrayed version of the CHELSA_TraCE21k data which was imported
# from raw files and processed as per the "full_analyses.R" script
# in the IMPORT CLIMATE DATA / local climate conditions section
localTempList <- readRDS("./outputs/localTempList.rds")

length(localTempList)

times <- sapply(localTempList, function(x){x$time[1]})
times <- (as.numeric(times)) * 100
times <- -1* (times - 1950)

localTMat <- sapply(localTempList, function(x){x$localTemp})
rownames(localTMat) <- localTempList[[1]]$site
colnames(localTMat) <- times

localTMat <- localTMat[,order(as.numeric(colnames(localTMat)))]
localTMat[1:10,1:10]

# okay now we need to do what we did above, which is model temperature as a
# function of time, interpolate at 50 year intervals for our bins (as the CHELSA
# bins are not quite the same) and standardize to a d1950 mean.

localTPred <- t(apply(localTMat, 1, function(x){
  
  temp <- as.numeric(colnames(localTMat))
  
  tGam <- gam(x ~ s(temp))
  tPred <- predict(tGam, newdata=data.frame(temp = temp.df$age)) 
  tPred - mean(tPred[2:3])
}))
rownames(localTPred) <- rownames(localTMat)
colnames(localTPred) <- temp.df$age

# ANALYSES ####
# Novelty through time ####
#           Model ####

all.nprob.models <- novel.prob.models(novel.list = all.novel,
                                      site.df = site.df,
                                      time.k = 40,
                                      test.model=TRUE,
                                      name = "all",
                                      time.age.limits = c(0,25000),
                                      factor.age.limits = c(0,1000),
                                      sauto.n = NA,
                                      sauto.iter = 999)

novel.prob.plot(prob.model.list = all.nprob.models, 
                env.data = temp.data,
                mod.env.data = mod.env.data,
                ylims=c(0,0.048), 
                regylims=c(0,0.075),
                name = "all",
                time.age.limits = c(0,25000),
                factor.age.limits = c(0,1000),
                group.letters=c("A","AB","B","B","B","B"))

saveRDS(all.nprob.models,
        date.wrap("./outputs/novel probability models (all)", ".rds"))

#           Sub-sample time-series to examine sampling bias ####
  
# when we look at sampling, the most abundant time points are those that
# are at the height of the post-ice age novelty peak (~10-12K years ago).
# We need to test that our novelty relationship is not due to simply having
# more time-series in more places.

# I'm going to do this by writing a function to sub-sample our time-series data
# to equal numbers of observations at each point of the time-series, much like
# rarefaction sub-sampling. We then run the reduced model, extract the predicted 
# probability over time curve, and then do it a lot of times. From this we get a 
# distribution of sub-sampled curves, which we can match up with the overall 
# predicted one from the all data model.

# sample function
sub.novels <- subsample.novel(novel.list=all.novel,
                              novel.models = all.nprob.models,
                              time.age.limits=c(1200,25000),
                              factor.age.limits = c(0,1000),
                              full.preds = all.nprob.models$time.pred,
                              iter=999,
                              plot.limits=c(0,0.15))

saveRDS(sub.novels, 
        date.wrap("./outputs/subsampling novel over times results (all)", ".rds"))

# Environmental correlation / projection ####
#             Data prep ####
# run models on same subset of data, for a number of different time lags used to calculate
# temperature change
env.lags <- expand.grid(c(1:25),
                        c(1:25))

lagList <- lapply(c(1:25), function(n){
  calcLag(novel.list = all.novel,
          env.data = temp.df,
          env.var = "fit",
          local.env.data = localTPred,
          global.lag=n,
          local.lag=n)
})

comb.df <- do.call("rbind", all.novel$novel)
comb.df$site.bin <- paste0(comb.df$site, ":", comb.df$bins)
comb.env <- do.call("cbind", lapply(lagList, function(x){
  x$site.bin <- paste0(x$site, ":", x$bin)
  temp <- x[match(comb.df$site.bin,
                  x$site.bin),3:4]
  
  return(temp)
}))
comb.filled <- complete.cases(cbind(comb.env, comb.df$novel))
comb.env <- cbind(comb.df, comb.env)[comb.filled,]

lagCors <- sapply(lagList, function(x){
  cor(x[,3:4], use="complete.obs")[2]
})
summary(lagCors)

library(usdm)
lagVifs <- sapply(lagList, function(x){
  vif(x[,3:4])[1,2]
})
summary(lagVifs)

# try log-transforming absolute temperature change
logLags <- sapply(comb.env[,grepl("Lag", colnames(comb.env))],
                  function(x){log(abs(x))})

comb.envLog <- comb.env
comb.envLog[,grepl("Lag", colnames(comb.envLog))] = logLags

# subset strongly human-influenced period
comb.envLog <- comb.envLog[as.numeric(comb.envLog$bins) > 600,]

#             Global model ####

env.model.list <- lapply(1:nrow(env.lags), function(n){

  print(paste0("global = ", env.lags[n,1], ": local = ", env.lags[n,2]))
  a <- novel.by.env.local(comb.df = comb.envLog,
                                      novel.list=all.novel,
                                      global.lag=env.lags[n,1],
                                      local.lag=env.lags[n,2])

  return(AIC(a))

})

targetLags <- env.lags[which.min(unlist(env.model.list)),]
#targetLags <- c(14, 25) # Upd 4/1/22

targetM <-  novel.by.env.local(comb.df = comb.envLog,
                               novel.list=all.novel,
                               global.lag=targetLags[1],
                               local.lag=targetLags[2])

envDiag <- modelDiagTests(model = targetM,
                          time=as.numeric(targetM$data$bins),
                          data=targetM$data,
                          spat.iter=999,
                          spat.sub.size=1)

#             plot pred ####

local.pred <- data.frame(local.lag = seq(min(exp(targetM$data$localLag24)),
                                         max(exp(targetM$data$localLag24)), len=200))
local.pred$local.lag.scale <- (log(local.pred$local.lag) - mean(targetM$data$localLag24)) / sd(targetM$data$localLag24)

global.pred <- data.frame(global.lag = seq(min(exp(targetM$data$globalLag9)),
                                                 max(exp(targetM$data$globalLag9)), len=200))
global.pred$global.lag.scale <- (log(global.pred$global.lag) - mean(targetM$data$globalLag9)) / sd(targetM$data$globalLag9)

grid.pred.df <- expand.grid(local.lag.scale=local.pred$local.lag.scale,
                         global.lag.scale=global.pred$global.lag.scale,
                         bin.lag.scale = 0,
                         bin.n.scale = 0,
                         tsLength.scale = 0,
                         tsRichness.scale = 0,
                         elev.scale = 0)

grid.pred <- predict(targetM, newdata=grid.pred.df)

#             Global plot ####

pdf("./plots/tempCorr (all).pdf", height=8, width=9, useDingbats = FALSE)

split.screen(rbind(c(0.1,0.475,0.47,0.88),
                   c(0.15,0.4,0.905, 0.93),
                   c(0.575,0.95,0.47,0.88),
                   c(0.625,0.9,0.905,0.93),
                   c(0.1,0.95,0.07,0.4)))

screen(1)
par(mar=c(0,0,0,0), ps=10, mgp=c(3,0.5,0), las=1, tcl=-0.25)

gridX <- seq(0,1, len=nrow(env.AIC)+2)[-c(1,nrow(env.AIC)+2)]
gridC <- expand.grid(gridX,
                     gridX)
gridW <- 0.5*diff(gridX[1:2])

plot(x=NULL, y=NULL, xlim=c(min(gridX) - gridW,
                            max(gridX) + gridW), 
     ylim=c(min(gridX) - gridW,
            max(gridX) + gridW), xaxs="i", yaxs="i",
     axes=FALSE, xlab="", ylab="")
box()

library(viridisLite)
env.AIC.plot <- (env.AIC - env.AIC[1,1])
#env.AIC.plot <- env.AIC.plot[nrow(env.AIC.plot):1,]
env.AIC.plot <- t(env.AIC.plot)
AICsc <- (env.AIC.plot - min(env.AIC.plot)) / (max(env.AIC.plot) - min(env.AIC.plot))
 
rect(xleft= gridC[,1] - gridW,
     xright = gridC[,1] + gridW,
     ybottom = gridC[,2] - gridW,
     ytop = gridC[,2] + gridW,
     col=rgb(colorRamp(viridis(10, option="F"))(AICsc)/255), lwd=0.5)

targetN <- unlist((targetLags[1]-1)*25 + targetLags[2])
rect(xleft= gridC[targetN,1] - gridW,
     xright = gridC[targetN,1] + gridW,
     ybottom = gridC[targetN,2] - gridW,
     ytop = gridC[targetN,2] + gridW, border="white", lwd=2)

axis(side=1, at=gridX, labels=NA, tcl=-0.125)
axis(side=1, at=gridX[seq(5,25,5)], mgp=c(3,0.2,0),
     labels=format(seq(1000,5000,1000), big.mark=","), tcl=-0.25)
mtext(side=1, line=1.25, text="Local temperature lag (years)")

axis(side=2, at=gridX, labels=NA, tcl=-0.125)
axis(side=2, at=gridX[seq(5,25,5)], mgp=c(3,0.2,0),
     labels=NA, tcl=-0.25)
mtext(side=2, line=2.5, text="Global temperature lag (years)", las=0)

par(xpd=NA)
text(y=gridX[seq(5,25,5)], x=relative.axis.point(-0.035, "x"),
    labels=format(seq(1000,5000,1000), big.mark=","), srt=45, adj=1)
par(xpd=FALSE)
mtext(side=3, at=par("usr")[1], text="(A)", font=2, adj=0)
box()
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=10, mgp=c(3,0.5,0), las=1, tcl=-0.25)
daic <- as.vector(env.AIC.plot) * -1
plot(x=NULL, y=NULL, xlim=range(daic), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="")
image(x=seq(min(daic), max(daic), len=200),
      y=c(0,1),
      z=matrix(seq(0, max(daic), len=200), ncol=1),
      col=rgb(colorRamp(viridis(10, option="F"))(rev(seq(0, 1, len=200)))/255),
      add=TRUE, useRaster=TRUE)
axis(side=3, at=seq(0,120,30), labels=seq(0,120,30)*-1, mgp=c(3,0.2,0))
mtext(side=3, line=1.25, text = expression(Delta*"AIC"))

box()
close.screen(2)

screen(3)
par(mar=c(0,0,0,0), ps=10, mgp=c(3,0.5,0), las=1, tcl=-0.25)

image(x=local.pred$local.lag,
      y=global.pred$global.lag,
      z=matrix(plogis(grid.pred), nrow=nrow(global.pred)),
      col=rev(colorRampPalette(rev(c("white", "#b3cedc", "#8a87c7", 
                             "#80388e", "#4d193c", "black")), bias=1.25)(200)),
      useRaster=TRUE, axes=FALSE, xlab="", ylab="")
contour(x=local.pred$local.lag,
      y=global.pred$global.lag,
      z=matrix(plogis(grid.pred), nrow=nrow(global.pred)), add=TRUE,
      levels=seq(0.01,0.06,0.01))
contour(x=local.pred$local.lag,
        y=global.pred$global.lag,
        z=matrix(plogis(grid.pred), nrow=nrow(global.pred)), add=TRUE,
        levels=seq(0.06,0.07,0.01), col="white")

axis(side=1, mgp=c(3,0.1,0))
axis(side=2)
mtext(side=1, line=1.25, text = expression("Local temperature change (5,000 year lag ("*Delta*degree*"C))"))
mtext(side=2, line=2, text = expression("Global temperature change (3,000 year lag ("*Delta*degree*"C))"), las=0)
mtext(side=3, at=par("usr")[1], text="(B)", font=2, adj=0)
box()
close.screen(3)

screen(4)
par(mar=c(0,0,0,0), ps=10, mgp=c(3,0.5,0), las=1, tcl=-0.25)
novelrange <- range(plogis(grid.pred))
plot(x=NULL, y=NULL, xlim=c(0,max(novelrange)), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="")
image(x=seq(0, max(novelrange), len=200),
      y=c(0,1),
      z=matrix(seq(0, max(novelrange), len=200), ncol=1),
      col=rev(colorRampPalette(rev(c("white", "#b3cedc", "#8a87c7", 
                                     "#80388e", "#4d193c", "black")), bias=1.25)(200)),
      add=TRUE, useRaster=TRUE)
axis(side=3, mgp=c(3,0.2,0))
mtext(side=3, line=1.25, text = "Novel community emergence probability")
box()
close.screen(4)

screen(5)
par(mar=c(0,0,0,0), ps=10, mgp=c(3,0.5,0), las=1, tcl=-0.25)
resM <- residuals(targetM)
res.df <- data.frame(res = residuals(targetM),
                novel = as.factor(targetM$data$novel),
                bins = targetM$data$bins)

resGam <- gam(res ~ -1 + s(bins) + s(novel, bs="re"), data=res.df)

pred.df <- data.frame(novel = rep(c("FALSE", "TRUE"),each=length(unique(res.df$bins))),
                      bins = rep(sort(unique(res.df$bins)), 2))
pred.df <- cbind(pred.df, as.data.frame(predict(resGam, newdata=pred.df, se.fit=TRUE)))

pred.df <- data.frame(bins = sort(unique(res.df$bins)),
                      novel= "a")
pred.df <- cbind(pred.df, as.data.frame(predict(resGam, newdata=pred.df, se.fit=TRUE, exclude="s(novel)")))

plot(y=pred.df$fit, x=pred.df$bins, type="n",
     xlim=c(17000,600), ylim=c(-0.42,0.15), xlab="", ylab="", xaxt="n")

axis(side=1, at=seq(5000,15000,5000), labels=format(seq(5000,15000,5000), big.mark=","),
     mgp=c(3,0.2,0))
axis(side=1, at=0,mgp=c(3,0.2,0))
axis(side=1, at=seq(0,17000,1000), labels=NA, tcl=-0.125)
mtext(side=1, line=1.5, text=expression("Years before present"))
mtext(side=2, line=2.25, text="Centered model residuals", las=0)

res.df$resCor <- res.df$res - resGam$coefficients[grepl("novel", names(resGam$coefficients))][res.df$novel]

# 95% predictive intervals
PI <- do.call("rbind", tapply(res.df$resCor, res.df$bins, function(x){quantile(x, probs=c(0.025,0.975))}))
polygon(y=c(PI[,1], rev(PI[,2])),
            x=c(pred.df$bins, rev(pred.df$bins)), border=NA, col="grey75")

quarts <- do.call("rbind", tapply(res.df$resCor, res.df$bins, function(x){quantile(x, probs=c(0.33,0.66))}))

polygon(y=c(quarts[,1], rev(quarts[,2])),
        x=c(pred.df$bins, rev(pred.df$bins)), border=NA, col="grey50")

resMean <- tapply(res.df$resCor, res.df$bins, mean)
lines(y=resMean, x=pred.df$bins, type="l", lwd=2)
abline(h=0, lty="31")

text(x=12500, y=-0.22, labels=expression(bold(""%+-%"2"*sigma)), col="grey50")
text(x=11000, y=-0.085, labels=expression(bold(""%+-%"1"*sigma)), col="grey20")
text(x=600, y=resMean[1], labels=expression(bold(mu)), col="black", pos=4, offset=0.25)

text(x=relative.axis.point(0.005, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(C)", font=2, adj=0)

close.screen(5)

close.screen(all.screens=TRUE)
dev.off()

#             Regional model ####

# so we mix combinations of lags for each region to see whether they follow
# the same lag, only in steps of 5 because otherwise it's 25^4 which is an
# impossible number of model combinations

env.lags.reg <- expand.grid(seq(5,25,5),
                            seq(5,25,5),
                            seq(5,25,5),
                            seq(5,25,5))

env.model.list.region <- lapply(1:nrow(env.lags.reg), function(n){
  
  print(paste0("E global = ", env.lags.reg[n,1], ": E local = ", env.lags.reg[n,2],
        ": NA global = ", env.lags.reg[n,3], ": NA local = ", env.lags.reg[n,4]))
  a <- novel.by.env.local.region(comb.df = comb.envLog,
                          novel.list=all.novel,
                          global.lagE=env.lags.reg[n,1],
                          local.lagE=env.lags.reg[n,2],
                          global.lagNA = env.lags.reg[n,3],
                          local.lagNA = env.lags.reg[n,4])
  a$data <- NULL
  return(AIC(a))
  
})

env.AIC <- sapply(env.model.list.region, function(x){x
})
round(env.AIC - env.AIC[1], 1)

# what's the AIC improvement for increasing along each axis?
AIC.array <- array(unlist(env.model.list.region), dim=c(5,5,5,5),
                   dimnames=list(paste0("Ge", seq(5,25,5)),
                                 paste0("Le", seq(5,25,5)),
                                 paste0("Gna", seq(5,25,5)),
                                 paste0("Lna", seq(5,25,5))))

GeAIC <- apply(AIC.array, c(2,3,4), function(x){
  x[-1] - x[1]
})
apply(GeAIC, 1, quantile, probs=c(0.025,0.975))
apply(GeAIC, 1, mean)

LeAIC <- apply(AIC.array, c(1,3,4), function(x){
  x[-1] - x[1]
})
apply(LeAIC, 1, quantile, probs=c(0.025,0.975))
apply(LeAIC, 1, mean)

GnaAIC <- apply(AIC.array, c(1,2,4), function(x){
  x[-1] - x[1]
})
apply(GnaAIC, 1, quantile, probs=c(0.025,0.975))
apply(GnaAIC, 1, mean)

LnaAIC <- apply(AIC.array, c(1,2,3), function(x){
  x[-1] - x[1]
})
apply(LnaAIC, 1, quantile, probs=c(0.025,0.975))
apply(LnaAIC, 1, mean)

# now we know the best model in units of five, we can zoom in by going 2 lags
# below and 2 lags above to try and find the overall best model

# targetLags = c(15, 25, 10, 25) # Upd 4/2/21 
targetLags <- unlist(env.lags.reg[which.min(AIC.array),])

targetMReg <- novel.by.env.local.region(comb.df = comb.envLog,
                                    novel.list=all.novel,
                                    global.lagE=targetLags[1],
                                    local.lagE=targetLags[2],
                                    global.lagNA = targetLags[3],
                                    local.lagNA = targetLags[4])

targetRegData <- targetMReg$data[,!colnames(targetMReg$data) %in% c("lat.y", "long.y")]
colnames(targetRegData)[colnames(targetRegData) %in% c("lat.x", "long.x")] = c("lat", "long")

envRegDiag <- modelDiagTests(model = targetMReg,
                          time=as.numeric(targetRegData$bins),
                          data=targetRegData,
                          spat.iter=999,
                          spat.sub.size=1)

#             Regional plot ####

Egtar <- targetMReg$data[targetMReg$data$REGION=="Europe",
                         paste0("globalLag", targetLags[1])]

Eltar <- targetMReg$data[targetMReg$data$REGION=="Europe",
                         paste0("localLag", targetLags[2])]

global.lagE = seq(min(exp(Egtar)), max(exp(Egtar)), len=200)
local.lagE = seq(min(exp(Eltar)), max(exp(Eltar)), len=200)
local.predE <- expand.grid(global.lag = global.lagE,
                           local.lag=local.lagE,
                          bin.lag.scale = 0, bin.n.scale = 0,
                          tsLength.scale = 0, tsRichness.scale = 0,
                          elev.scale = 0, REGION="Europe")
local.predE$local.lag.scale <- (log(local.predE$local.lag) - mean(targetMReg$data$modelLocal)) / sd(targetMReg$data$modelLocal)
local.predE$global.lag.scale <- (log(local.predE$global.lag) - mean(targetMReg$data$modelGlobal)) / sd(targetMReg$data$modelGlobal)

NAgtar <- targetMReg$data[targetMReg$data$REGION=="North America",
                         paste0("globalLag", targetLags[3])]
NAltar <- targetMReg$data[targetMReg$data$REGION=="North America",
                         paste0("localLag", targetLags[4])]

global.lagNA = seq(min(exp(NAgtar)), max(exp(NAgtar)), len=200)
local.lagNA = seq(min(exp(NAltar)), max(exp(NAltar)), len=200)

local.predNA <- expand.grid(global.lag = global.lagNA,
                          local.lag = local.lagNA,
                          bin.lag.scale = 0, bin.n.scale = 0,
                          tsLength.scale = 0, tsRichness.scale = 0,
                          elev.scale = 0, REGION="North America")
local.predNA$local.lag.scale <- (log(local.predNA$local.lag) - mean(targetMReg$data$modelLocal)) / sd(targetMReg$data$modelLocal)
local.predNA$global.lag.scale <- (log(local.predNA$global.lag) - mean(targetMReg$data$modelGlobal)) / sd(targetMReg$data$modelGlobal)

grid.predE <- predict(targetMReg, newdata=local.predE)
grid.predNA <- predict(targetMReg, newdata=local.predNA)

pdf("./plots/tempCorrRegion.pdf", height=8, width=9, useDingbats = FALSE)

split.screen(rbind(c(0.1,0.475,0.47,0.88),
                   c(0.15,0.4,0.905, 0.93),
                   c(0.575,0.95,0.47,0.88),
                   c(0.35,0.65,0.905,0.93),
                   c(0.1,0.95,0.07,0.4)))

screen(1)
par(mar=c(0,0,0,0), ps=10, mgp=c(3,0.5,0), las=1, tcl=-0.25)

image(x=local.lagE,
      y=global.lagE,
      z=matrix(plogis(grid.predE), nrow=length(local.lagE)),
      col=rev(colorRampPalette(rev(c("white", "#b3cedc", "#8a87c7", 
                                     "#80388e", "#4d193c", "black")), bias=1.25)(200)),
      zlim=c(0,0.06),
      useRaster=TRUE, axes=FALSE, xlab="", ylab="")

contour(x=local.lagE,
        y=global.lagE,
        z=matrix(plogis(grid.predE), nrow=length(local.lagE)), add=TRUE,
        levels=seq(0.01,0.08,0.01))

axis(side=1, mgp=c(3,0.1,0))
axis(side=2)
mtext(side=1, line=1.25, text = expression("Local temperature change (5,000 year lag ("*Delta*degree*"C))"))
mtext(side=2, line=2, text = expression("Global temperature change (3,000 year lag ("*Delta*degree*"C))"), las=0)
mtext(side=3, at=par("usr")[1], text="(A) Europe", font=2, adj=0)
box()
close.screen(1)
# 
# screen(2)
# par(mar=c(0,0,0,0), ps=10, mgp=c(3,0.5,0), las=1, tcl=-0.25)
# close.screen(2)

screen(3)
par(mar=c(0,0,0,0), ps=10, mgp=c(3,0.5,0), las=1, tcl=-0.25)

image(x=local.lagNA,
      y=global.lagNA,
      z=matrix(plogis(grid.predNA), nrow=length(local.lagNA)),
      col=rev(colorRampPalette(rev(c("white", "#b3cedc", "#8a87c7", 
                                     "#80388e", "#4d193c", "black")), bias=1.25)(200)),
      zlim=c(0,0.06),
      useRaster=TRUE, axes=FALSE, xlab="", ylab="")

contour(x=local.lagNA,
        y=global.lagNA,
        z=matrix(plogis(grid.predNA), nrow=length(local.lagNA)), add=TRUE,
        levels=seq(0.01,0.08,0.01))

axis(side=1, mgp=c(3,0.1,0))
axis(side=2)
mtext(side=1, line=1.25, text = expression("Local temperature change (5,000 year lag ("*Delta*degree*"C))"))
mtext(side=2, line=2, text = expression("Global temperature change (2,000 year lag ("*Delta*degree*"C))"), las=0)
mtext(side=3, at=par("usr")[1], text="(B) North America", font=2, adj=0)
box()
close.screen(3)

screen(4)
par(mar=c(0,0,0,0), ps=10, mgp=c(3,0.5,0), las=1, tcl=-0.25)
novelrange <- range(c(0,0.06))
plot(x=NULL, y=NULL, xlim=c(0,max(novelrange)), ylim=c(0,0.06), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="")
image(x=seq(0, max(novelrange), len=200),
      y=c(0,1),
      z=matrix(seq(0, max(novelrange), len=200), ncol=1),
      col=rev(colorRampPalette(rev(c("white", "#b3cedc", "#8a87c7", 
                                     "#80388e", "#4d193c", "black")), bias=1.25)(200)),
      add=TRUE, useRaster=TRUE)
axis(side=3, mgp=c(3,0.2,0))
mtext(side=3, line=1.25, text = "Novel community emergence probability")
box()
close.screen(4)

screen(5)
par(mar=c(0,0,0,0), ps=10, mgp=c(3,0.5,0), las=1, tcl=-0.25)
resM <- residuals(targetMReg)
res.df <- data.frame(res = residuals(targetMReg),
                     REGION = as.factor(targetMReg$data$REGION),
                     novel = as.factor(targetMReg$data$novel),
                     bins = targetMReg$data$bins)

resGam <- gam(res ~ -1 + s(bins) + s(novel, bs="re"), data=res.df)

pred.df <- data.frame(novel = rep(c("FALSE", "TRUE"),each=length(unique(res.df$bins))),
                      bins = rep(sort(unique(res.df$bins)), 2))
pred.dfE <- cbind(pred.df, 
                  as.data.frame(predict(resGam, newdata=cbind(pred.df, REGION="Europe"), se.fit=TRUE)))
pred.dfNA <- cbind(pred.df, 
                  as.data.frame(predict(resGam, newdata=cbind(pred.df, REGION="North America"), se.fit=TRUE)))

pred.df <- data.frame(bins = sort(unique(res.df$bins)),
                      novel= "a")
pred.df <- cbind(pred.df, as.data.frame(predict(resGam, newdata=pred.df, se.fit=TRUE, exclude="s(novel)")))

plot(y=pred.df$fit, x=pred.df$bins, type="n",
     xlim=c(17000,0), ylim=c(-0.3,0.2), xlab="", ylab="", xaxt="n")

axis(side=1, at=seq(5000,15000,5000), labels=format(seq(5000,15000,5000), big.mark=","),
     mgp=c(3,0.2,0))
axis(side=1, at=0,mgp=c(3,0.2,0))
axis(side=1, at=seq(0,17000,1000), labels=NA, tcl=-0.125)
mtext(side=1, line=1.5, text=expression("Years before present"))
mtext(side=2, line=2.25, text="Centered model residuals", las=0)

res.df$resCor <- res.df$res - resGam$coefficients[grepl("novel", names(resGam$coefficients))][res.df$novel]

# 95% predictive intervals
PI <- tapply(res.df$resCor,
             list(res.df$bins, res.df$REGION),
             function(x){quantile(x, probs=c(0.025,0.975))}, 
             simplify=FALSE)
PIE <- do.call("rbind", PI[,1])
PINA <- do.call("rbind", PI[,2])

# polygon(y=c(PIE[,1], rev(PIE[,2])),
#         x=c(as.numeric(rownames(PIE)), rev(as.numeric(rownames(PIE)))), border=NA, 
#         col=rgb(0.5,0.5,0.9,0.5))
# polygon(y=c(PINA[,1], rev(PINA[,2])),
#         x=c(as.numeric(rownames(PINA)), rev(as.numeric(rownames(PINA)))), border=NA, 
#         col=rgb(0.5,0.7,0.5,0.5))

quarts <- tapply(res.df$resCor, list(res.df$bins, res.df$REGION), 
                         function(x){quantile(x, probs=c(0.33,0.66))}, simplify=FALSE)
qE <- do.call("rbind", quarts[,1])
qNA <- do.call("rbind", quarts[,2])

polygon(y=c(qE[,1], rev(qE[,2])),
        x=c(as.numeric(rownames(qE)), rev(as.numeric(rownames(qE)))), border=NA, 
        col=rgb(0.2,0.2,0.7,0.5))
polygon(y=c(qNA[,1], rev(qNA[,2])),
        x=c(as.numeric(rownames(qNA)), rev(as.numeric(rownames(qNA)))), border=NA, 
        col=rgb(0.2,0.5,0.2,0.5))

resMean <- tapply(res.df$resCor, list(res.df$bins, res.df$REGION), mean)
lines(y=resMean[,1], x=pred.df$bins, type="l", lwd=2, col="blue")
lines(y=resMean[,2], x=pred.df$bins, type="l", lwd=2, col="darkgreen")
abline(h=0, lty="31")

text(x=relative.axis.point(0.005, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(C)", font=2, adj=0)
text(x=relative.axis.point(0.14, "x"),
     y=relative.axis.point(0.95, "y"),
     labels=expression(bold("North America"*phantom(" & Europe"))), 
     font=2, col="#629D26")
text(x=relative.axis.point(0.14, "x"),
     y=relative.axis.point(0.95, "y"),
     labels=expression(bold(phantom("North America")*" & "*phantom("Europe"))), 
     font=2, col="black")
text(x=relative.axis.point(0.14, "x"),
     y=relative.axis.point(0.95, "y"),
     labels=expression(bold(phantom("North America & ")*"Europe")), 
     font=2, col="blue")

close.screen(5)

close.screen(all.screens=TRUE)
dev.off()

# Novelty by latitude ####

latitude.model <- regions.novel.prob.model.weird(novel.list = all.novel,
                                                 name = "all",
                                                 site.df = site.df,
                                                 all.prob.model = all.nprob.models,
                                                 time.age.limits = c(1200, 25000),
                                                 factor.age.limits = c(0,1000),
                                                 lat.lims = c(25,70),
                                                 time.k = 5,
                                                 fact.k = 5,
                                                 test.model=TRUE,
                                                 sauto.n = 999,
                                                 sauto.size=NA)

saveRDS(latitude.model, 
        date.wrap("./outputs/novel probability by latitude model (all)", ".rds"))

latitude.model.plot(lat.model = latitude.model,
                    all.nprob.models = all.nprob.models,
                    name = "all",
                    major.xlims = c(24000,1200),
                    lat.lims=c(25,70),
                    zlim=c(0,0.08))

gc()

# latitude.model <- readRDS("./outputs/novel probability by latitude model (all) 2021-08-31.rds")

#               Latitude modern bin sub-sampling ####
  
# what is the minimum sampling we should aim for?
table(cut(latitude.model$data$abs.lat[latitude.model$data$bins=="0"], 
          breaks=seq(0,80,10)))

sub.novels  <- subsample.novel.latitude(novel.geo.df = latitude.model$data[latitude.model$data$abs.lat <= 70,],
                                        lat.fact.model = latitude.model$fact.model,
                                        iter=999,
                                        sample = 22,
                                        ylims=c(0,0.25))
saveRDS(sub.novels, "./outputs/latitude_novel_subset.rds")


sub.edf <- do.call("rbind", lapply(sub.novels$model.list, function(x){x$edf}))
table(round(sub.edf[,1],3) == 1)

sort(latitude.model$data$abs.lat[latitude.model$data$bins == "0" &
                                 latitude.model$data$novel])

table(unlist(sapply(1:length(sub.novels$model.list), function(n){
  print(n)
  
  x<-sub.novels$model.list[[n]]

  if(is.null(x$preds)){return(NULL)}
  if(mean(plogis(x$preds$fit)) < 0.001){return(NULL)}
  
  (rev(plogis(x$preds$fit))[1] - rev(plogis(x$preds$fit))[10]) > 0
})))

#               Latitude models by region ####

latitude.model.region <- regions.novel.prob.model.continent(novel.list = all.novel,
                                                            name = "all (continent)",
                                                            site.df = site.df,
                                                            all.prob.model = all.nprob.models,
                                                            lat.lims = c(25,70),
                                                            time.age.limits = c(1000, 25000),
                                                            factor.age.limits = c(0,1200),
                                                            time.k = 5,
                                                            fact.k = 5,
                                                            test.model=FALSE,
                                                            sauto.size=999,
                                                            sauto.n = 5)

latitude.model.plot.region(lat.model = latitude.model.region,
                           all.nprob.models = all.nprob.models,
                           name = "all regions (continent)",
                           major.xlims = c(25000,1000),
                           novelLims=c(0,0.11),
                           ylims=c(25,70))

# FIGURE 1 (didactic local novelty figure) ####

pdf("./plots/localNoveltyExp.pdf", height=4, width=10)

split.screen(rbind(c(0.03,0.2,0.3,0.75),
                   c(0.25,0.8,0.125,0.85),
                   c(0.85,0.975,0.125,0.85)))

screen(1)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=c(0.15,0.85), ylim=c(0.15,0.85), asp=1, axes=FALSE,
     xlab="", ylab="")
library(plotrix)
circ.ang <- seq(0, 2*pi, len=7)[-1]
circ.cents <- cbind(c(0.5, 0.5 + 0.22 * cos(circ.ang)),
                    c(0.5, 0.5 + 0.22 * sin(circ.ang)))
sapply(1:nrow(circ.cents), function(n){
draw.circle(x=circ.cents[n,1],y=circ.cents[n,2],0.1)
text(x=circ.cents[n,1], y=circ.cents[n,2], labels=n)
})
mtext(side=1, line=0, text="Longitude")
mtext(side=2, line=0, text="Latitude", las=0)
box()
axis(side=3, at=circ.cents[3,1] + c(-0.1,0.1),
     line=-1, tcl=0.125, labels=NA)
text(x=circ.cents[3,1], y=circ.cents[3,2] +0.15,
     labels="Time series")

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(A)", font=2, adj=0)
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, ylim=c(0.5,nrow(circ.cents)+0.5), xlim=c(1,10),
     axes=FALSE, xlab="", ylab="", yaxs="i")
axis(side=2, at=nrow(circ.cents):1, labels=1:nrow(circ.cents), las=1)
mtext(side=2, line=1.25, text="Time series", las=0)

states <- read.csv("./raw.datafiles/novelDiagStates.csv",
                   row.names = 1)
cols <- cbind(LETTERS[1:5],
              c("darkred", "darkblue", "darkgreen", "purple", "yellow"),
              c("white","white","white","white", "black"))

novMat <- do.call("rbind", lapply(1:nrow(states), function(n){
  
  xPos <- nrow(circ.cents)-(n-1)
  
  segments(x0=1, x1=10, 
           y0=xPos, y1=xPos)
  
  do.call("rbind", lapply(1:ncol(states), function(n1){
    
    # local novel
    if(n1 > 1){
    firstMatch <- which(states[n,1:n1] == states[n,n1])[1]
    
    # overall novel
    pastMatch <- which(states[,1:n1] == states[n,n1], arr.ind=TRUE)[1,2]
    } else {
    firstMatch=Inf
    pastMatch=Inf
    }
    
    if(pastMatch == n1 & pastMatch > 1){
      rect(xleft=n1-0.35,
           xright=n1+0.35,
           ybottom=xPos-0.5,
           ytop=xPos+0.5,col="darkorange", lwd=2)
    }
    
    if(firstMatch == n1 & firstMatch > 1){
      draw.circle(x=n1,y=xPos, 0.3, col="orange", lwd=2)
    }
    
    draw.circle(x=n1,y=xPos, 0.2, col=cols[cols[,1]==states[n,n1],2])
    text(x=n1, y=xPos, labels=states[n,n1], col=cols[cols[,1]==states[n,n1],3])
    
    return(cbind(localNov = ifelse(firstMatch == n1 & firstMatch > 1, TRUE, FALSE),
                 trueNov = ifelse(pastMatch == n1 & pastMatch > 1, TRUE, FALSE)))
  }))
  
}))
box()
mtext(side=3, at=relative.axis.point(0.01, "x"), text="(B)", font=2)

par(lheight=0.9)
par(xpd=NA)
mtext(side=1, at=relative.axis.point(0, "x"),
      text="Past", adj=0)
Arrows(x0=relative.axis.point(0.06, "x"),
       x1=relative.axis.point(0.9, "x"),
       y0=relative.axis.point(-0.04, "y"),
       y1=relative.axis.point(-0.04, "y"),
       arr.type="triangle", arr.width=0.1, arr.length=0.1)

mtext(side=1, at=relative.axis.point(1, "x"),
      text="Present", adj=1)

draw.circle(x=relative.axis.point(0.1,'x'),
            y=relative.axis.point(1.1,"y"), 
            0.2, col=cols[1,2])
text(x=relative.axis.point(0.1,'x'),
     y=relative.axis.point(1.1,"y"), 
     labels="A", col="white")
text(x=relative.axis.point(0.115,'x'),
     y=relative.axis.point(1.09,"y"), 
     labels="Community\nstate", pos=4)

draw.circle(x=relative.axis.point(0.3,'x'),
            y=relative.axis.point(1.1,"y"),
            0.3, col="orange", lwd=2)
draw.circle(x=relative.axis.point(0.3,'x'),
            y=relative.axis.point(1.1,"y"),
            0.2, col="white")
text(x=relative.axis.point(0.325,'x'),
     y=relative.axis.point(1.09,"y"), 
     labels="Novel\nlocally", pos=4)

rect(xleft=relative.axis.point(0.5,'x')-0.35,
     xright=relative.axis.point(0.5,'x')+0.35,
     ybottom=relative.axis.point(1.1,"y")-0.5,
     ytop=relative.axis.point(1.1,"y")+0.5,
     col="darkorange", lwd=2)
draw.circle(x=relative.axis.point(0.5,'x'),
            y=relative.axis.point(1.1,"y"),
            0.3, col="white")
text(x=relative.axis.point(0.525,'x'),
     y=relative.axis.point(1.08,"y"), 
     labels="Novel\nacross\nspace", pos=4)

par(xpd=FALSE)
par(lheight=1)

close.screen(2)

screen(3)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=c(0.5,2.5), ylim=c(0,max(colSums(novMat))+1),
     axes=FALSE, xlab="", ylab="", yaxs="i")

# true novels
barPos <- 1:2
barWidth = 0.35
rect(xleft=barPos - barWidth, xright = barPos + barWidth,
     ybottom=0, ytop = colSums(novMat)[2],
     col="darkorange")

rect(xleft=barPos[2] - barWidth, xright = barPos[2] + barWidth,
     ybottom=0, ytop = colSums(novMat)[1],
     col="orange")
axis(side=2, las=1)
mtext(side=2, line=1.25, text="Count", las=0)

par(lheight=0.9)
axis(side=1, at=barPos[1], labels=c("Novel\nacross\nspace"), mgp=c(3,1.29,0))
axis(side=1, at=barPos[2], labels=c("Novel\nlocally"), mgp=c(3,0.7,0))
par(lheight=1)
box()
text(x=relative.axis.point(0.035, "x"),
     y=relative.axis.point(0.965, "y"),
     labels="(C)", font=2, adj=0)

close.screen(3)
dev.off()

# FIGURE 2 (pollen zone comparison) ####

# we need to find a decent length time series with a mixture of novelty cats
# and variation in expectations

siteStats <- do.call("rbind", lapply(all.novel$novel, function(x){
  
  data.frame(site = x$site[1],
             len = nrow(x),
             novC = sum(x$novel, na.rm=TRUE),
             instC = sum(x$instant, na.rm=TRUE),
             cumuC = sum(x$cumul, na.rm=TRUE),
             instVar = var(x$seq.exp, na.rm=TRUE),
             cumuVar = var(x$min.exp, na.rm=TRUE))
  
}))

siteStats[order(siteStats$len, decreasing=TRUE),]
summary(siteStats$instVar)

siteStats[siteStats$len > 55 &
            siteStats$instVar > 0.005 & 
            siteStats$cumuVar > 1e-3,]

subsite <- "23911"
subGenus <- droplevels(plant.with.genus[plant.with.genus$site == subsite,])
subGenus[1,]
subSsmat <- all.novel$raw.ssmats[[subsite]]

subSsmat <- subSsmat[,!colnames(subSsmat) %in% c("Myriophyllum", "Potamogeton")]

subdata <- identify.novel.gam( prop.table(subSsmat, 1),
                               alpha=0.05, metric="bray",
                               site="1", plot=FALSE, plot.data=TRUE)

propMat <- prop.table(subSsmat, margin=1)

# okay, let's get each genus overall proportion, order and drop genera with
# less than 10% representation in at least one sample bin
keepgen <- colSums(propMat > 0.1) > 0

propMat <- cbind(propMat[,keepgen], rowSums(propMat[,!keepgen]))

cumul.col <- rgb(0.373,0.651,0.765)
cumul.col.rgb <- col2rgb(cumul.col)/255

areaplot <- function(mat, cols){
  
  mat <- mat[,order(colMeans(mat), decreasing=TRUE)]
  cumMat <- rbind(0, apply(mat, 1, cumsum))
  xPoints <- as.numeric(rownames(mat))
  
  sapply(1:ncol(mat), function(n){
    
    polygon(x = c(cumMat[n,], rev(cumMat[n+1,])), 
            y = c(xPoints, rev(xPoints)), col=cols[n])
    
  })
  
}

ylims=c(0,max(as.numeric(rownames(propMat))))
pdf("./plots/example pollen.pdf", height=5, width=7, useDingbats = FALSE)

split.screen(rbind(c(0.10,0.32,0.1,0.96),
                   c(0.3935,0.55,0.1,0.96),
                   c(0.55,0.7065,0.1,0.96),
                   c(0.10,0.99,0.1,0.96)))

screen(1)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25)
plot(x=NULL, y=NULL, xlim=c(0,1), ylim=rev(ylims), xaxs="i", yaxs="i", 
     axes=FALSE, xlab="", ylab="")
axis(side=4, at=seq(0,20000,2000), labels=NA)
axis(side=4, at=seq(0,20000,500), labels=NA, tcl=-0.125)
mtext(side=4, at=seq(0,12000,2000), 
      adj=0.5, padj=0.5, line=1.25, text=format(seq(0,12000,2000), big.mark=","), las=1)
mtext(side=4, at=relative.axis.point(1.1, "y"), line=0, text="Years before present", las=1)


axis(side=1, at=seq(0,1,0.25), mgp=c(3,0,0))
mtext(side=1, line=1, text="Relative abundance")
mtext(side=2, line=2, text="Pollen zones")

library(viridis)
areaCols <- viridis(ncol(propMat)+1, option="V")
areaCols <- areaCols[rep(1:round(length(areaCols)/2), each=2) + c(0, round(length(areaCols)/2))]

areaplot(mat=propMat, cols=areaCols)

# zones from Fossitt 1994
foss <- read.csv("./raw.datafiles/Fossitt_zones.csv")
segments(x0=0,x1=1, y0=foss$end, col="white", lty="31")
mtext(side=2, line=0.1, at=rowMeans(foss[,2:3]), las=1, text=foss$zone)

mtext(side=3, line=0.05, at=par("usr")[1], text="(A)", font=2, adj=0)

box()
close.screen(1)

screen(4)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25)

novel.col <- col2rgb("orange")/255

plot(x=NULL, y=NULL, xlim=c(0,1), ylim=rev(ylims), xaxs="i", yaxs="i", axes=FALSE)

text(x=0.875, y=relative.axis.point(0.98, "y"), 
     labels=expression("Novelty framework result"), adj=0.5)

rect(xleft=0, xright=1, ybottom=12000, ytop=12400, border=NA, 
     col=rgb(novel.col[1],novel.col[2],novel.col[3],0.5))
text(x=0.69, y=12200, labels=expression("(1) Emergence of novel state"), adj=0)

rect(xleft=0, xright=1, ybottom=10000, ytop=10400, border=NA,
     col=rgb(0.5,0.5,0.5,0.25))
text(x=0.69, y=10800, labels=expression("(2) Rapid return to previously observed\nstate: no emergence of novelty as not\nsufficiently different to past states"), adj=c(0,1))

rect(xleft=0, xright=1, ybottom=8800, ytop=9400, border=NA,
     col=rgb(novel.col[1],novel.col[2],novel.col[3],0.5))
text(x=0.69, y=9400, labels=expression("(3) Two successive novel states during\nrapid ecological turnover to new state"), adj=0)

rect(xleft=0, xright=1, ybottom=8400, ytop=8000, border=NA,
     col=rgb(0.5,0.5,0.5,0.25))
text(x=0.69, y=8200, labels=expression("(4) Continuation of the novel state in (3)\nbut no emergence of new novel state"), adj=c(0,1))


rect(xleft=0, xright=1, ybottom=c(1800, 4800, 7200), ytop=c(2200, 4400, 6800), border=NA,
     col=rgb(0.5,0.5,0.5,0.25))
text(x=0.69, y=3000, labels=expression("(5) Gradual vegetational turnover does\nnot result in the emergence of ecological\nnovelty, even if ending composition is\nvastly different from the start of the time\nseries"), adj=0)

close.screen(4)

screen(2)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25)
plot(x=NULL, y=NULL, xlim=c(0,0.8), ylim=rev(ylims),  yaxs="i", axes=FALSE)
axis(side=2, at=seq(0,20000,2000), labels=NA)
axis(side=2, at=seq(0,20000,500), labels=NA, tcl=-0.125)

axis(side=1, at=seq(0,1,0.25), mgp=c(3,0,0))
mtext(side=1, line=1.25, text="Instantaneous\ndissimilarity")

novbins <- as.numeric(subdata[[1]]$bins)
subnovel <- subdata[[1]]

polygon(x=c(subdata[[2]]$lwr, rev(subdata[[2]]$upr)),
        y=c(novbins, rev(novbins)), border="red", col=rgb(1,0,0,0.3))
lines(as.numeric(subnovel$bins) ~ subnovel$seq.exp, col="red", lwd=2)
lines(as.numeric(subnovel$bins) ~ subnovel$seq.dist)
with(subnovel[subnovel$instant,], points(as.numeric(bins) ~ seq.dist, pch=21, bg="red"))

mtext(side=3, line=0.05, at=par("usr")[1], text="(B)", font=2, adj=0)

box()
close.screen(2)

screen(3)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25)
plot(x=NULL, y=NULL, xlim=c(0,0.8), ylim=rev(ylims), yaxs="i", axes=FALSE)
axis(side=1, at=seq(0,1,0.25), mgp=c(3,0,0))
mtext(side=1, line=1.25, text="Cumulative\ndissimilarity")

polygon(x=c(subdata[[3]]$lwr, rev(subdata[[3]]$upr)),
        y=c(novbins, rev(novbins)), border=cumul.col, col=rgb(cumul.col.rgb[1],
                                                              cumul.col.rgb[2],
                                                              cumul.col.rgb[3],0.3))
lines(as.numeric(subnovel$bins) ~ subnovel$min.exp, col=cumul.col, lwd=2)
lines(as.numeric(subnovel$bins) ~ subnovel$raw.min.dist)
with(subnovel[subnovel$cumul,], points(as.numeric(bins) ~ raw.min.dist, pch=21, bg=cumul.col))

mtext(side=3, line=0.05, at=par("usr")[1], text="(C)", font=2, adj=0)


box()
close.screen(3)

dev.off()
# ####
# SUPPLEMENTARY ANALYSES ####
# Novelty over space ####

# number of rows needed in mat
novel.rows <- lapply(all.novel$prop.ssmats, rownames)
novel.length <- data.frame(site = unlist(sapply(1:length(novel.rows), 
                             function(n){rep(names(novel.rows)[n], length(novel.rows[[n]]))})),
                           bin = unlist(novel.rows))
novel.length <- novel.length[as.numeric(novel.length$bin) <= 25000,]
head(novel.length)

# number of cols
novel.tax <- sort(unique(unlist(sapply(all.novel$prop.ssmats, colnames))))

# empty matrix
full.mat <- matrix(0, nrow=nrow(novel.length), ncol=length(novel.tax))
colnames(full.mat) = novel.tax

# fill with values
for(i in 1:length(all.novel$prop.ssmats)){
  print(i)
  temp.data <- all.novel$prop.ssmats[[i]]
  temp.data <- temp.data[as.numeric(rownames(temp.data)) <= 25000, ]
  
  full.mat[which(paste0(novel.length$site, ".", novel.length$bin) %in%
                 paste0(all.novel$sites[i], ".", rownames(temp.data))), 
           match(colnames(temp.data), colnames(full.mat))] = temp.data

}
rownames(full.mat) = paste0(novel.length$site, ".",  novel.length$bin)

# now we can start to build up our I and C scores

# reorder into time order blocks
full.mat <- full.mat[order(novel.length$bin, as.numeric(novel.length$site)),]
novel.length <- novel.length[order(novel.length$bin, as.numeric(novel.length$site)),]

# now pull out each block one at a time to do novelty calcs
library(parallel)
cross.space.novelty <- do.call("rbind", lapply(seq(1200,24800,200), function(bin){
  
  print(bin)
  
  target.mat <- full.mat[novel.length$bin == bin,]
  ref.mat <- full.mat[novel.length$bin == (bin + 200),]
  
  comb.mat <- rbind(target.mat, ref.mat)
  
  # comb.inst1 <- apply(target.mat, 1, function(x){
  #   min(apply(ref.mat, 1, function(y){vegdist(rbind(x,y))}))
  #   })
  comb.dist <- as.matrix(vegdist(comb.mat, metric="bray"))
  
  #instantaneous calculations
  comb.inst <- apply(comb.dist[1:nrow(target.mat), 
                               (nrow(target.mat)+1):ncol(comb.dist)], 1, min)
  
  # cumulative calculations
  ref.mat <- full.mat[as.numeric(novel.length$bin) > bin,]
  
  # doing the cumul distance matrix turns out to be incredibly cpu intensive for
  # a lot of values we don't need (between the ref mat), so we'll do all comparisons
  # one by one in a loop
  comb.cumul <- unlist(mclapply(1:nrow(target.mat), function(n){
    
  x <- target.mat[n,]
  a <- min(apply(ref.mat, 1, function(y){ vegdist(rbind(x, y))}))

  return(a)
      }, mc.cores=4))
  
  plot(comb.inst ~ comb.cumul)
  
  saveRDS(data.frame(CSinstant = comb.inst,
                     CScumul = comb.cumul),
          paste0("./outputs/crossspaceNovelty:", bin, ".rds"))
  
  return(data.frame(CSinstant = comb.inst,
                    CScumul = comb.cumul))
}))

crossSpaceFiles <- list.files("./outputs", pattern="crossspaceNovelty")
crossSpaceorder <- as.numeric(gsub("crossspaceNovelty\\:|\\.rds", "", crossSpaceFiles))
crossSpaceFiles <- crossSpaceFiles[order(crossSpaceorder)]

cross.space.novelty <- do.call("rbind", lapply(crossSpaceFiles, function(x){
  readRDS(paste0("./outputs/",x))
  }))

cross.space.novelty$id <- rownames(cross.space.novelty)
cross.space.novelty$site <- substr(cross.space.novelty$id, 1, regexpr("\\.", cross.space.novelty$id)-1)
cross.space.novelty$bin <- as.numeric(substr(cross.space.novelty$id, 
                                             regexpr("\\.", cross.space.novelty$id)+1,
                                             nchar(cross.space.novelty$id)))

novel.df <- do.call("rbind", all.novel$novel)

novel.df$id <- paste0(novel.df$site, ".", novel.df$bins)

novel.df <- merge(novel.df, cross.space.novelty[,1:3],
                  by.x="id", by.y="id", all.x=TRUE, all.y=TRUE, sort=FALSE)
novel.df <- merge(novel.df, site.df[,c("REGION", "site", "long", "lat")],
                  by.x="site", by.y="site", all.x=TRUE, all.y=FALSE, sort=FALSE)

novel.df$bins <- as.numeric(novel.df$bins)

pdf("./plots/spacetimeNovcomp.pdf", height=3.5, width=6, useDingbats = FALSE)
par(mfrow=c(1,2), mar=c(1,0,1,0), oma=c(2,3,0,1), ps=10, tcl=-0.25, mgp=c(3,0,0), las=1)
plot(novel.df$seq.dist ~ novel.df$CSinstant, pch=16, col=rgb(1,0,0,0.02), axes=FALSE, xlab="", ylab="",
     xlim=c(0,1), ylim=c(0,1))
axis(side=2, mgp=c(3,0.5,0))
axis(side=1)
box()

icor <- cor(cbind(novel.df$seq.dist,novel.df$CSinstant), use="complete.obs")
isame <- sum(round(novel.df$seq.dist, 3) == round(novel.df$CSinstant, 3), na.rm=TRUE) / nrow(novel.df)

text(x=relative.axis.point(0.98,"x"),
     y=relative.axis.point(0.05,"y"),
     labels=paste0("r = ", round(icor[2],3), "\n",
                   round(isame*100,2),"% identical"), pos=2, offset=0.25)

text(x=relative.axis.point(0,"x"),
     y=relative.axis.point(0.95,"y"),
     labels="(A) Instantaneous dissimilarity", pos=4, offset=0.25, font=2)
mtext(side=1, line=1, text="Dissimilarity across space", at=par("usr")[2])
mtext(side=2, line=2, text="Dissimilarity constrained to time series", las=0)

plot(novel.df$raw.min.dist ~ novel.df$CScumul, pch=16, col=rgb(0,0,1,0.02), axes=FALSE, xlab="", ylab="",
     xlim=c(0,1), ylim=c(0,1))
axis(side=2, mgp=c(3,0.5,0), labels=NA)
axis(side=1)

box()

ccor <- cor(cbind(novel.df$raw.min.dist,novel.df$CScumul), use="complete.obs")
csame <- sum(round(novel.df$raw.min.dist, 3) == round(novel.df$CScumul, 3), na.rm=TRUE) / nrow(novel.df)

text(x=relative.axis.point(0,"x"),
     y=relative.axis.point(0.95,"y"),
     labels="(B) Cumulative dissimilarity", pos=4, offset=0.25, font=2)

text(x=relative.axis.point(1,"x"),
     y=relative.axis.point(0.05,"y"),
     labels=paste0("r = ", round(ccor[2],3), "\n",
                   round(csame*100,2),"% identical"), pos=2, offset=0.25)
dev.off()

#           model through time cross space ####

novel.df$binNum <- as.numeric(novel.df$bins)

# mean models
novel.df <- novel.df[novel.df$binNum <= 25000 &
                    rowSums(is.na(novel.df)) < 15,]
novel.df$binFact <- as.factor(novel.df$binNum)

csNovI <- gam(beta.tr(CSinstant) ~ s(binNum), family=betar, data=novel.df)
csNovC <- gam(beta.tr(CScumul) ~ s(binNum), family=betar, data=novel.df)
csNovIFact <- gam(beta.tr(CSinstant) ~ binFact, family=betar, data=droplevels(novel.df[novel.df$binNum <= 1000,]))
csNovCFact <- gam(beta.tr(CScumul) ~ binFact, family=betar, data=droplevels(novel.df[novel.df$binNum <= 1000,]))

NovI <- gam(beta.tr(seq.dist) ~ s(binNum), family=betar, data=novel.df)
NovC <- gam(beta.tr(raw.min.dist) ~ s(binNum), family=betar, data=novel.df)
NovIFact <- gam(beta.tr(seq.dist) ~ binFact, family=betar, data=droplevels(novel.df[novel.df$binNum <= 1000,]))
NovCFact <- gam(beta.tr(raw.min.dist) ~ binFact, family=betar, data=droplevels(novel.df[novel.df$binNum <= 1000,]))

# mean predictions
pred.df <-  data.frame(binNum = seq(1200,25000,200), site=NA)
IpredC <- cbind(pred.df, predict(csNovI, newdata=pred.df, se.fit=TRUE))
CpredC <- cbind(pred.df, predict(csNovC, newdata=pred.df, se.fit=TRUE))
Ipred <- cbind(pred.df, predict(NovI, newdata=pred.df, se.fit=TRUE))
Cpred <- cbind(pred.df, predict(NovC, newdata=pred.df, se.fit=TRUE))

fact.pred.df <- data.frame(binFact=factor(seq(0,1000,200)))
IpredCFact <- cbind(fact.pred.df, predict(csNovIFact, newdata=fact.pred.df, se.fit=TRUE))
CpredCFact <- cbind(fact.pred.df, predict(csNovCFact, newdata=fact.pred.df, se.fit=TRUE))
IpredFact <- cbind(fact.pred.df, predict(NovIFact, newdata=fact.pred.df, se.fit=TRUE))
CpredFact <- cbind(fact.pred.df, predict(NovCFact, newdata=fact.pred.df, se.fit=TRUE))

# cross-space novel community classification

# find outliers using 95% prediction intervals

seq.mu <- csNovI$fitted.values
phi <- as.numeric(substr(csNovI$family$family,
                         regexpr("\\(", csNovI$family$family)+1,
                         nchar(csNovI$family$family)-1))
A = seq.mu * phi
B = phi - A
pCI <- do.call("rbind", lapply(1:length(A), function(n){
  data.frame(lwr=qbeta(0.05, shape1=A[n], shape2=B[n]),
             upr=qbeta(0.95, shape1=A[n], shape2=B[n]),
             seq.p = pbeta(novel.df$CSinstant[n], shape1=A[n], shape2=B[n], lower.tail=FALSE))
}))

seq.mu <- csNovC$fitted.values
phi <- as.numeric(substr(csNovC$family$family,
                         regexpr("\\(", csNovC$family$family)+1,
                         nchar(csNovC$family$family)-1))
A = seq.mu * phi
B = phi - A
pCC <- do.call("rbind", lapply(1:length(A), function(n){
  data.frame(lwr=qbeta(0.05, shape1=A[n], shape2=B[n]),
             upr=qbeta(0.95, shape1=A[n], shape2=B[n]),
             seq.p = pbeta(novel.df$CScumul[n], shape1=A[n], shape2=B[n], lower.tail=FALSE))
}))

novel.df$CIN <- pCI$seq.p < 0.05
novel.df$CCN <- pCC$seq.p < 0.05
novel.df$CNovel <- novel.df$CIN & novel.df$CCN 

cross.novel.list <- lapply(all.novel$sites, function(site){
  
  # get novel object
   temp <- all.novel$novel[[which(all.novel$sites==site)]]
  
  # get cross-space calculations
   tempCS <- novel.df[novel.df$site == site,]
   tempCS <- tempCS[order(tempCS$bins, decreasing=TRUE),]
   
   temp$instant = tempCS$CIN
   temp$cumul = tempCS$CCN
   temp$novel = tempCS$CNovel
   
   return(temp[,c("site", "bins", "instant", "cumul", "novel")])
   
})

cross.novel <- all.novel
cross.novel$novel <- cross.novel.list

cross.nprob.models <- novel.prob.models.space(novel.list = cross.novel,
                                      site.df = site.df,
                                      time.k = 40,
                                      test.model=FALSE,
                                      name = "allcrossspace",
                                      time.age.limits = c(1200,25000),
                                      factor.age.limits = c(0,1000),
                                      sauto.n = 5000,
                                      sauto.iter = 999)

# version of Fig 1 with different sized bins ####

plant.with.genus <- read.csv(paste0(large_file_path, "/processed_genus_records.csv"))
#plant.with.genus$site = plant.with.genus$site.id

plant.slimmed <- plant.with.genus[!is.na(plant.with.genus$age),
                                  colnames(plant.with.genus) %in%
                                    c("family", "site", "age", "count")]

five.novel <- neotoma.novelty(dataset = plant.slimmed,
                              ssmat.type = "abund",
                              bins = seq(-100,max(plant.slimmed$age), 500),
                              rich.cutoff = c(100, 5000),
                              age.limits = c(-150, 25100),
                              taxon.res = "family",
                              bin.cutoff = 10,
                              taxa.cutoff = 5,
                              novel.alpha = 0.05,
                              novel.metric = "bray",
                              sqrt.mat=TRUE)

ten.novel <- neotoma.novelty(dataset = plant.slimmed,
                             ssmat.type = "abund",
                             bins = seq(-100,max(plant.slimmed$age), 1000),
                             rich.cutoff = c(100, 5000),
                             age.limits = c(-150, 25100),
                             taxon.res = "family",
                             bin.cutoff = 10,
                             taxa.cutoff = 5,
                             novel.alpha = 0.05,
                             novel.metric = "bray",
                             sqrt.mat = TRUE)

# prob models

five.nprob.models <- novel.prob.models(novel.list = five.novel,
                                       site.df = site.df,
                                       time.k = 20,
                                       test.model=FALSE,
                                       name = "500bins",
                                       time.age.limits = c(1000,25000),
                                       factor.age.limits = c(0,1000),
                                       sauto.n = 5000,
                                       sauto.iter=999)

novel.prob.plot(prob.model.list = five.nprob.models, 
                env.data = temp.data,
                mod.env.data = mod.env.data,
                ylims=c(0,0.11), 
                regylims=c(0,0.175),
                name = "500bins",
                time.age.limits = c(0,25000),
                factor.age.limits = c(0,1000),
                group.letters=c("A","AB","B","B","B","B"))

saveRDS(five.nprob.models,
        date.wrap("./outputs/novel probability models (500 yr bins)", ".rds"))

ten.nprob.models <- novel.prob.models(novel.list = ten.novel,
                                      site.df = site.df,
                                      time.k = 10,
                                      test.model=FALSE,
                                      name = "1000bins",
                                      time.age.limits = c(1000,25000),
                                      factor.age.limits = c(0,1000),
                                      sauto.n = 5000,
                                      sauto.iter=999)

novel.prob.plot(prob.model.list = ten.nprob.models, 
                env.data = temp.data,
                mod.env.data = mod.env.data,
                ylims=c(0,0.11), 
                regylims=c(0,0.175),
                name = "1000bins",
                time.age.limits = c(0,25000),
                factor.age.limits = c(0,1000),
                group.letters=c("A","AB","B","B","B","B"))

saveRDS(ten.nprob.models,
        date.wrap("./outputs/novel probability models (1000 yr bins)", ".rds"))



# analyses with different alpha ####

test.alpha <- c(0.01,0.025,0.075,0.1)

alpha.list <- lapply(test.alpha, function(alpha){
  print(alpha)
  # re-calculate novelty using p-values from framework
  
  alpha.novel <- all.novel
  alpha.novel$novel <- lapply(alpha.novel$novel, function(x){
    
    x$novel = x$seq.p <= alpha &
      x$min.p <= alpha
    
    return(x)
  })
  
  temp.model <- novel.prob.models(novel.list = alpha.novel,
                                  site.df = site.df,
                                  time.k = 40,
                                  name = paste0("all - ", alpha),
                                  test.model=FALSE,
                                  time.age.limits = c(1000,25000),
                                  factor.age.limits = c(0,1000))
  
  return(temp.model)
})

saveRDS(alpha.list, "./outputs/alpha difference nprob models.rds")

sapply(1:length(alpha.list), function(n){
  print(n)
  
  novel.prob.plot(alpha.list[[n]], 
                  env.data = temp.data,
                  mod.env.data=mod.env.data,
                  ylims=c(0, 
                          c(0.03,0.045,0.09,0.11)[n]),
                  regylims=c(0,
                             c(0.04, 0.06, 0.12, 0.12)[n]),
                  name = paste0("all alpha = ", c(0.01,0.025,0.075,0.1)[n]),
                  time.age.limits = c(0,25000),
                  factor.age.limits = c(0,1000),
                  group.letters=c(""))
})




# Test the reliability of age modelling (using NAPD subset) ####

NAPDSites <- unique(plant.with.genus$site[plant.with.genus$chronology.name=="NAPD 1"])
NAsites <- site.df$site[site.df$REGION == "North America"]
NAPDSites <- NAPDSites[NAPDSites %in% NAsites]

NAPDNovList <- all.novel
NAPDNovList$novel = NAPDNovList$novel[NAPDNovList$sites %in% NAPDSites]
NAPDNovList$prop.ssmats = NAPDNovList$prop.ssmats[NAPDNovList$sites %in% NAPDSites]
NAPDNovList$raw.ssmats = NAPDNovList$raw.ssmats[NAPDNovList$sites %in% NAPDSites]
NAPDNovList$sites = NAPDNovList$sites[NAPDNovList$sites %in% NAPDSites]

NAPDProb <-  novel.prob.models(novel.list = NAPDNovList,
                               site.df = droplevels(site.df[site.df$site %in% NAPDSites,]),
                               time.k = 40,
                               test.model=FALSE,
                               name = "NAPD",
                               time.age.limits = c(1200,25000),
                               factor.age.limits = c(0,1000),
                               sauto.n = 5000,
                               sauto.iter = 999)

novel.prob.plot.comp(prob.model.list = all.nprob.models, 
                     prob.model.list2 = NAPDProb,
                     env.data = temp.data,
                     mod.env.data = mod.env.data,
                     name = "NAPD",
                     time.age.limits = c(0,25000),
                     factor.age.limits = c(0,1000),
                     time.ylims = c(0,0.08),
                     factor.ylims=c(0, 0.08))


# Genus-level analyses ####

# process raw genus records (as per genus records)
plant.record.df <- read.csv(paste0(large_file_path, "/processed_neotoma_records.csv"))
plant.record.df <- droplevels(plant.record.df[plant.record.df$dataset.type == "pollen",])
plant.with.genus <- plant.record.df[!is.na(plant.record.df$genus) &
                                       plant.record.df$phylum %in% c("A", "G"),]
plant.with.genus <- droplevels(plant.with.genus[!is.na(plant.with.genus$age), ])
plant.with.genus$site = plant.with.genus$site.id
plant.with.genus <- plant.with.genus[order(plant.with.genus$site),]

# some duplicate sample IDs
dupeSamps <- paste(plant.with.genus$sample.id, plant.with.genus$taxon, 
                   plant.with.genus$count, plant.with.genus$ContactName, sep=".")

plant.with.genus <- droplevels(plant.with.genus[!duplicated(dupeSamps),])

write.csv(plant.with.genus, paste0(large_file_path, "/processed_genus_records.csv"))

# novelty analysis at genus level
plant.with.genus <- read.csv(paste0(large_file_path, "/processed_genus_records.csv"))
plant.with.genus <- droplevels(plant.with.genus)
plant.genus.slimmed <- plant.with.genus[,colnames(plant.with.genus) %in%
                                            c("genus", "site", "age", "count")]

genus.novel <- neotoma.novelty(dataset = plant.genus.slimmed,
                                ssmat.type = "abund",
                                bins = seq(-100,max(plant.genus.slimmed$age), 200),
                                rich.cutoff = c(100, 5000),
                                age.limits = c(-150, Inf),
                                taxon.res = "genus",
                                bin.cutoff = 10,
                                taxa.cutoff = 2,
                                novel.alpha = 0.05,
                                novel.metric = "bray",
                                novel.plot=TRUE,
                               sqrt.mat = TRUE)

saveRDS(genus.novel, "./outputs/genus neotoma novelty.rds")
#genus.novel <- readRDS("./outputs/genus neotoma novelty.rds")

#         genus over time model ####

genus.nprob.models <- novel.prob.models(novel.list = genus.novel,
                                         site.df = site.df,
                                         time.k = 20,
                                         test.model=FALSE,
                                         name = "genus",
                                         time.age.limits = c(0,25000),
                                         factor.age.limits = c(0,1000),
                                         sauto.n = 5000)

novel.prob.plot(prob.model.list = genus.nprob.models, 
                env.data = temp.data,
                mod.env.data = mod.env.data,
                ylims=c(0,0.049), 
                regylims=c(0,0.07),
                name = "genus",
                time.age.limits = c(0,25000),
                factor.age.limits = c(0,1000),
                group.letters=c("A","AB","B","B","B","B"))

saveRDS(genus.nprob.models,
        date.wrap("./outputs/novel probability models (genus)", ".rds"))
#genus.nprob.models <- readRDS()

#         genus latitudinal model ####

latitude.model.genus <- regions.novel.prob.model.weird(novel.list = genus.novel,
                                                        name = "genus",
                                                        site.df = site.df,
                                                        all.prob.model = genus.novel,
                                                        time.age.limits = c(0, 25000),
                                                        factor.age.limits = c(0,1200),
                                                        time.k = 5,
                                                        fact.k = 5,
                                                        test.model=FALSE,
                                                        sauto.size=999,
                                                        sauto.n = 5)

saveRDS(latitude.model.genus, 
        date.wrap("./outputs/novel probability by latitude model (genus)", ".rds"))
#latitude.model.genus <- readRDS("./outputs/novel probability by latitude model (genus) 2020-12-15.rds")

latitude.model.plot(lat.model = latitude.model.genus,
                    all.nprob.models = genus.nprob.models,
                    name = "genus",
                    major.xlims = c(25000,1000),
                    lat.lims = c(25,70),
                    zlim=c(0,0.08))

# plot novel expectations ####

comb.df <- do.call("rbind", lapply(all.novel$novel, function(x){
  
  x$seq.exp.std <- x$seq.exp / max(x$seq.exp, na.rm=TRUE)
  x$min.exp.std <- x$min.exp / max(x$min.exp, na.rm=TRUE)
  
  return(x)
  
}))
comb.df$bin.num <- as.numeric(as.character(comb.df$bins))

plot(comb.df$seq.exp.std ~ comb.df$bin.num, xlim=c(25000,0), pch=16, cex=0.2)

inst.exp.model <- gam(beta.tr(seq.exp.std) ~ s(bin.num, k=40),
                      data=comb.df[comb.df$bin.num <= 25000,],
                      family=betar)
plot(inst.exp.model)

cumul.exp.model <- gam(beta.tr(min.exp.std) ~ s(bin.num, k=40), 
                       data=comb.df[comb.df$bin.num <= 25000,],
                       family=betar)
plot(cumul.exp.model)

pred.df <- data.frame(bin.num = seq(0,25000, 200))
inst.pred <- cbind(pred.df,
                   as.data.frame(predict(inst.exp.model, newdata=pred.df, se.fit=TRUE)))
inst.pred$upper <- plogis(inst.pred$fit + 1.96 * inst.pred$se.fit)
inst.pred$lower <- plogis(inst.pred$fit - 1.96 * inst.pred$se.fit)
inst.pred$fit <- plogis(inst.pred$fit)

cumul.pred <- cbind(pred.df,
                    as.data.frame(predict(cumul.exp.model, newdata=pred.df, se.fit=TRUE)))
cumul.pred$upper <- plogis(cumul.pred$fit + 1.96 * cumul.pred$se.fit)
cumul.pred$lower <- plogis(cumul.pred$fit - 1.96 * cumul.pred$se.fit)
cumul.pred$fit <- plogis(cumul.pred$fit)

pdf("./plots/expectations through time.pdf", height=4, width=6, useDingbats = FALSE)

par(mfrow=c(2,1), mar=c(0,0,0,0), oma=c(3,3.5,1,1), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(x=NULL, y=NULL, xlim=rev(c(-100, 25000)), ylim=c(0.66,0.95), 
     axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
axis(side=1, mgp=c(3,0,0), at=seq(5000, 35000, 5000), 
     labels=NA)
axis(side=1, mgp=c(3,0,0), at=0, labels=NA)
axis(side=1, at=seq(0,40000,1000), tcl=-0.125, labels=NA)

axis(side=2)

mtext(side=2, line=2, las=0,
      text="Mean expected instantaneous\ndissimilarity (standardized)")
polygon(y = c(inst.pred$upper, rev(inst.pred$lower)),
        x = c(inst.pred$bin.num, rev(inst.pred$bin.num)),
        border = NA, col=rgb(1,0,0,0.2))
lines(y = inst.pred$fit, x = inst.pred$bin.num, col="red", lwd=2)

text(x = relative.axis.point(0.005, "x"),
     y = relative.axis.point(0.925, "y"),
     labels="(A)", font=2, adj=0)
box()

plot(x=NULL, y=NULL, xlim=rev(c(-100, 25000)), ylim=c(0.6,1), 
     axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
axis(side=1, mgp=c(3,0,0), at=seq(5000, 35000, 5000), 
     labels=NA)
axis(side=1, mgp=c(3,0,0), at=0, labels=NA)
axis(side=1, at=seq(0,40000,1000), tcl=-0.125, labels=NA)

axis(side=2)
axis(side=2, at=seq(0,0.1,0.005), labels=NA, tcl=-0.125)

mtext(side=2, line=2, las=0,
      text="Mean expected cumulative\ndissimilarity (standardized)")

polygon(y = c(cumul.pred$upper, rev(cumul.pred$lower)),
        x = c(cumul.pred$bin.num, rev(cumul.pred$bin.num)),
        border = NA, col=rgb(0,0,1,0.2))
lines(y = cumul.pred$fit, x = cumul.pred$bin.num, col="blue", lwd=2)

axis(side=1, mgp=c(3,0,0), at=seq(5000, 35000, 5000), 
     labels=format(seq(5000, 35000, 5000), big.mark=","))
axis(side=1, mgp=c(3,0,0), at=0)
mtext(side=1, line=1, text="Years before present")
axis(side=1, at=seq(0,40000,1000), tcl=-0.125, labels=NA)
box()
text(x = relative.axis.point(0.005, "x"),
     y = relative.axis.point(0.925, "y"),
     labels="(B)", font=2, adj=0)

dev.off()

# novelty by human pop density ####

comb.df <- do.call("rbind", all.novel$novel)
comb.df <- comb.df[comb.df$bins == "0",]

comb.df <- merge(comb.df, site.df[,c("site", "long", "lat", "REGION")],
                 by.x = "site", by.y="site",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)


library(raster)
dens.raster <- raster("./raw.datafiles/gpw_v4_population_density_rev11_2000_1_deg.tif")

world <- countriesLow
submap <- world[world$REGION %in% c("Europe", "North America"),]
subraster <- mask(dens.raster.log1, submap)
subraster <- crop(subraster, extent(-180, 180, 30, 80))
subcoasts <- countriesCoarse

point.coords <- comb.df[,c("long", "lat")]
coordinates(point.coords) <- c("long", "lat")

comb.df$pop.dens <- extract(dens.raster, point.coords)

boxplot(log(comb.df$pop.dens+1) ~ comb.df$novel)

comb.df$pop.densM <- log(comb.df$pop.dens+1)
comb.df$REGION <- as.factor(comb.df$REGION)
dens.model <- glmer(novel ~ pop.densM + (1|REGION), data=droplevels(comb.df), family=binomial)
summary(dens.model)

nov.pred <- data.frame(pop.dens = seq(min(comb.df$pop.dens, na.rm=TRUE),
                                      max(comb.df$pop.dens, na.rm=TRUE), len=200))
nov.pred$pop.densM <- log(nov.pred$pop.dens + 1)

nov.pred <- cbind(nov.pred,
                  mer.ci(dens.model, newdata=nov.pred, sims=99, parallel="multicore",
                         cores=8))
nov.pred$fit <- plogis(nov.pred$fit)
nov.pred$lower <- plogis(nov.pred$lower)
nov.pred$upper <- plogis(nov.pred$upper)

pdf(date.wrap("./plots/modern novel comms by pop density", ".pdf"),
    height=4, width=8, useDingbats=FALSE)

split.screen(rbind(c(0.04,0.99,0.6,0.98),
                   c(0.04,0.99,0.6,0.98),
                   c(0.04,0.6,0.1,0.55),
                   c(0.65,0.85,0.255,0.295)))

screen(1)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
rast.ramp <- colorRampPalette(c("white", "pink", "red", "darkred", "black"))(200)
plot(x=NULL, y=NULL, xlab="", ylab="", axes=FALSE, xaxs="i",
     xlim=c(-180,180), ylim=c(25,80))
plot(subraster, legend = FALSE, col = rast.ramp, add=TRUE)
pars <- par("usr")

screen(2)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlab="", ylab="", xaxs="i",
     xlim=c(-180,180), ylim=c(25,80), axes=FALSE)
plot(coastsCoarse, add=TRUE, col="grey70")
plot(submap, add=TRUE)
with(subdf[!subdf$novel,], points(lat ~ long, pch=21, bg="grey50", cex=0.3, lwd=0.5))
with(subdf[subdf$novel,], points(lat ~ long, pch=21, bg="orange", cex=0.65, lwd=0.5))
# box()

axis(side=1, at=seq(-180,180,30), mgp=c(3,0,0),
     labels=parse(text=paste0(seq(-180,180,30), "*degree")))
axis(side=1, at=seq(-180,180,10), mgp=c(3,0,0), labels=NA, tcl=-0.125)

axis(side=2, at=seq(-90,90,30),
     labels=parse(text=paste0(seq(-90,90,30), "*degree")))
axis(side=2, at=seq(-90,90,10), labels=NA, tcl=-0.125)

text(x=relative.axis.point(0.015, "x"),
     y=relative.axis.point(0.955, "y"),
     labels="(A)", font=2)
box()
close.screen(1)

screen(3)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(x=NULL, y=NULL, ylim=c(0,0.15),
     xlim=range(comb.df$pop.densM, na.rm=TRUE),
     axes=FALSE, xlab="", ylab="")

boxplot(comb.df$pop.densM ~ comb.df$novel, add=TRUE, horizontal=TRUE, 
        at=relative.axis.point(c(0.095, 0.905), "y"),
        boxwex=0.02, axes=FALSE)

novel.col <- col2rgb("orange")/255

polygon(x=c(nov.pred$pop.densM, rev(nov.pred$pop.densM)),
        y=c(nov.pred$upper, rev(nov.pred$lower)), 
        col=rgb(novel.col[1], novel.col[2], novel.col[3], 0.3), border=NA)
lines(nov.pred$fit ~ nov.pred$pop.densM, lwd=2, col="orange")

axis(side=2)
axis(side=2, at=seq(0,0.2,0.01), labels=NA, tcl=-0.125)

axis(side = 1,
     at = log(c(0,10,100,1000) + 1),
     labels=c(0,10,100,1000), mgp=c(3,0,0))

axis(side = 1,
     at = log(c(seq(1,10,1),
                seq(10,100,10),
                seq(100,2000,100),
                seq(1000,5000,1000))+1), 
     labels=NA, tcl = -0.125)

mtext(side=1, line=1, adj=0.5,
      text=expression("Population density (people km"^-2*")"))

mtext(side=2, line=2, text="Novel community emergence probability", las=0)

text(y = relative.axis.point(c(0.195, 0.805), "y"),
     x = c(median(subdf$pop.densM[!subdf$novel], na.rm=TRUE),
           median(subdf$pop.densM[subdf$novel], na.rm=TRUE)),
     adj=0.5,
     labels=c("Non-novel", "Novel"), font=2)

text(x=relative.axis.point(0.035, "x"),
     y=relative.axis.point(0.85, "y"),
     labels="(B)", font=2)
box()
close.screen(3)

screen(4)
rast.seq <- seq(cellStats(subraster, "min"),
                cellStats(subraster, "max"), len=200)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=c(0,max(rast.seq)), ylim=c(0,1), 
     xaxs="i", yaxs="i", xlab="", ylab="", axes=FALSE)
box()

image(x=rast.seq,
      y=0,
      z=matrix(1:200, ncol=1),
      col=rast.ramp, add=TRUE, useRaster=TRUE)

axis(side = 1, at = log(c(0,10,100,1000,10000)+1), 
     labels=c(0,10,100,1000,10000), mgp=c(3,0,0))

axis(side = 1,
     at = log(c(seq(1,10,1),
                seq(10,100,10),
                seq(100,2000,100),
                seq(1000,10000,1000))+1), 
     labels=NA, tcl = -0.125)

mtext(side=1, line = 0.75, 
      las=0, adj=0.5,
      text=expression("Population density (people km"^-2*")"))
close.screen(4)
close.screen(all.screens=TRUE)
dev.off()


# test for 'pull of the recent effects' ####

plant.with.genus <- read.csv(paste0(large_file_path, "/processed_genus_records.csv"))

# divide pollen into pre or post last 1000 years
recent.divide <- cut(plant.with.genus$age, breaks = c(0,1000,1e6))

table(recent.divide)

temp <- tapply(plant.with.genus$family,
               recent.divide,
               table)

temp <- lapply(temp, function(x){x[x>0]})

names(temp[[2]])[!names(temp[[2]]) %in% names(temp[[1]])]

recentonly <- names(temp[[1]])[!names(temp[[1]]) %in% names(temp[[2]])]

sum(sort(table(plant.with.genus$family[plant.with.genus$family %in% recentonly]),
         decreasing=TRUE))

plant.old <- plant.with.genus[!plant.with.genus$family %in% recentonly,]

plant.old.slimmed <- plant.old[complete.cases(plant.old[,c("family", "site", "age", "count")]),
                               c("family", "site", "age", "count")]

old.novel <- neotoma.novelty(dataset = plant.old.slimmed,
                             ssmat.type = "abund",
                             bins = seq(-100,max(plant.old.slimmed$age), 200),
                             rich.cutoff = c(100, 5000),
                             age.limits = c(-150, Inf),
                             taxon.res = "family",
                             bin.cutoff = 10,
                             taxa.cutoff = 5,
                             novel.alpha = 0.05,
                             novel.metric = "bray",
                             sqrt.mat=TRUE)

saveRDS(old.novel, 
        "./outputs/all neotoma novelty (pull of the recent).rds")
#old.novel <- readRDS("./outputs/all neotoma novelty (pull of the recent).rds")

# compare old versus all novelty
all.comb <- do.call("rbind", all.novel$novel)
all.comb$siteBin <- paste0(all.comb$site, ".", all.comb$bins)
old.comb <- do.call("rbind", old.novel$novel)
old.comb$siteBin <- paste0(old.comb$site, ".", old.comb$bins)

comb <- merge(all.comb[,c("siteBin", "cat")],
              old.comb[,c("siteBin", "cat")],
              by.x="siteBin", by.y="siteBin",
              all.x=TRUE, all.y=TRUE, sort=FALSE)

table(comb$cat.x, comb$cat.y)
table(comb$cat.x)
table(comb$cat.y)

comb[comb$cat.x=="novel" & comb$cat.y!="novel",]

old.nprob.models <- novel.prob.models(novel.list = old.novel,
                                      site.df = site.df,
                                      name = "potr",
                                      time.k = 40,
                                      test.model=FALSE,
                                      time.age.limits = c(0,25000),
                                      factor.age.limits = c(0,1000),
                                      sauto.n = 5000,
                                      sauto.iter = 999)

novel.prob.plot(prob.model.list = old.nprob.models, 
                env.data = temp.data,
                mod.env.data = mod.env.data,
                ylims=c(0,0.065), 
                regylims=c(0,0.09),
                name = "potr",
                time.age.limits = c(0,25000),
                factor.age.limits = c(0,1000),
                group.letters=c("A","AB","B","B","B","B"))

saveRDS()


# compare varying alphas ####

thresholds <- c(seq(0.01,0.1,0.005))

plant.record.df <- read.csv(paste0(large_file_path, "/processed_neotoma_records.csv"))
plant.record.df <- droplevels(plant.record.df[plant.record.df$dataset.type == "pollen",])

length(unique(plant.record.df$site.id))

summary(is.na(plant.record.df$species))
summary(is.na(plant.record.df$genus))
summary(is.na(plant.record.df$family))

# subset pollen counts at family or higher taxonomic resolution
plant.with.genus <- plant.record.df[!is.na(plant.record.df$genus),]

# next we need to examine chronology
summary(!is.na(plant.with.genus$age))

#~10,000 records with no taxonomy
plant.with.genus <- plant.with.genus[!is.na(plant.with.genus$age), ]

plant.with.genus$site = plant.with.genus$site.id
plant.with.genus <- plant.with.genus[order(plant.with.genus$site),]

plant.with.genus <- droplevels(plant.with.genus)
plant.genus.slimmed <- plant.with.genus[,colnames(plant.with.genus) %in%
                                          c("genus", "site", "age", "count")]

threshold.list <- lapply(thresholds, function(cutoff){
  
  print(paste0("------- ", cutoff))
  
  temp.novel <- neotoma.novelty(dataset = plant.genus.slimmed,
                                ssmat.type = "abund",
                                bins = seq(-100,max(plant.genus.slimmed$age), 200),
                                rich.cutoff = c(10, 10000),
                                age.limits = c(-150, Inf),
                                taxon.res = "genus",
                                bin.cutoff = 10,
                                taxa.cutoff = 2,
                                novel.alpha = cutoff,
                                novel.metric = "jaccard")
  
  temp.novel$novel <- cut.novel(temp.novel$novel, 5)
  
  return(do.call("rbind", temp.novel$novel))
  
})

saveRDS(threshold.list, "./outputs/alpha_threshold_data_list.rds")

threshold.prob.list <- lapply(threshold.list, function(x){
  
  novel.ts.freq <- do.call("rbind", lapply(unique(x$site),
                                           function(site){
                                             
                                             temp <- x[x$site == site,]
                                             
                                             temp <- data.frame(site = site,
                                                                t(sapply(c("back", "instant", "cumul", "novel"),
                                                                         function(y){sum(temp$cat == y)})))
                                             temp$lat <- x$data$lat[x$site == site][1]
                                             temp$long <- x$data$long[x$site == site][1]
                                             
                                             return(temp)
                                           }))
  novel.ts.freq <- as.data.frame(novel.ts.freq)
  
  # modelling the probability of each classification occurring.
  novel.ts.freq$non.novel <- rowSums(novel.ts.freq[,c("back", "cumul", "instant")])
  novel.ts.freq$non.back <- rowSums(novel.ts.freq[,c("novel", "cumul", "instant")])
  novel.ts.freq$non.instant <- rowSums(novel.ts.freq[,c("novel", "cumul", "back")])
  novel.ts.freq$non.cumul <- rowSums(novel.ts.freq[,c("novel", "back", "instant")])
  
  ## modeling the probability of our two GAM tests.
  novel.ts.freq$all.instant <- rowSums(novel.ts.freq[,c("novel", "instant")])
  novel.ts.freq$all.cumul <- rowSums(novel.ts.freq[,c("novel", "cumul")])
  novel.ts.freq$non.all.instant <- rowSums(novel.ts.freq[,c("cumul", "back")])
  novel.ts.freq$non.all.cumul <- rowSums(novel.ts.freq[,c("instant", "back")])
  
  require(multcomp)
  require(lme4)
  require(merTools)
  
  prob.models <- lapply(1:6, function(n){
    
    success.var = c("all.instant", "all.cumul", "back", "instant", "cumul", "novel")[n]
    failure.var = c("non.all.instant", "non.all.cumul", "non.back", 
                    "non.instant", "non.cumul", "non.novel")[n]
    
    print(success.var)
    
    temp.df <- novel.ts.freq
    temp.df$success = temp.df[, success.var]
    temp.df$failure = temp.df[, failure.var]
    
    taxa.prob.m <- glm(cbind(success, failure) ~ 1, data=temp.df, family=binomial)
    
    # group predictions
    pred.df <- summary(taxa.prob.m)$coefficients
    
    return(list(model=taxa.prob.m,
                pred.df = pred.df))
    
  })
  
  return(list(data = novel.ts.freq,
              prob.models = prob.models))
  
})

saveRDS(threshold.prob.list, "./outputs/alpha_threshold_prob_list.rds")

threshold.prob.list <- readRDS("./outputs/alpha_threshold_prob_list.rds")

varying.alpha.test(probs = thresholds,
                   threshold.prob.list = threshold.prob.list,
                   circle.1 = 0.01,
                   circle.2 = 0.05,
                   circle.3 = 0.1,
                   ylims = c(0.005,max(thresholds)+0.0035))

# Bray vs chord^2 ####

# create square-chord novelty and probability models

plant.with.genus <- droplevels(plant.with.genus)
plant.genus.slimmed <- plant.with.genus[,colnames(plant.with.genus) %in%
                                          c("family", "site", "age", "count")]

chord.nov <- neotoma.novelty(dataset = plant.genus.slimmed,
                             ssmat.type = "abund",
                             bins = seq(-100,max(plant.genus.slimmed$age), 200),
                             rich.cutoff = c(100, 10000),
                             age.limits = c(-150, Inf),
                             taxon.res = "family",
                             bin.cutoff = 10,
                             taxa.cutoff = 5,
                             novel.alpha = 0.05,
                             novel.metric = "SQchord",
                             sqrt.mat=TRUE)

saveRDS(chord.nov, 
        "./outputs/chord neotoma novelty (sub-sampled).rds")

chord.nov.df <- do.call('rbind', chord.nov$novel)

table(chord.nov.df$cat[as.numeric(chord.nov.df$bins) <= 25000])

# now model chord-squared novel probability over time

# how many novel communities did framework detect?

chord.nprob.models <- novel.prob.models(novel.list = chord.nov,
                                        site.df = site.df,
                                        time.k = 40,
                                        test.model=FALSE,
                                        name = "chord",
                                        time.age.limits = c(1200,25000),
                                        factor.age.limits = c(0,1000),
                                        sauto.n = 5000,
                                        sauto.iter = 999)

novel.prob.plot(prob.model.list = chord.nprob.models, 
                env.data = temp.data,
                mod.env.data = mod.env.data,
                ylims=c(0,0.022), 
                regylims=c(0,0.022),
                name = "chord",
                time.age.limits = c(0,25000),
                factor.age.limits = c(0,1000),
                group.letters=c("A","AB","B","B","B","B"))

saveRDS(chord.nprob.models,
        date.wrap("./outputs/novel probability models (chord)", ".rds"))

# Burke et al Taxa novelty ####

plant.genus.slimmed <- plant.with.genus[,colnames(plant.with.genus) %in%
                                          c("burkeTaxa", "site", "age", "count")]

burke.nov <- neotoma.novelty(dataset = plant.genus.slimmed,
                             ssmat.type = "abund",
                             bins = seq(-100,max(plant.genus.slimmed$age, na.rm=TRUE), 200),
                             rich.cutoff = c(100, 10000),
                             age.limits = c(-150, Inf),
                             taxon.res = "burkeTaxa",
                             bin.cutoff = 10,
                             taxa.cutoff = 5,
                             novel.alpha = 0.05,
                             novel.metric = "bray",
                             sqrt.mat=TRUE)

saveRDS(burke.nov, 
        "./outputs/Burke neotoma novelty (sub-sampled).rds")

burke.nov.df <- do.call('rbind', burke.nov$novel)

table(burke.nov.df$cat[as.numeric(burke.nov.df$bins) <= 25000])

# now model chord-squared novel probability over time

# how many novel communities did framework detect?
burke.nprob.models <- novel.prob.models(novel.list = burke.nov,
                                        site.df = site.df,
                                        time.k = 40,
                                        test.model=FALSE,
                                        name = "burke",
                                        time.age.limits = c(1200,25000),
                                        factor.age.limits = c(0,1000),
                                        sauto.n = 5000,
                                        sauto.iter = 999)

novel.prob.plot(prob.model.list = burke.nprob.models, 
                env.data = temp.data,
                mod.env.data = mod.env.data,
                ylims=c(0,0.065), 
                regylims=c(0,0.065),
                name = "burke",
                time.age.limits = c(0,25000),
                factor.age.limits = c(0,1000),
                group.letters=c("A","AB","B","B","B","B"))

saveRDS(burke.nprob.models,
        date.wrap("./outputs/novel probability models (burke)", ".rds"))

# Adding some families into genus-level analyses ####

# find which families are most often missing genus-level data
plant.with.family <- read.csv(paste0(large_file_path, "/processed_family_records.csv"))
missingGen <- table(plant.with.family$family[is.na(plant.with.family$genus)])

# include families that have >0.1% records
sort(missingGen, decreasing=TRUE)/nrow(plant.with.family)

famToInclude <- names(missingGen[(missingGen / nrow(plant.with.family)) > 1e-3])

# now override the NA genus IDs for these families with a family-level Id so they
# are included in our novelty detection analyses
plant.with.family$genus <- as.character(plant.with.family$genus)
plant.with.family$genus[plant.with.family$family %in% famToInclude &
                          is.na(plant.with.family$genus)] = as.character(plant.with.family$family[plant.with.family$family %in% famToInclude &
                                                                                                    is.na(plant.with.family$genus)])

plant.family.slimmed <- plant.with.family[,colnames(plant.with.family) %in%
                                            c("genus", "site", "age", "count")]
all.novel.subsamp <- neotoma.novelty(dataset = plant.family.slimmed,
                                     ssmat.type = "abund",
                                     bins = seq(-100,max(plant.family.slimmed$age), 200),
                                     rich.cutoff = c(100, 10000),
                                     age.limits = c(-150, Inf),
                                     taxon.res = "genus",
                                     bin.cutoff = 10,
                                     taxa.cutoff = 5,
                                     novel.alpha = 0.05,
                                     novel.metric = "bray",
                                     sqrt.mat=TRUE)

famgenhybrid.model <- novel.prob.models(novel.list = all.novel.subsamp,
                                        site.df = site.df,
                                        time.k = 40,
                                        test.model=FALSE,
                                        name = "famgen-hybrid",
                                        time.age.limits = c(1200,25000),
                                        factor.age.limits = c(0,1000),
                                        sauto.n = 5000,
                                        sauto.iter = 999)

novel.prob.plot(prob.model.list = famgenhybrid.model, 
                env.data = temp.data,
                mod.env.data = mod.env.data,
                ylims=c(0,0.065), 
                regylims=c(0,0.09),
                name = "famgen-hybrid",
                time.age.limits = c(0,25000),
                factor.age.limits = c(0,1000),
                group.letters=c("A","AB","B","B","B","B"))

famgen.df <- do.call('rbind', all.novel.subsamp$novel)
with(famgen.df[as.numeric(as.character(famgen.df$bins)) <= 25000,],
     table(cat))

all.df <- do.call('rbind', all.novel$novel)
with(all.df[as.numeric(as.character(all.df$bins)) <= 25000,],
     table(cat))

all.base <- all.df[,c("site", "bins", "cat")]
all.base$siteID <- paste0(all.base$site, all.base$bins)

famgen.df$siteID <- paste0(famgen.df$site, famgen.df$bins)


saveRDS(famgenhybrid.model,
        date.wrap("./outputs/novel probability models (famgen-hybrid)", ".rds"))


test <- merge(plant.with.genus, site.df[,c("site", "REGION")],
              by.x="site", by.y="site", all.x=TRUE, all.y=FALSE, sort=FALSE)
head(test)
table(test$taxon == "Cyperaceae", test$REGION)

# pollen sampling over time ####

plantSplit <- split(plant.with.genus,
                    f=paste0(plant.with.genus$site.id, ".",
                             plant.with.genus$sample.id))

pollen.sampling <- do.call("rbind", lapply(plantSplit,
                                           function(x){
                                             
                                             genSum <- tapply(x$count, factor(!is.na(x$genus), 
                                                                              levels=c("FALSE", "TRUE")), 
                                                              sum, na.rm=TRUE)
                                             genSum[is.na(genSum)] = 0
                                             genProp <- genSum[2] / sum(genSum)
                                             return(data.frame(genProp = genProp,
                                                               age = x$age[1],
                                                               site = x$site[1]))
                                           }))

# add in continent
pollen.sampling <- merge(x=pollen.sampling, y=site.df[,c("site.id", "REGION")],
                         by.x="site", by.y="site.id",)

# model pollen genus fraction over time for each region
pollenGam <- gam(genProp ~ s(age, by=REGION) + REGION, family=betar, data=pollen.sampling)

pred.df <- data.frame(age=rep(seq(0,25000,len=200),2),
                      REGION=rep(c("Europe", "North America"), each=200))
pred.df <- cbind(pred.df,
                 as.data.frame(predict(pollenGam, newdata=pred.df, se.fit=TRUE)))

pollTab <- table(cut(pollen.sampling$genProp, breaks=seq(0,1,0.1)),
                 cut(pollen.sampling$age, breaks=seq(0, 25000, 5000)))
pollTab <- apply(pollTab, 2, function(x){x / sum(x)})

pdf("./plots/pollen genus patterns.pdf", height=7, width=6, useDingbats = FALSE)
par(mfrow=c(2,1), mar=c(3,0,0,0), oma=c(1,4,1,1), ps=10, tcl=-0.25, las=1, mgp=c(3,0.5,0))

plot(x=NULL, y=NULL, xlim=c(25000,0), ylim=c(0,1), xaxt="n")
with(pred.df[pred.df$REGION=="Europe",],
     {
       polygon(x=c(age, rev(age)),
               y=plogis(c(fit + 1.96 * se.fit, rev(fit - 1.96*se.fit))),
               border=NA, col=rgb(0,0,1,0.5))
       lines(plogis(fit) ~ age, lwd=2, col="blue")
     })

dgrgb <- col2rgb("darkgreen")/255
with(pred.df[pred.df$REGION=="North America",],
     {
       polygon(x=c(age, rev(age)),
               y=plogis(c(fit + 1.96 * se.fit, rev(fit - 1.96*se.fit))),
               border=NA, col=rgb(dgrgb[1],dgrgb[2],dgrgb[3],0.5))
       lines(plogis(fit) ~ age, lwd=2, col="darkgreen")
     })

mtext(side=2, line=2, las=0, text="Fraction of pollen with genus-level identification")
axis(side=1, mgp=c(3,0.2,0))
mtext(side=1, line=1.5, las=0, text="Age (years before 1950AD)")

text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(A)", font=2)

text(x=relative.axis.point(0.215, "x"),
     y=relative.axis.point(0.95, "y"),
     labels=expression(bold("North America"*phantom(" & Europe"))), 
     font=2, col="#629D26")
text(x=relative.axis.point(0.215, "x"),
     y=relative.axis.point(0.95, "y"),
     labels=expression(bold(phantom("North America")*" & "*phantom("Europe"))), 
     font=2, col="black")
text(x=relative.axis.point(0.215, "x"),
     y=relative.axis.point(0.95, "y"),
     labels=expression(bold(phantom("North America & ")*"Europe")), 
     font=2, col="blue")

barWidth = 0.015
plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,0.6), yaxs="i", xaxt="n", xlab="", ylab="",
     xaxs="i")
rect(xleft=rep(seq(0.05,0.95,0.1), 4) + rep(seq(-(2*barWidth),(2*barWidth),len=5), each=10)- 0.5*barWidth,
     xright=rep(seq(0.05,0.95,0.1), 4) + rep(seq(-(2*barWidth),(2*barWidth),len=5), each=10) + 0.5*barWidth,
     ybottom=0, ytop=pollTab,
     col=rep(colorRampPalette(c("black", "white"))(5), each=10))
axis(side=1, at=seq(0.05,0.95,0.1), labels=NA)
mtext(side=2, line=2, las=0, text="Proportion of pollen samples")
mtext(side=1, line=2.5, las=0, text="Fraction of pollen with genus-level identification")
par(xpd=NA)
text(x=seq(0.05,0.95,0.1),
     y=relative.axis.point(-0.035, "y"),
     labels=paste0(seq(0,90,10),
                   "-",
                   seq(10,100,10), "%"), srt=30, adj=1)
text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(B)", font=2)

legend(x=0, y=0.5, fill=colorRampPalette(c("black", "white"))(5),
       legend=c("0-5,000 ybp",
                "5,001-10,000 ybp",
                "10,001-15,000 ybp",
                "15,001-20,000 ybp",
                "20,000-25,000 ybp"), bty="n",
       pt.cex=1.25, y.intersp=0.8, x.intersp=0.7)

dev.off()
# Sampling map ####

library(sp)
library(maps)
library(maptools)
library(raster)
library(rworldmap)

# fix some of the annoying regions
world <- countriesLow
allWorld <- crop(x=world, y=extent(-180,180,20,90))
world <- allWorld[allWorld$REGION %in% c("Europe", "North America"),]

pdf("./plots/neotomaMap.pdf", height=2.4, width=8)

split.screen(rbind(c(0.05,0.95,0.02,0.99)))

screen(1)
par(mar=c(2,0,0,2), ps=8, tcl=-0.25)
plot(allWorld, yaxs="i", col="grey90", border="grey90", xaxs="i")
plot(world, add=TRUE, col=c("red","white", "goldenrod", "brown", "blue", "#629D26", "purple")[world$REGION], 
     border=c("red","white", "goldenrod", "brown", "blue", "#629D26", "purple")[world$REGION])
points(site.df$lat ~ site.df$long, pch=21, cex=0.3, add=TRUE, bg="white", lwd=0.5)

par(xpd=NA)
latCut <- table(cut(site.df$lat, breaks=seq(20,80,5)))
latCut <- latCut / sum(latCut)

axis(side=2, at=seq(20,90,10), mgp=c(3,0.35,0),
     labels=paste0(seq(20,90,10), "°"), las=1)

axis(side=1, at=seq(-180,180,30), mgp=c(3,0.35,0),
     labels=paste0(seq(-180,180,30), "°"), las=1, mgp=c(3,0.1,0))

par(xpd=NA)
rect(xleft=181.5, xright= 181.5 + 50*latCut,
     ybottom=seq(20,75,5),
     ytop=seq(25,80,5), col="grey")
text(x= 181.5 + 50*latCut, y=seq(22.5,77.5,5),
     labels=paste0(round(latCut*100, 1),"%"), pos=4)
par(xpd=FALSE)
box()
close.screen(1)

dev.off()
