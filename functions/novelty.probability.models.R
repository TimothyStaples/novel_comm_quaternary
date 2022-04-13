novel.prob.models <- function(novel.list,
                              test.model=FALSE,
                              time.age.limits,
                              factor.age.limits,
                              name,
                              site.df,
                              time.k,
                              sauto.n = 100,
                              sauto.iter = 5){

require(gamm4)
require(DHARMa)
require(merTools)
require(multcomp)
require(lmtest)

print("Preparing data")  
comb.df <- do.call("rbind", novel.list$novel)

comb.df$bins <- as.numeric(as.character(comb.df$bins))

# add in site-level info
comb.df <- droplevels(merge(comb.df, site.df[,c("REGION", "site", "elev", "long", "lat")],
                 by.x="site", by.y="site", all.x=TRUE, all.y=FALSE, sort=FALSE))

comb.df$REGION <- as.factor(comb.df$REGION)

# drop sites < -10 msl
comb.df <- comb.df[comb.df$elev > -10, ]

tsRichness <- sapply(novel.list$prop.ssmats, ncol)
tsRichness <- data.frame(site=names(tsRichness),
                         tsRichness = tsRichness,
                         tsLength = sapply(novel.list$prop.ssmats, nrow))

comb.df <- merge(comb.df, tsRichness,
                 by.x="site", by.y="site",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)

# get bin position within each time-series
comb.df <- do.call("rbind", lapply(split(comb.df, f=comb.df$site), function(x){
         bin.raw <- as.numeric(as.character(x$bins))
         bin.sort <- sort(as.numeric(as.character(x$bins)))
         return(cbind(x, bin.n = nrow(x)+1 - match(bin.raw, bin.sort)))
         }))

comb.df.full <- comb.df


# remove time points outside of time.age.limits
comb.df <- comb.df[comb.df$bins >= time.age.limits[1] &
                     comb.df$bins <= time.age.limits[2],]

comb.df$bin.lag.scale = as.vector(scale(log(comb.df$bin.lag)))
comb.df$bin.n.scale = as.vector(scale(log(comb.df$bin.n)))
comb.df$tsLength.scale = as.vector(scale(log(comb.df$tsLength)))
comb.df$tsRichness.scale = as.vector(scale(log(comb.df$tsRichness)))

# ln-transform elevation, adding +11 so -10 msl scores are positive
comb.df$elev.scale <- as.vector(scale(log(comb.df$elev + 11)))

comb.df <- comb.df[complete.cases(comb.df$REGION),]

# Time model ####

#                         model ####
print("Running continuous time model")
time.m <- gam(novel ~ s(bins, k = time.k) + 
                bin.lag.scale + bin.n.scale + tsLength.scale + 
                tsRichness.scale + elev.scale,
                family=binomial, data=comb.df)

write.csv(summary(time.m)$s.table,
          date.wrap(paste0("./outputs/time model novel smooth coefs (", name, ")"), ".csv"))

write.csv(summary(time.m)$p.table,
          date.wrap(paste0("./outputs/time model novel parametric coefs (", name, ")"), ".csv"))

#                   diagnostics ####

if(test.model){
print("Running model diagnostics")

  timeDiag <- modelDiagTests(model = time.m,
                 time=time.m$model$bins,
                 data = comb.df,
                 sautoSubsample = TRUE,
                 spat.iter = sauto.iter,
                 spat.sub.size = sauto.n)

}

#                   predictions ####

pred.df <- data.frame(bins=seq(min(comb.df$bins, na.rm=TRUE), 
                               max(comb.df$bins, na.rm=TRUE), 200),
                      bin.lag.scale = 0,
                      bin.n.scale = 0,
                      tsLength.scale = 0,
                      tsRichness.scale = 0,
                      elev.scale=0)

pred.df <- cbind(pred.df,
                 as.data.frame(predict(time.m, newdata=pred.df, se.fit=TRUE)))

write.csv(pred.df, date.wrap(paste0("./outputs/novel comms through time model coefs (",
                                        name, ")"), ".csv"))


# Time REGION model ####

if(length(unique(comb.df$REGION)) > 1){
#                         model ####
print("Running continuous time (REGION) model")
time.reg.m <- gam(novel ~ s(bins, k = time.k, by=REGION) + REGION + 
                    bin.lag.scale + bin.n.scale + tsLength.scale + tsRichness.scale + elev.scale,
                    family=binomial, data=comb.df)

write.csv(summary(time.reg.m)$s.table,
          date.wrap(paste0("./outputs/time REGION model novel smooth coefs (", name, ")"), ".csv"))

write.csv(summary(time.reg.m)$p.table,
          date.wrap(paste0("./outputs/time REGION model novel parametric coefs (", name, ")"), ".csv"))

#                   diagnostics ####

if(test.model){
  print("Running model diagnostics")
  
  timeRegionDiag <- modelDiagTests(model = time.reg.m,
                             time=time.reg.m$model$bins,
                             data = comb.df,
                             sautoSubsample = TRUE,
                             spat.iter = sauto.iter,
                             spat.sub.size = sauto.n)
  
}

#                   predictions ####

reg.pred.list <- lapply(1:length(levels(comb.df$REGION)),
                        function(n){
print(n)                          
pred.df <- data.frame(bins=seq(min(comb.df$bins, na.rm=TRUE),
                               max(comb.df$bins, na.rm=TRUE), 200),
                      bin.lag.scale = 0,
                      bin.n.scale = 0,
                      tsLength.scale = 0,
                      tsRichness.scale=0,
                      elev.scale=0)

pred.df$REGION <- levels(comb.df$REGION)[n]

regRange <- range(as.numeric(as.character(comb.df$bins[comb.df$REGION==levels(comb.df$REGION)[n]])),
                  na.rm=TRUE)
pred.df <- pred.df[pred.df$bins >= regRange[1] &
                   pred.df$bins <= regRange[2],]

return(cbind(pred.df,
             as.data.frame(predict(time.reg.m, newdata=pred.df, se.fit=TRUE))))
                          
                        })

} else {
  time.reg.m = NULL
  reg.pred.list=NULL
}

# Factor model ####

comb.df.fact = comb.df.full[comb.df.full$bins >= factor.age.limits[1] &
                              comb.df.full$bins <= factor.age.limits[2],]
comb.df.fact$bin.fact <- as.factor(comb.df.fact$bins)
comb.df.fact$bin.lag.scale = as.vector(scale(log(comb.df.fact$bin.lag)))
comb.df.fact$bin.n.scale = as.vector(scale(log(comb.df.fact$bin.n)))
comb.df.fact$tsLength.scale = as.vector(scale(log(comb.df.fact$tsLength)))
comb.df.fact$tsRichness.scale = as.vector(scale(log(comb.df.fact$tsRichness)))
comb.df.fact$elev.scale <- as.vector(scale(log(comb.df.fact$elev + 11)))

comb.df.fact <- comb.df.fact[!is.na(comb.df.fact$novel),]

if(length(unique(comb.df.fact$bin.fact)) > 1){
print("Running factor model...")

time.fact.m <- glm(novel ~ bin.fact + bin.lag.scale + bin.n.scale + 
                       tsLength.scale + tsRichness.scale + elev.scale,
                     family=binomial(), data=comb.df.fact) 
  
write.csv(summary(time.fact.m)$coefficients, 
          date.wrap(paste0("./outputs/time factor novel coefs (", name, ")"), ".csv"))

if(test.model){
  print("Running post-model tests...")
  
  FactDiag <- modelDiagTests(model = time.fact.m,
                                   time=as.numeric(as.character(time.fact.m$model$bin.fact)),
                                   data = comb.df.fact,
                                   spat.iter = sauto.iter,
                                   spat.sub.size = sauto.n)
}

print("Calculating CIs")

fact.pred.df <- data.frame(bin.fact=factor(levels(comb.df.fact$bin.fact),
                                           levels = levels(comb.df.fact$bin.fact)),
                           bin.lag.scale = 0,
                           bin.n.scale = 0,
                           tsLength.scale = 0,
                           tsRichness.scale = 0,
                           REGION="A",
                           elev.scale=0)

fact.pred.df <- cbind(fact.pred.df,
                 as.data.frame(predict(time.fact.m,
                                   newdata=fact.pred.df,
                                   se.fit=TRUE)))

fact.pred.df$upper <- plogis(fact.pred.df$fit + 1.96 * fact.pred.df$se.fit)
fact.pred.df$lower <- plogis(fact.pred.df$fit - 1.96 * fact.pred.df$se.fit)
fact.pred.df$fit = plogis(fact.pred.df$fit)
fact.pred.df$bin.num <- as.numeric(as.character(fact.pred.df$bin.fact))

} else {
  time.fact.m = NULL
  fact.pred.df = NULL
  }

# Factor REGION model ####

if(length(unique(comb.df.fact$bin.fact)) > 1 &
   length(unique(comb.df.fact$REGION)) > 1){
  print("Running factor region model...")
  
  time.fact.reg.m <- glm(novel ~ bin.fact * REGION + 
                        bin.lag.scale + bin.n.scale + tsLength.scale + tsRichness.scale + elev.scale,
                       family=binomial(), data=comb.df.fact) 
  
  write.csv(summary(time.fact.reg.m)$coefficients, 
            date.wrap(paste0("./outputs/time factor region novel coefs (", name, ")"), ".csv"))
  
  if(test.model){
    print("Running post-model tests...")
    FactRegDiag <- modelDiagTests(model = time.fact.reg.m,
                               time=as.numeric(as.character(time.fact.reg.m$model$bin.fact)),
                               data = comb.df.fact,
                               spat.iter = sauto.iter,
                               spat.sub.size = sauto.n)
  }
  
  print("Calculating CIs")
  fact.reg.pred.list <- do.call("rbind", lapply(1:length(levels(comb.df$REGION)),
                          function(n){
                            
                            fact.pred.df <- data.frame(bin.fact=factor(levels(comb.df.fact$bin.fact),
                                                                       levels = levels(comb.df.fact$bin.fact)),
                                                       bin.lag.scale = 0,
                                                       bin.n.scale = 0,
                                                       tsLength.scale = 0,
                                                       tsRichness.scale = 0,
                                                       REGION=levels(comb.df$REGION)[n],
                                                       elev.scale=0)
                            
                            return(cbind(fact.pred.df,
                                         as.data.frame(predict(time.fact.reg.m,
                                                                    newdata=fact.pred.df,
                                                                    re.form=NA,
                                                                    se.fit=TRUE))))
                            
                          }))
  
  fact.reg.pred.list$upper <- plogis(fact.reg.pred.list$fit + 1.96 * fact.reg.pred.list$se.fit)
  fact.reg.pred.list$lower <- plogis(fact.reg.pred.list$fit - 1.96 * fact.reg.pred.list$se.fit)
  fact.reg.pred.list$fit = plogis(fact.reg.pred.list$fit)
  fact.reg.pred.list$bin.num <- as.numeric(as.character(fact.reg.pred.list$bin.fact))
  
} else {
  time.fact.reg.m = NULL
  fact.reg.pred.list = NULL
  }

# Return data ####

return.list <- list(time.data = comb.df,
                    fact.data = comb.df.fact,
                    time.model = time.m,
                    time.reg.model = time.reg.m,
                    fact.model = time.fact.m,
                    fact.reg.model = time.fact.reg.m,
                    time.pred = pred.df,
                    time.reg.pred = reg.pred.list,
                    fact.pred = fact.pred.df,
                    fact.reg.pred = fact.reg.pred.list)

if(test.model){
  return.list <- list(time.data = comb.df,
                      fact.data = comb.df.fact,
                      time.model = time.m,
                      time.reg.model = time.reg.m,
                      fact.model = time.fact.m,
                      fact.reg.model = time.fact.reg.m,
                      time.pred = pred.df,
                      time.reg.pred = reg.pred.list,
                      fact.pred = fact.pred.df,
                      fact.reg.pred = fact.reg.pred.list,
                      timeDiag = timeDiag,
                      timeRegionDiag = timeRegionDiag,
                      FactDiag = FactDiag,
                      FactRegDiag = FactRegDiag)
  
}

return(return.list)

}
