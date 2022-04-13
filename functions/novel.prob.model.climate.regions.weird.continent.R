regions.novel.prob.model.continent <- function(novel.list,
                                     site.df,
                                     all.prob.model,
                                     name,
                                     lat.lims=c(0,90),
				                             time.k,
			 	                             fact.k,
                                     time.age.limits,
                                     factor.age.limits,
                                     test.model=FALSE,
                                     sauto.size = 1000,
                                     sauto.n = 3){

  require(merTools)
  require(DHARMa)
  require(gamm4)
  require(lme4)
  
  # Prepare data ####
  print("prepare data")
  comb.df <- do.call("rbind", novel.list$novel)
  comb.df$bins <- as.numeric(as.character(comb.df$bins))
  
  # get bin position within each time-series
  comb.df <- do.call("rbind", lapply(split(comb.df, f=comb.df$site), function(x){
    bin.raw <- as.numeric(as.character(x$bins))
    bin.sort <- sort(as.numeric(as.character(x$bins)))
    return(cbind(x, bin.n = nrow(x)+1 - match(bin.raw, bin.sort)))
  }))
  
  comb.df.geo <- merge(comb.df, site.df[,c("site", "long", "lat", "REGION", "elev")],
                       by.x="site", by.y="site",
                       all.x=TRUE, all.y=FALSE, sort=FALSE)

  comb.df.geo$abs.lat <- abs(comb.df.geo$lat)
  
  # for this analysis, just run on Nth America & Europe
  comb.df.geo <- droplevels(comb.df.geo[comb.df.geo$REGION %in% c("Europe", "North America"),])
  comb.df.geo$REGION <- as.factor(comb.df.geo$REGION)
  
  comb.df.geo <- comb.df.geo[comb.df.geo$elev > -10 &
                             comb.df.geo$abs.lat >= lat.lims[1] &
                             comb.df.geo$abs.lat <= lat.lims[2],]
  
  tsRichness <- sapply(novel.list$prop.ssmats, ncol)
  tsRichness <- data.frame(site=names(tsRichness),
                           tsRichness = tsRichness,
                           tsLength = sapply(novel.list$prop.ssmats, nrow))
  
  comb.df.geo <- merge(comb.df.geo, tsRichness,
                        by.x="site", by.y="site",
                        all.x=TRUE, all.y=FALSE, sort=FALSE)
  
  comb.df.time <- droplevels(comb.df.geo[comb.df.geo$bins >= time.age.limits[1] &
                                       comb.df.geo$bins <= time.age.limits[2],])
  

  tslength <- table(comb.df.time$site)
  comb.df.time$tslength <- tslength[match(comb.df.time$site, names(tslength))]
  comb.df.time$ts.scale <- as.vector(scale(log(comb.df.time$tslength)))
  comb.df.time$bins <- as.numeric(as.character(comb.df.time$bins))
  comb.df.time$bin.lag.scale = as.vector(scale(log(comb.df.time$bin.lag)))
  comb.df.time$bin.n.scale = as.vector(scale(log(comb.df.time$bin.n)))
  comb.df.time$tsRichness.scale = as.vector(scale(log(comb.df.time$tsRichness)))
  comb.df.time$elev.scale <- as.vector(scale(log(comb.df.time$elev + 11)))
  
  # Run latitude model ####
  print("Running full latitudinal model")
  novel.geo.m <- gam(novel ~ t2(bins, abs.lat, k=time.k, by=REGION) + 
                     REGION + bin.lag.scale + bin.n.scale + ts.scale + tsRichness.scale + elev.scale,
                     family=binomial(), data=comb.df.time)
  
  write.csv(summary(novel.geo.m)$p.table,
            date.wrap(paste0("./outputs/time latitude novel parametric coefs (", name, ")"), ".csv"))
  
  write.csv(summary(novel.geo.m)$s.table,
            date.wrap(paste0("./outputs/time latitude novel smooth coefs (", name, ")"), ".csv"))
  
  # Factor spline model ####
  comb.geo.sub <- droplevels(comb.df.geo[comb.df.geo$bins >= factor.age.limits[1] &
                                           comb.df.geo$bins <= factor.age.limits[2],])
  comb.geo.sub$bin.fact <- as.factor(comb.geo.sub$bins)
  comb.geo.sub$bin.lag.scale = as.vector(scale(log(comb.geo.sub$bin.lag)))
  comb.geo.sub$bins <- as.numeric(as.character(comb.geo.sub$bins))
  comb.geo.sub$tslength <- tslength[match(comb.geo.sub$site, names(tslength))]
  comb.geo.sub$ts.scale <- as.vector(scale(log(comb.geo.sub$tslength)))
  comb.geo.sub$bins <- as.numeric(as.character(comb.geo.sub$bins))
  comb.geo.sub$bin.n.scale = as.vector(scale(log(comb.geo.sub$bin.n)))
  comb.geo.sub$tsRichness.scale = as.vector(scale(log(comb.geo.sub$tsRichness)))
  comb.geo.sub$elev.scale <- as.vector(scale(log(comb.geo.sub$elev + 11)))
  
  # now model over latitude explicitly for each bin
  print("Running last 100 years latitudinal model")
  novel.geo.fact.m <- gam(novel ~ t2(bins, abs.lat, k=fact.k, by=REGION) + 
                            REGION + bin.lag.scale + bin.n.scale + ts.scale + tsRichness.scale + elev.scale,
                            family=binomial(), data=comb.geo.sub)
  
  write.csv(summary(novel.geo.fact.m)$p.table,
            date.wrap(paste0("./outputs/modern time latitude novel parametric coefs (", name, ")"), ".csv"))
  
  write.csv(summary(novel.geo.fact.m)$s.table,
            date.wrap(paste0("./outputs/modern time latitude novel smooth coefs (", name, ")"), ".csv"))
  
  
# Model diagnostics ####
if(test.model){
  print("Running model diagnostics")
  gamm.binom.sim <- function(gam.object, n){
    pred <- predict(gam.object)
    replicate(n, rbinom(length(pred), 1, plogis(pred)))
  }
  
  time.gam.sim <- gamm.binom.sim(novel.geo.m, 250)
  time.gam.res <- createDHARMa(simulatedResponse = time.gam.sim,
                               observedResponse = novel.geo.m$model$novel,
                               fittedPredictedResponse = predict(novel.geo.m, type = "response"),
                               integerResponse = TRUE)
  time.gam.res$refit=FALSE
  time.gam.res$simulatedResponse = time.gam.sim
  
  print("Time model...")
  time.unif.test = testUniformity(time.gam.res)
  print(ifelse(abs(time.unif.test$statistic) > 0.1, 
               "Uniformity not okay :(",
               "Uniformity okay :)"))
  
  time.disp.test = testDispersion(time.gam.res)
  print(ifelse(time.disp.test$p.value <=0.05, 
               "Dispersal not okay :(",
               "Dispersal okay :)"))
  
  time.data <- novel.geo.m$model$bins
  time.offset <- time.data
  while(sum(table(time.offset) > 2) > 0){
    time.offset <- time.data + rnorm(length(time.offset), 0, 1)
  }
  
  time.tauto.test <- testTemporalAutocorrelation(time.gam.res,
                                                 time=time.offset)
  
  print(ifelse(time.tauto.test$p.value <=0.05, 
               "Temporal autocorrelation not okay :(",
               "Temporal autocorrelation okay :)"))
  
  time.sauto.test <- subsample.spat.auto.test.binom.gamm(model = novel.geo.m,
                                                         coords = comb.df.geo[,c("long","lat")],
                                                         n = sauto.size,
                                                         iter = sauto.n)
  print(ifelse(mean(time.sauto.test$p) <=0.05, 
               "Spatial autocorrelation not okay :(",
               "Spatial autocorrelation okay :)"))
  
  print("Recent model...")
  
  fact.gam.sim <- gamm.binom.sim(novel.geo.fact.m, 250)
  fact.gam.res <- createDHARMa(simulatedResponse = fact.gam.sim,
                               observedResponse = novel.geo.fact.m$model$novel,
                               fittedPredictedResponse = predict(novel.geo.fact.m, type = "response"),
                               integerResponse = TRUE)
  fact.gam.res$refit=FALSE
  fact.gam.res$simulatedResponse = fact.gam.sim
  
  fact.unif.test = testUniformity(fact.gam.res)
  print(ifelse(abs(fact.unif.test$statistic) > 0.1, 
               "Uniformity not okay :(",
               "Uniformity okay :)"))
  
  fact.disp.test = testDispersion(fact.gam.res)
  print(ifelse(fact.disp.test$p.value <=0.05, 
               "Dispersal not okay :(",
               "Dispersal okay :)"))
  
  time.data <- novel.geo.fact.m$model$bins
  time.offset <- time.data
  while(sum(table(time.offset) > 2) > 0){
    time.offset <- time.data + rnorm(length(time.offset), 0, 1)
  }
  
  fact.tauto.test <- testTemporalAutocorrelation(fact.gam.res,
                                                 time=time.offset)
  
  print(ifelse(fact.tauto.test$p.value <=0.05, 
               "Temporal autocorrelation not okay :(",
               "Temporal autocorrelation okay :)"))

  fact.sauto.test <- gc.spat.auto.test.binom.gamm4(model = novel.geo.fact.m,
                                                   data.sub = comb.geo.sub)
  
  print(ifelse(fact.sauto.test$p <=0.05, 
               "Spatial autocorrelation not okay :(",
               "Spatial autocorrelation okay :)"))
  
}

# Create return list ####
return.list = list(data = comb.df.geo,
                   time.model = novel.geo.m,
                   fact.model = novel.geo.fact.m)

if(test.model){
    return.list = list(data = comb.df.geo,
                       time.model = novel.geo.m,
                       fact.model = novel.geo.fact.m,
                       time.unif.test = time.unif.test,
                       time.disp.test = time.disp.test,
                       time.tauto.test = time.tauto.test,
                       time.sauto.test = time.sauto.test,
                       fact.unif.test = fact.unif.test,
                       fact.disp.test = fact.disp.test,
                       fact.tauto.test = fact.tauto.test,
                       fact.sauto.test = fact.sauto.test)
  }

return(return.list)

}