regions.novel.prob.model.weird <- function(novel.list,
                                     site.df,
                                     all.prob.model,
                                     name,
				                             time.k,
			 	                             fact.k,
				                             lat.lims = c(0,70),
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
  
  comb.df$bin.lag.scale = as.vector(scale(log(comb.df$bin.lag)))
  comb.df$bin.n.scale = as.vector(scale(log(comb.df$bin.n)))
  
  comb.df.geo <- merge(x=comb.df, y=site.df[,c("site", "long", "lat", "REGION", "elev")],
                       by.x="site", by.y="site",
                       all.x=TRUE, all.y=FALSE, sort=FALSE)

  comb.df.geo$abs.lat <- abs(comb.df.geo$lat)

  # drop sites < -10 msl
  comb.df.geo <- comb.df.geo[comb.df.geo$elev > -10 & 
                             comb.df.geo$abs.lat <= lat.lims[2] &
                             comb.df.geo$abs.lat >= lat.lims[1], ]
  
  tsRichness <- sapply(novel.list$prop.ssmats, ncol)
  tsRichness <- data.frame(site=names(tsRichness),
                           tsRichness = tsRichness,
                           tsLength = sapply(novel.list$prop.ssmats, nrow))
  
  comb.df.geo <- merge(comb.df.geo, tsRichness,
                   by.x="site", by.y="site",
                   all.x=TRUE, all.y=FALSE, sort=FALSE)
  
  comb.df.geo$bin.lag.scale = as.vector(scale(log(comb.df.geo$bin.lag)))
  comb.df.geo$bin.n.scale = as.vector(scale(log(comb.df.geo$bin.n)))
  comb.df.geo$tsLength.scale = as.vector(scale(log(comb.df.geo$tsLength)))
  comb.df.geo$tsRichness.scale = as.vector(scale(log(comb.df.geo$tsRichness)))
 
  # ln-transform elevation, adding +11 so -10 msl scores are positive
  comb.df.geo$elev.scale <- as.vector(scale(log(comb.df.geo$elev + 11)))

  comb.df.geo <- comb.df.geo[!is.na(comb.df.geo$novel),]

  comb.df.time <- droplevels(comb.df.geo[comb.df.geo$bins >= time.age.limits[1] &
                                           comb.df.geo$bins <= time.age.limits[2],])
  
  # Run latitude model ####
  print("Running full latitudinal model")
  novel.geo.m <- gam(novel ~ t2(bins, abs.lat, k=time.k) + bin.lag.scale + bin.n.scale +
                        tsLength.scale + tsRichness.scale + elev.scale,
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
  comb.geo.sub$bin.n.scale = as.vector(scale(log(comb.geo.sub$bin.n)))
  comb.geo.sub$tsLength.scale = as.vector(scale(log(comb.geo.sub$tsLength)))
  comb.geo.sub$tsRichness.scale = as.vector(scale(log(comb.geo.sub$tsRichness)))
  
  # ln-transform elevation, adding +11 so -10 msl scores are positive
  comb.geo.sub$elev.scale <- as.vector(scale(log(comb.geo.sub$elev + 11)))
  
  # now model over latitude explicitly for each bin
  print("Running last 1000 years latitudinal model")
  novel.geo.fact.m <- gam(novel ~ t2(bins, abs.lat, k=fact.k) + bin.lag.scale + bin.n.scale +
                            tsLength.scale + tsRichness.scale + elev.scale,
                            family=binomial(), data=comb.geo.sub)
  
  write.csv(summary(novel.geo.fact.m)$p.table,
            date.wrap(paste0("./outputs/modern time latitude novel parametric coefs (", name, ")"), ".csv"))
  
  write.csv(summary(novel.geo.fact.m)$s.table,
            date.wrap(paste0("./outputs/modern time latitude novel smooth coefs (", name, ")"), ".csv"))
  
  
# Model diagnostics ####
if(test.model){
  print("Running model diagnostics")
  
  timeDiag <- modelDiagTests(model  = novel.geo.m,
                 time = novel.geo.m$model$bins,
                 data = comb.df.time,
                 spat.iter = sauto.n,
                 spat.sub.size = sauto.size,
                 sautoSubsample = TRUE)
  
  factDiag <- modelDiagTests(model  = novel.geo.fact.m,
                             time = novel.geo.fact.m$model$bins,
                             data=comb.geo.sub,
                             spat.iter=sauto.n,
                             sautoSubsample = FALSE)
  
}

# Create return list ####
return.list = list(data = comb.df.geo,
                   time.model = novel.geo.m,
                   fact.model = novel.geo.fact.m)

if(test.model){
    return.list = list(data = comb.df.geo,
                       time.model = novel.geo.m,
                       fact.model = novel.geo.fact.m,
                       timeDiag = timeDiag,
                       factDiag = factDiag)
  }

return(return.list)

}