novel.by.env.local <- function(comb.df, 
                               novel.list,
                               local.env.data,
                               local.lag,
                               global.lag){

comb.df$bins <- as.numeric(as.character(comb.df$bins))

# get bin position within each time-series
comb.df <- do.call("rbind", lapply(split(comb.df, f=comb.df$site), function(x){
  bin.raw <- as.numeric(as.character(x$bins))
  bin.sort <- sort(as.numeric(as.character(x$bins)))
  return(cbind(x, bin.n = nrow(x)+1 - match(bin.raw, bin.sort)))
}))

# add in site-level data
comb.df <- merge(comb.df, site.df[,c("site","long","lat", "REGION", "elev")],
                 by.x="site", by.y="site",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)

comb.df$modelGlobal <- comb.df[,paste0("globalLag", global.lag)]
comb.df$modelLocal <- comb.df[,paste0("localLag", local.lag)]


print(paste0("Local global temp diff correlation = ", 
             round(cor(cbind(comb.df$modelGlobal, comb.df$modelLocal), 
                 use="complete.obs")[2], 3)))

# drop Antarctica sites
comb.df <- droplevels(comb.df[comb.df$REGION!= "Antartica",])

# drop sites < -10 msl
comb.df <- comb.df[comb.df$elev > -10, ]

tsRichness <- sapply(novel.list$prop.ssmats, ncol)
tsRichness <- data.frame(site=names(tsRichness),
                         tsRichness = tsRichness,
                         tsLength = sapply(novel.list$prop.ssmats, nrow))

comb.df <- merge(comb.df, tsRichness,
                 by.x="site", by.y="site",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)

comb.df.full <- comb.df

comb.df.fit <- comb.df[complete.cases(comb.df[,c("modelGlobal",
                                                 "modelLocal",
                                                 "bin.lag",
                                                 "bin.n",
                                                 "site")]), ]

comb.df.fit$local.lag.scale <- as.vector(scale(comb.df.fit$modelLocal))
comb.df.fit$global.lag.scale <- as.vector(scale(comb.df.fit$modelGlobal))
comb.df.fit$bin.lag.scale = as.vector(scale(log(comb.df.fit$bin.lag)))
comb.df.fit$bin.n.scale = as.vector(scale(log(comb.df.fit$bin.n)))
comb.df.fit$tsLength.scale = as.vector(scale(log(comb.df.fit$tsLength)))
comb.df.fit$tsRichness.scale = as.vector(scale(log(comb.df.fit$tsRichness)))

# ln-transform elevation, adding +11 so -10 msl scores are positive
comb.df.fit$elev.scale <- as.vector(scale(log(comb.df.fit$elev + 11)))

comb.df.fit$bin.num <- as.numeric(as.character(comb.df.fit$bins))

# correlation model ####
env.m <- glm(novel ~ global.lag.scale * local.lag.scale + bin.lag.scale + bin.n.scale + tsLength.scale + 
               tsRichness.scale + elev.scale, 
             family=binomial, data=comb.df.fit)
 return(env.m)
}