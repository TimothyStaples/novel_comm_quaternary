subsample.novel <- function(novel.list,
                            full.preds,
                            novel.models,
                            time.age.limits, 
                            factor.age.limits,
                            min.sample, 
                            iter,
                            plot.limits){
  
  comb.df <- do.call("rbind", novel.list$novel)
  comb.df$bins <- as.numeric(as.character(comb.df$bins))
  
  # add in site-level info
  comb.df <- droplevels(merge(comb.df, site.df[,c("REGION", "site", "elev", "long", "lat")],
                              by.x="site", by.y="site", all.x=TRUE, all.y=FALSE, sort=FALSE))
  
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
  
  comb.df <- comb.df[complete.cases(comb.df$REGION),]

  s.size <- tapply(comb.df$bins, comb.df$bins, length)
  
  # we'll run reduced models from the year 20K BP to present, wihch gives us a min
  # sample of
  min.sample <- min(s.size)
  print(paste0("Minimum bin sample size is ", min.sample))
  
  print("Subsampling data")
  rownames(comb.df) <- 1:nrow(comb.df)
  
  time.bins <- sort(unique(comb.df$bins))
  # sample all rows of data for each iter
  
  samp.list <- do.call("rbind", lapply(time.bins, function(bin){
    
    sub.rows <- rownames(comb.df)[comb.df$bins == bin]
      
     replicate(iter, sample(sub.rows, min.sample))
     
      }))
  
  # get x-positions of prediction curve
  pred.df <- data.frame(bins = sort(c(seq(time.age.limits[1], time.age.limits[2], len=200),
                                 time.bins)))
  pred.df <- data.frame(bins = pred.df[!duplicated(pred.df$bins),])
  pred.df$bin.lag.scale = 0
  pred.df$bin.n.scale = 0
  pred.df$tsLength.scale = 0
  pred.df$tsRichness.scale = 0
  pred.df$elev.scale = 0
  
  # update model with each set of data
  print("Running sub-sampled models")
  sub.list <- lapply(1:iter, function(n){
    
    print(n)
    
    #subset data
    head(comb.df)
    sub.df <- comb.df[as.numeric(samp.list[,n]),]
    
    sub.df$bin.lag.scale = as.vector(scale(log(sub.df$bin.lag)))
    sub.df$bin.n.scale = as.vector(scale(log(sub.df$bin.n)))
    sub.df$tsLength.scale = as.vector(scale(log(sub.df$tsLength)))
    sub.df$tsRichness.scale = as.vector(scale(log(sub.df$tsRichness)))
    sub.df$elev.scale <- as.vector(scale(log(sub.df$elev + 11)))
    
    temp.m <- try(gam(novel ~ s(bins, k=10) + bin.lag.scale + bin.n.scale +
                        tsLength.scale + tsRichness.scale + elev.scale, 
          family=binomial, data=sub.df))
    
    # if(class(temp.m)[1] != "gamm"){
    #   print("convergence error")
    #   return(NULL)}
    # 
    temp.preds <- cbind(pred.df,
                        as.data.frame(predict(temp.m, 
                                              newdata=pred.df, 
                                              se.fit=TRUE)))
    
    temp.coefs <- list(summary(temp.m)$p.table,
                        summary(temp.m)$s.table)
    
    return(preds = temp.preds)
  })
  
  sub.preds <- do.call("rbind", lapply(1:length(sub.list), function(n){
    temp.preds <- sub.list[[n]]
    temp.preds$iter = n
    return(temp.preds)
  }))
  
  pred.stats <- data.frame(bins = as.numeric(levels(as.factor(sub.preds$bins))),
                                mean = tapply(sub.preds$fit,
                                              sub.preds$bins,
                                              function(x){mean(plogis(x))}),
                                upper = tapply(sub.preds$fit,
                                               sub.preds$bins,
                                               function(x){quantile(plogis(x), prob=0.975)}),
                                lower = tapply(sub.preds$fit,
                                               sub.preds$bins,
                                               function(x){quantile(plogis(x), prob=0.025)}))
  
  # now factor model subsettings
  # get x-positions of prediction curve
  fact.pred.df <- data.frame(bin.fact = seq(factor.age.limits[1], factor.age.limits[2], 200))
  fact.pred.df$bin.lag.scale = 0
  fact.pred.df$bin.n.scale = 0
  fact.pred.df$tsLength.scale = 0
  fact.pred.df$tsRichness.scale = 0
  fact.pred.df$elev.scale = 0
  
  
  # update model with each set of data
  print("Running sub-sampled factor models")
  sub.list.fact <- lapply(1:iter, function(n){
    
    print(n)
    
    #subset data
    sub.df <- comb.df.full[as.numeric(samp.list[,n]),]
    sub.df <- sub.df[sub.df$bins >= factor.age.limits[1] &
                     sub.df$bins <= factor.age.limits[2],]
    
    sub.df$bin.lag.scale = as.vector(scale(log(sub.df$bin.lag)))
    sub.df$bin.n.scale = as.vector(scale(log(sub.df$bin.n)))
    sub.df$tsLength.scale = as.vector(scale(log(sub.df$tsLength)))
    sub.df$tsRichness.scale = as.vector(scale(log(sub.df$tsRichness)))
    sub.df$elev.scale <- as.vector(scale(log(sub.df$elev + 11)))
    
    sub.df$bin.fact <- as.factor(sub.df$bins)
    
    temp.m <- try(gam(novel ~ bin.fact + bin.lag.scale + bin.n.scale +
                        tsLength.scale + tsRichness.scale + elev.scale, 
                      family=binomial, data=sub.df))
    
    # if(class(temp.m)[1] != "gamm"){
    #   print("convergence error")
    #   return(NULL)}
    # 
    temp.preds <- cbind(fact.pred.df,
                        as.data.frame(predict(temp.m, 
                                              newdata=fact.pred.df, 
                                              se.fit=TRUE)))
    
    temp.coefs <- list(summary(temp.m)$p.table,
                       summary(temp.m)$s.table)
    
    return(preds = temp.preds)
  })
  
  sub.preds.fact <- do.call("rbind", lapply(1:length(sub.list.fact), function(n){
    temp.preds <- sub.list.fact[[n]]
    temp.preds$iter = n
    return(temp.preds)
  }))
  
  fact.pred.stats <- data.frame(bin.fact = as.numeric(levels(as.factor(sub.preds.fact$bin.fact))),
                           mean = tapply(sub.preds.fact$fit,
                                         sub.preds.fact$bin.fact,
                                         function(x){mean(plogis(x))}),
                           upper = tapply(sub.preds.fact$fit,
                                          sub.preds.fact$bin.fact,
                                          function(x){quantile(plogis(x), prob=0.975)}),
                           lower = tapply(sub.preds.fact$fit,
                                          sub.preds.fact$bin.fact,
                                          function(x){quantile(plogis(x), prob=0.025)}))
                           
  # Plot ####

  pdf(date.wrap("./plots/subsample novel comms through time", ".pdf"), 
      height=2.5, width=7, useDingbats = FALSE)
    
  split.screen(rbind(c(0.1,0.65,0.15,0.95),
                     c(0.65,0.99,0.15,0.95)))
  
  screen(1)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25,  mgp=c(3,0.5,0), las=1)
  
  plot(x=NULL, y=NULL, xlim=rev(time.age.limits), 
       ylim=plot.limits, xlab="", ylab="", xaxt="n", xaxs="i")
  
  mtext(side=2, line=2, las=0, text="Probability of novel community emergence")
  mtext(side=1, line=1, text="Years before present")
  axis(side=1, mgp=c(3,0,0))
  axis(side=1, at=seq(0,40000,1000), tcl=-0.125, labels=NA)
  
  test <- do.call("rbind", lapply(1:length(sub.list), function(n){
  temp.preds <- sub.list[[n]]
    
  lines(plogis(temp.preds$fit) ~ temp.preds$bins,
          col=rgb(0.65,0.65,0.65,0.075), lwd=0.75)
    
    temp.preds$iter = n
    return(temp.preds)
    
  }))
  
  # predictions from full model
  novel.col <- col2rgb("orange")/255
  
  full.preds$upper <- plogis(full.preds$fit + 1.96*full.preds$se.fit)
  full.preds$lower <- plogis(full.preds$fit - 1.96*full.preds$se.fit)
  
  polygon(y = c(full.preds$upper,
                rev(full.preds$lower)),
          x = c(full.preds$bins,
                rev(full.preds$bins)),
          border=NA, col=rgb(novel.col[1], novel.col[2], novel.col[3], 0.5))
  lines(plogis(full.preds$fit) ~ full.preds$bins, col="orange", lwd=1.5)
  
  # mean predictions from subset models
  lines(pred.stats$upper ~ pred.stats$bins, lty="31")
  lines(pred.stats$lower ~ pred.stats$bins, lty="31")
  lines(pred.stats$mean ~ pred.stats$bins)
  
  text(x=relative.axis.point(0.01, "x"),
       y=relative.axis.point(0.95, "y"),
       labels="(A)", font=2, adj=0)
  
  
  box()
  close.screen(1)
  
  screen(2)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25,  mgp=c(3,0.5,0), las=1)
  
  point.space = 35
  
  plot(x=NULL, y=NULL, xlim=c(1100,-100), ylim=plot.limits, xlab="", ylab="", axes=FALSE)
  axis(side=2, labels=NA)
  axis(side=1, mgp=c(3,0,0))
  # points(plogis(sub.preds.fact$fit) ~ jitter(as.numeric(as.character(sub.preds.fact$bin.fact))+point.space,
  #                                            amount=25),
  #       pch=16, cex=0.8, col=rgb(0.65, 0.65, 0.65, 0.075))
  # segments(x0=as.numeric(as.character(fact.pred.stats$bin.fact)) + point.space,
  #          x1=as.numeric(as.character(fact.pred.stats$bin.fact)) + point.space,
  #          y0=fact.pred.stats$lower,
  #          y1=fact.pred.stats$upper)
  boxplot(plogis(sub.preds.fact$fit) ~ sub.preds.fact$bin.fact,
          at= as.numeric(levels(as.factor(sub.preds.fact$bin.fact)))+point.space,
  add=TRUE, boxwex=40, cex=0.5, axes=FALSE, border="grey")
  points(y=fact.pred.stats$mean, x=as.numeric(as.character(fact.pred.stats$bin.fact)) + point.space,
         pch=16)
  
  # all data models
  all.sub <- novel.models$fact.pred
  segments(x0=as.numeric(as.character(all.sub$bin.fact)) - point.space,
           x1=as.numeric(as.character(all.sub$bin.fact)) - point.space,
           y0=all.sub$upper,
           y1=all.sub$lower)
  points(y=all.sub$fit, x=as.numeric(as.character(all.sub$bin.fact)) - point.space,
         pch=21, bg="orange")
  box()
  text(x=relative.axis.point(0.025, "x"),
       y=relative.axis.point(0.95, "y"),
       labels="(B)", font=2, adj=0)
  
  close.screen(2)
  
  close.screen(all.screens=TRUE)
  dev.off()
  
  
  return(data.frame(model.list=sub.list,
                    mean.preds = pred.stats))
}
