subsample.novel.latitude <- function(novel.geo.df,
                                     lat.fact.model,
                                     sample,
                                     iter,
                                     ylims=c(0,0.1)){
  
  print("Subsampling data")
  novel.geo.df$lat.reg <- cut(novel.geo.df$abs.lat, breaks=seq(20,70,10))
  novel.modern <- novel.geo.df[novel.geo.df$bins == "0" &
                               novel.geo.df$abs.lat <=70,]
  
  min.sample = sample
  print(paste0("Minimum bin sample size is ", min.sample))
  
  rownames(novel.modern) <- 1:nrow(novel.modern)
  
  reg.split <- split(rownames(novel.modern), novel.modern$lat.reg)
  
  samp.list <- lapply(1:iter, function(n){
    unlist(lapply(reg.split, function(x){
      sample(x, size = min.sample, replace=TRUE)
    }))
    
  })
  
  # get x-positions of prediction curve
  pred.df <- data.frame(abs.lat = seq(min(novel.modern$abs.lat),
                                      max(novel.modern$abs.lat),
                                      len=200))
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
    sub.df <- novel.modern[as.numeric(samp.list[[n]]), 
                           c("site", "novel", "abs.lat", "bin.lag", "bin.n")]
    sub.df$bin.lag.scale <- as.vector(scale(log(sub.df$bin.lag)))
    sub.df$bin.n.scale <- as.vector(scale(log(sub.df$bin.n)))
    
    temp <- tryCatch(gam(novel ~ s(abs.lat, k=5) + bin.lag.scale + bin.n.scale, 
                         family=binomial, data=sub.df), 
                     warning = function(w) w)
    
    if("warning" %in% class(temp)){
      print("Non-convergence")
      return(NULL)
}
    
    temp.preds <- cbind(pred.df,
                        as.data.frame(predict(temp, 
                                              newdata=pred.df, 
                                              se.fit=TRUE)))
    
    temp.coefs <- list(summary(temp)$p.table,
                        summary(temp)$s.table)
    
    return(list(preds = temp.preds,
           model.class = class(temp),
           edf = temp.coefs[[2]][1,]))
  })
  
  sub.edf <- do.call("rbind", lapply(sub.list, function(x){x$edf}))
  table(sub.edf[,1] <= 1.05)
  
  pred.df.full <- do.call("rbind", lapply(sub.list, function(x){x$preds}))
  
  #pred.df.full <- pred.df.full[pred.df.full$se.fit < 2000,]

  pred.stats <- data.frame(abs.lat = as.numeric(levels(as.factor(pred.df$abs.lat))),
                           mean = plogis(tapply(pred.df.full$fit,
                                         pred.df.full$abs.lat,
                                         mean)),
                           upper = plogis(tapply(pred.df.full$fit,
                                                 pred.df.full$abs.lat,
                                                 quantile, prob=0.975)),
                           lower = plogis(tapply(pred.df.full$fit,
                                                 pred.df.full$abs.lat,
                                                 quantile, prob=0.025)))

  pred.stats <- data.frame(abs.lat = as.numeric(levels(as.factor(pred.df$abs.lat))),
                           mean = tapply(plogis(pred.df.full$fit),
                                                pred.df.full$abs.lat,
                                                mean),
                           upper = tapply(plogis(pred.df.full$fit),
                                                 pred.df.full$abs.lat,
                                                 quantile, prob=0.975),
                           lower = tapply(plogis(pred.df.full$fit),
                                                 pred.df.full$abs.lat,
                                                 quantile, prob=0.025))
  
  
# Plot ####
  sub.dens <- lapply(samp.list, function(sub){
    density(novel.modern$abs.lat[rownames(novel.modern) %in% sub],
            from=min(novel.modern$abs.lat), to=max(novel.modern$abs.lat))
  })
  
  dens.stats <- do.call("cbind", lapply(sub.dens, function(x){x$y}))
  dens.mean <- rowMeans(dens.stats)
  dens.upper <- apply(dens.stats, 1, function(x){quantile(x, 0.975)})
  dens.lower <- apply(dens.stats, 1, function(x){quantile(x, 0.025)})
  
  pdf(date.wrap("./plots/subsample novel comms through latitude", ".pdf"), 
      height=2.5, width=8, useDingbats = FALSE)
  
  split.screen(rbind(c(0.69,0.99,0.15,0.95),
                     c(0.075,0.6,0.15,0.95)))
  
  # density sampling over latitude plot
  screen(1)
  par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=8, tcl=-0.25, 
      mgp=c(3,0.5,0), las=1)
  plot(density(novel.geo.df$abs.lat), col="red", lwd=2,
       xaxt="n", main="")
  
  mtext(side=1, line=0.75, text="Latitude")
  axis(side=1, at=seq(0,90,10), labels=parse(text=paste0(seq(0,90,10), "*degree")),
       mgp=c(3,0,0))
  axis(side=1, at=seq(0,90,5), labels=NA, tcl=-0.125)
  mtext(side=2, line=2, text="Density", las=0)
  
  # apply(dens.stats, 2, function(sub){
  #   lines(sub ~ sub.dens[[1]]$x,
  #         col=rgb(0.5,0.5,0.5,0.075))
  # })
  
  lines(dens.mean ~ sub.dens[[1]]$x, lwd=2)
  lines(dens.upper ~ sub.dens[[1]]$x, lty="31")
  lines(dens.lower ~ sub.dens[[1]]$x, lty="31")
  
  text(x=relative.axis.point(0.015, "x"),
       y=relative.axis.point(0.95, "y"),
       labels="(B)", font=2, adj=0)
  close.screen(1)
  
  screen(2)
  par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=8, tcl=-0.25, 
      mgp=c(3,0.5,0), las=1)
  plot(x=NULL, y=NULL, xlim=c(25,70), 
       ylim=ylims, xlab="", ylab="", xaxt="n", xaxs="i")
  
  mtext(side=2, line=2, las=0, text="Probability of novel community")
  mtext(side=1, line=0.75, text="Latitude")
  axis(side=1, at=seq(0,90,10), labels=parse(text=paste0(seq(0,90,10), "*degree")),
       mgp=c(3,0,0))
  axis(side=1, at=seq(0,90,5), labels=NA, tcl=-0.125)
  
  sub.preds <- do.call("rbind", lapply(1:length(sub.list), function(n){
    temp.preds <- sub.list[[n]]$preds
    
    if(is.null(temp.preds)){return(NULL)}
    lines(plogis(temp.preds$fit) ~ temp.preds$abs.lat,
          col=rgb(0.5,0.5,0.5,0.075), lwd=0.75)
    
    temp.preds$iter = n
    return(temp.preds)
    
  }))
  
  # predictions from full model
  abs.lat.seq = seq(25,70, len=2000)
  fact.pred <- data.frame(abs.lat = rep(seq(25,70, len=length(abs.lat.seq)), 6),
                          bins = rep(seq(0,1000,200), each=length(abs.lat.seq)),
                          bin.lag.scale = 0,
                          bin.n.scale = 0,
                          tsLength.scale = 0,
                          tsRichness.scale = 0,
                          elev.scale = 0)
  
  fact.pred <- cbind(fact.pred,
                     as.data.frame(predict(lat.fact.model, newdata=fact.pred, se.fit=TRUE)))
  fact.pred$raw.fit <- plogis(fact.pred$fit)
  fact.pred$upper <- plogis(fact.pred$fit + 1.96*fact.pred$se.fit)
  fact.pred$lower <- plogis(fact.pred$fit - 1.96*fact.pred$se.fit)
  
  full.preds <- fact.pred[fact.pred$bins == 0 &
                          fact.pred$abs.lat <= 70,]
  full.preds$upper <- plogis(full.preds$fit + 1.96*full.preds$se.fit)
  full.preds$lower <- plogis(full.preds$fit - 1.96*full.preds$se.fit)
  full.preds$raw.fit <- plogis(full.preds$fit)
  
  novel.col <- col2rgb("orange")/255
  polygon(y = c(full.preds$upper,
                rev(full.preds$lower)),
          x = c(full.preds$abs.lat,
                rev(full.preds$abs.lat)),
          border="orange", col=rgb(novel.col[1], novel.col[2], novel.col[3], 0.3))
  lines(plogis(full.preds$fit) ~ full.preds$abs.lat, col="orange", lwd=1.5)
  
  # mean predictions from subset models
  lines(pred.stats$upper ~ pred.stats$abs.lat, lty="31")
  lines(pred.stats$lower ~ pred.stats$abs.lat, lty="31")
  lines(pred.stats$mean ~ pred.stats$abs.lat)
  
  box()
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95, "y"),
       labels="(A)", font=2, adj=0)
  
  close.screen(2)
  
  dev.off()
  
  return(list(model.list=sub.list,
                    mean.preds = pred.stats))
}
