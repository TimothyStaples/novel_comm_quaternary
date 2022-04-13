latitude.model.plot.region <- function(all.nprob.models,
                                lat.model,
                                major.xlims,
                                ylims=c(0,80),
                                name,
                                novelLims=c(0,0.06)){

novel.geo.m <- lat.model$time.model
novel.geo.fact.m <- lat.model$fact.model
comb.df.geo <- lat.model$time.model$model
abs.lat.seq = seq(0,80, len=2000)
bin.seq = seq(1000,max(major.xlims), len=2000)

print("Model predictions")

grid.pred <- expand.grid(abs.lat = abs.lat.seq,
                         bins = bin.seq)
grid.pred$bin.lag.scale = 0
grid.pred$bin.n.scale = 0
grid.pred$ts.scale = 0
grid.pred$tsRichness.scale = 0
grid.pred$elev.scale = 0

grid.pred <- lapply(1:2, function(n){

  temp.pred <- grid.pred
  temp.pred$REGION = c("Europe", "North America")[n]
  temp.pred <- cbind(temp.pred,
         as.data.frame(predict(novel.geo.m, newdata=temp.pred, se.fit=TRUE)))
  temp.pred$raw.fit = plogis(temp.pred$fit)
  temp.pred$upper <- plogis(temp.pred$fit + 1.96*temp.pred$se.fit)
  temp.pred$lower <- plogis(temp.pred$fit - 1.96*temp.pred$se.fit)
  return(temp.pred)
})

# factor spline preds ####
fact.pred <- lapply(1:2, function(n){
fact.pred <- data.frame(abs.lat = rep(seq(0,80, len=length(abs.lat.seq)), 6),
                        bins = rep(seq(0,1000,200), each=length(abs.lat.seq)),
                        bin.lag.scale = 0,
                        bin.n.scale = 0,
                        ts.scale = 0,
                        tsRichness.scale = 0,
                        elev.scale = 0,
                        REGION = c("Europe", "North America")[n])

fact.pred <- cbind(fact.pred,
                   as.data.frame(predict(novel.geo.fact.m, newdata=fact.pred, se.fit=TRUE)))
fact.pred$raw.fit <- plogis(fact.pred$fit)
fact.pred$upper <- plogis(fact.pred$fit + 1.96*fact.pred$se.fit)
fact.pred$lower <- plogis(fact.pred$fit - 1.96*fact.pred$se.fit)
return(fact.pred)
})

print("Generating plot")

# Plot ####
color.grad <- colorRampPalette(c("white", "#b3cedc", "#8a87c7", 
                                 "#80388e", "#4d193c", "black"))(200)

pdf(date.wrap(paste0("./plots/latitude region tensor smooth plot (", name, ")"),".pdf"),
    height=9, width=8, useDingbats = FALSE)
  
  split.screen(rbind(c(0.075,0.55,0.5,0.9), # Euro Time
                     c(0.55,0.85,0.5,0.9), # Euro Fact
                     c(0.075, 0.55, 0.1,0.5), # Nth Am Time
                     c(0.55, 0.85, 0.1,0.5), # Nth Am Fact
                     c(0.90,0.93,0.3125,0.6875))) # Euro Legend))
  
  # Europe 
  screen(1)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  
  plot(x=NULL, y=NULL, xlim = rev(major.xlims), ylim=ylims, xaxs="i", yaxs="i",
       xlab="", ylab="", axes=FALSE)
  
  axis(side=1, at = par("usr")[2] - seq(5000,25000,5000) + par("usr")[1], mgp=c(3,0,0), 
       labels=NA)
  axis(side=1, at=seq(0,25000,1000), tcl=-0.125, labels=NA)
  
  axis(side=2, at=seq(0,90,10), labels=parse(text=paste0(seq(0,90,10), "*degree")))
  axis(side=2, at=seq(0,90,5), labels=NA, tcl=-0.125)
  mtext(side=2, line=1.75, las=0, text="Latitude")
  
  line.col <- "grey50"
  
  axis(side=3, at=par("usr")[2]-19000, tcl=-1.5, labels=NA, col=line.col)
  mtext(side=3, at=par("usr")[2]-19350, line=0.2, text="Start of\nglacial retreat", 
        adj=1, col=line.col, cex=0.9)

  axis(side=3, line=0.35, tcl=0.125, 
       at=par("usr")[2] - c(14700,12700), labels=NA, col=line.col)
  mtext(side=3, at = par("usr")[2]-mean(c(14700,12700)), line=0.4,
       text="Bølling-\nAllerød\ninterstadial",
       adj=c(0.5,0.5), col=line.col, cex=0.9)

  axis(side=3, line=0.35, tcl=0.125, at=par("usr")[2]-c(11000,5000), labels=NA, col=line.col)
  mtext(side=3, at =par("usr")[2]- mean(c(11000,5000)), line=0.4,
        text="Holocene\nThermal\nMaximum",
        adj=c(0.5,0.5), col=line.col, cex=0.9)

  novel.mat <- matrix(grid.pred[[1]]$raw.fit, 
                      nrow=length(bin.seq), 
                      ncol=length(abs.lat.seq), byrow=TRUE)
  novel.mat <- novel.mat[nrow(novel.mat):1,]
  
  # remove predictions beyond scope of observed data
  tempdata <- lat.model$data
  nov.range <- range(tempdata$abs.lat[tempdata$REGION=="Europe"])
  
  novel.mat[, abs.lat.seq < nov.range[1]] = NA
  novel.mat[, abs.lat.seq > nov.range[2]] = NA
  
  image(x=bin.seq,
        y=abs.lat.seq,
        z=novel.mat,
        col=color.grad, add=TRUE,
        useRaster=TRUE,
        zlim=novelLims)

    text(x=relative.axis.point(0.02, "x"),
         y=relative.axis.point(0.97, "y"),
         labels = "(A) Europe", font=2, adj=0)
    
    contour(x=bin.seq,
            y=abs.lat.seq,
            z=novel.mat,
            levels=seq(0.01,round(0.5*max(novelLims),2),0.01), drawlabels=TRUE,
            labels=seq(0.01,round(0.5*max(novelLims), 2),0.01), lwd=0.5,
            add=TRUE, col="black", method="edge")
    
    contour(x=bin.seq,
            y=abs.lat.seq,
            z=novel.mat,
            levels=seq(round(0.5*max(novelLims)+0.01, 2), round(max(novelLims), 2),0.01), drawlabels=TRUE,
            labels=seq(round(0.5*max(novelLims)+0.01, 2), round(max(novelLims), 2),0.01), lwd=0.5,
            add=TRUE, col="white", method="flattest")
    
  box()
  close.screen(1)
  
  screen(3)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  
  plot(x=NULL, y=NULL, xlim = rev(major.xlims), ylim=ylims, xaxs="i", yaxs="i",
       xlab="", ylab="", axes=FALSE)
  
  axis(side=1, at = par("usr")[2] - seq(5000,25000,5000) + par("usr")[1], mgp=c(3,0,0), 
       labels=format(seq(5000,25000,5000), big.mark=","))
  axis(side=1, at=0, labels=0, mgp=c(3,0,0))
  axis(side=1, at=seq(0,25000,1000), tcl=-0.125, labels=NA)
  mtext(side=1, line=1, text="Years before present")
  
  axis(side=2, at=seq(30,90,10), labels=parse(text=paste0(seq(30,90,10), "*degree")))
  axis(side=2, at=seq(30,90,5), labels=NA, tcl=-0.125)
  mtext(side=2, line=1.75, las=0, text="Latitude")
  
  novel.mat <- matrix(grid.pred[[2]]$raw.fit, 
                      nrow=length(bin.seq), 
                      ncol=length(abs.lat.seq), byrow=TRUE)
  novel.mat <- novel.mat[nrow(novel.mat):1,]
  
  # remove predictions beyond scope of observed data
  tempdata <- lat.model$data
  nov.range <- range(tempdata$abs.lat[tempdata$REGION=="North America"])
  
  novel.mat[, abs.lat.seq < nov.range[1]] = NA
  novel.mat[, abs.lat.seq > nov.range[2]] = NA
  
  image(x=bin.seq,
        y=abs.lat.seq,
        z=novel.mat,
        col=color.grad, add=TRUE,
        useRaster=TRUE,
        zlim=novelLims)
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.97, "y"),
       labels = "(C) North America", font=2, adj=0)
  
  contour(x=bin.seq,
          y=abs.lat.seq,
          z=novel.mat,
          levels=seq(0.01,round(0.5*max(novelLims),2),0.01), drawlabels=TRUE,
          labels=seq(0.01,round(0.5*max(novelLims),2),0.01), lwd=0.5,
          add=TRUE, col="black", method="edge")
  
  contour(x=bin.seq,
          y=abs.lat.seq,
          z=novel.mat,
          levels=seq(round(0.5*max(novelLims)+0.01,2),round(max(novelLims),2),0.01), drawlabels=TRUE,
          labels=seq(round(0.5*max(novelLims)+0.01,2),round(max(novelLims),2),0.01), lwd=0.5,
          add=TRUE, col="white", method="flattest")
  
  box()
  close.screen(3)
  
  screen(2)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot(x=NULL, y=NULL, xlim=rev(c(1100,-100)), ylim=ylims, xaxs="i", yaxs="i",
       xlab="", ylab="", axes=FALSE)
  
  axis(side=4, at=seq(0,90,10), labels=parse(text=paste0(seq(0,90,10), "*degree")))
  axis(side=4, at=seq(0,90,5), labels=NA, tcl=-0.125)
  
  axis(side=1, mgp=c(3,0,0), at=0, labels=NA)
  axis(side=1, mgp=c(3,0,0), at=seq(200,1000,200), labels=NA)
  
  tempdata <- lat.model$data
  nov.range <- range(tempdata$abs.lat[tempdata$REGION=="Europe"])
  temp <- do.call("rbind", lapply(1:6, function(x){
    
    bin.center <- unique(fact.pred[[1]]$bins)[x]
    rawfit <- fact.pred[[1]][fact.pred[[1]]$bins == bin.center,]$raw.fit
    latrange <- range(tempdata$abs.lat[tempdata$bins == bin.center &
                                         tempdata$REGION == "Europe"])
    rawfit[abs.lat.seq < latrange[1]] = NA
    rawfit[abs.lat.seq > latrange[2]] = NA
    return(rawfit)
    
  }))
  
  rev.mat <- temp[nrow(temp):1,]
  
  image(x=c(unique(fact.pred[[1]]$bins)-100, 1100),
        y=abs.lat.seq,
        z=matrix(rev.mat, nrow=6, ncol=length(abs.lat.seq)),
        col=color.grad, add=TRUE,
        useRaster=TRUE,
        zlim=novelLims)
  
    sapply(1:6, function(x){
      bin.center <- unique(fact.pred[[1]]$bins)[x] + c(-100, 100)
      rect(xleft=bin.center[1], xright=bin.center[2],
           ybottom=par("usr")[3], ytop=par("usr")[4], lwd=0.5)
      
    })
    
    nov.above <- do.call("cbind", lapply(1:6, function(x){
      
      bin.center <- unique(fact.pred[[1]]$bins)[x]
      plot.bin.pos <- unique(fact.pred[[1]]$bins)[(6:1)[x]]
      
      sub.pred <- fact.pred[[1]][fact.pred[[1]]$bins == bin.center,]
      sub.pred$lower <- plogis(sub.pred$fit - 1.96*sub.pred$se.fit)
      
      contour(x=as.numeric(plot.bin.pos) + c(-100, 100),
              y=sub.pred$abs.lat,
              z=matrix(rbind(temp[x,], temp[x,]), 
                       nrow=2, ncol=length(abs.lat.seq)),
              levels=seq(0.01,round(0.5*max(novelLims),2),0.01), drawlabels=TRUE,
              labels=seq(0.01,round(0.5*max(novelLims),2),0.01), lwd=0.5,
              add=TRUE, col="black", method="edge")
      
      contour(x=as.numeric(plot.bin.pos) + c(-100, 100),
              y=sub.pred$abs.lat,
              z=matrix(rbind(rev.mat[x,], rev.mat[x,]), 
                       nrow=2, ncol=length(abs.lat.seq)),
              levels=seq(round(0.5*max(novelLims)+0.01,2),round(max(novelLims),2),0.01), drawlabels=TRUE,
              labels=seq(round(0.5*max(novelLims)+0.01,2),round(max(novelLims),2),0.01), lwd=0.5,
              add=TRUE, col="white", method="edge")
      
    }))
    
    text(x=relative.axis.point(0.02, "x"),
         y=relative.axis.point(0.97, "y"),
         labels = "(B)", font=2, adj=0)

  box()
  
  close.screen(2)
  
  screen(4)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot(x=NULL, y=NULL, xlim=rev(c(1100,-100)), ylim=ylims, xaxs="i", yaxs="i",
       xlab="", ylab="", axes=FALSE)
  
  axis(side=4, at=seq(0,90,10), labels=parse(text=paste0(seq(0,90,10), "*degree")))
  axis(side=4, at=seq(0,90,5), labels=NA, tcl=-0.125)
  
  axis(side=1, mgp=c(3,0,0), at=0, labels=format(1000, big.mark=","))
  axis(side=1, mgp=c(3,0,0), at=seq(200,1000,200), labels=rev(seq(0,800,200)))
  
  axis(side=1, mgp=c(3,0.6,0),
       at=seq(0,1000,200)[c(1,3,5)],
       labels=paste0("(",rev(c(1950,1750,1550,1350,1150,950))," AD)")[c(1,3,5)], cex.axis=0.75)
  
  axis(side=1, mgp=c(3,0.6,0),
       at=seq(0,1000,200)[c(2,4,6)],
       labels=paste0("(",rev(c(1950,1750,1550,1350,1150,950))," AD)")[c(2,4,6)], cex.axis=0.75)
  
  mtext(side=1, line=1.5, las=0,
        text="Years before present")
  mtext(side=1, line=2, las=0,
        text="(center of 200 year bin)", cex=0.75)
  
  tempdata <- lat.model$data
  nov.range <- range(tempdata$abs.lat[tempdata$REGION=="North America"])
  temp <- do.call("rbind", lapply(1:6, function(x){
    
    bin.center <- unique(fact.pred[[2]]$bins)[x]
    rawfit <- fact.pred[[2]][fact.pred[[2]]$bins == bin.center,]$raw.fit
    latrange <- range(tempdata$abs.lat[tempdata$bins == bin.center &
                                       tempdata$REGION == "North America"])
    print(latrange)
    rawfit[abs.lat.seq < latrange[1]] = NA
    rawfit[abs.lat.seq > latrange[2]] = NA
    return(rawfit)
    
  }))
  
  rev.mat <- temp[nrow(temp):1,]
  
  image(x=c(unique(fact.pred[[2]]$bins)-100, 1100),
        y=abs.lat.seq,
        z=matrix(rev.mat, nrow=6, ncol=length(abs.lat.seq)),
        col=color.grad, add=TRUE,
        useRaster=TRUE,
        zlim=novelLims)
  
  sapply(1:6, function(x){
    bin.center <- unique(fact.pred[[2]]$bins)[x] + c(-100, 100)
    rect(xleft=bin.center[1], xright=bin.center[2],
         ybottom=par("usr")[3], ytop=par("usr")[4], lwd=0.5)
    
  })
  
  nov.above <- do.call("cbind", lapply(1:6, function(x){
    
    bin.center <- unique(fact.pred[[2]]$bins)[x]
    plot.bin.pos <- unique(fact.pred[[2]]$bins)[(6:1)[x]]
    
    sub.pred <- fact.pred[[2]][fact.pred[[2]]$bins == bin.center,]
    sub.pred$lower <- plogis(sub.pred$fit - 1.96*sub.pred$se.fit)
    
    contour(x=as.numeric(plot.bin.pos) + c(-100, 100),
            y=sub.pred$abs.lat,
            z=matrix(rbind(temp[x,], temp[x,]), 
                     nrow=2, ncol=length(abs.lat.seq)),
            levels=seq(0.01,round(0.5*max(novelLims),2),0.01), drawlabels=TRUE,
            labels=seq(0.01,round(0.5*max(novelLims),2),0.01), lwd=0.5,
            add=TRUE, col="black", method="edge")
    
    contour(x=as.numeric(plot.bin.pos) + c(-100, 100),
            y=sub.pred$abs.lat,
            z=matrix(rbind(temp[x,], temp[x,]), 
                     nrow=2, ncol=length(abs.lat.seq)),
            levels=seq(round(0.5*max(novelLims)+0.01,2),round(max(novelLims),2),0.01), drawlabels=TRUE,
            labels=seq(round(0.5*max(novelLims)+0.01,2),round(max(novelLims),2),0.01), lwd=0.5,
            add=TRUE, col="white", method="edge")
    
  }))
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.97, "y"),
       labels = "(D)", font=2, adj=0)
  
  box()
  
  close.screen(4)
  
  screen(5)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=novelLims, axes=FALSE, xlab="", ylab="",
       xaxs="i", yaxs="i")
  image(x=0,
        y=seq(min(novelLims), max(novelLims), len=200),
        z=matrix(1:200, nrow=1),
        col=color.grad, add=TRUE, useRaster=TRUE)
  
  axis(side=4, at=seq(0,max(novelLims),0.01))
  axis(side=4, at=c(seq(0,max(novelLims),0.005)), tcl=-0.125, labels=NA)
  mtext(side=4, line=1.75, text="Novel community emergence probability", las=0)
  box()
  
  # par(xpd=NA)
  # rect(xleft=relative.axis.point(0.5, "x"),
  #      xright = relative.axis.point(1.55, "x"),
  #      ybottom = relative.axis.point(-0.25, "y"),
  #      ytop = relative.axis.point(-0.1, "y"),
  #      border="red", lwd=2.5)
  # text(x=relative.axis.point(1, "x"),
  #      y = relative.axis.point(-0.25, "y"),
  #      pos = 1, offset = 0.5, adj=0,
  #      labels = "> background\nprobability")
  # par(xpd=FALSE)
  
  close.screen(5)
  
  close.screen(all.screens=TRUE)
  dev.off()

}