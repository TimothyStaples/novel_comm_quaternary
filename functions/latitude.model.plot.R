latitude.model.plot <- function(all.nprob.models,
                                lat.model,
                                major.xlims,
                                lat.lims = c(0,70),
                                zlim=c(0,0.06),
                                name){

novel.geo.m <- lat.model$time.model
novel.geo.fact.m <- lat.model$fact.model
comb.df.geo <- lat.model$time.model$model
comb.df.full <- lat.model$data
abs.lat.seq = seq(lat.lims[1], lat.lims[2], len=2000)
bin.seq = seq(major.xlims[2],major.xlims[1], len=2000)

# tensor smoother preds
grid.pred <- expand.grid(abs.lat = abs.lat.seq,
                         bins = bin.seq)
grid.pred$bin.lag.scale = 0
grid.pred$bin.n.scale = 0
grid.pred$tsLength.scale = 0
grid.pred$tsRichness.scale = 0
grid.pred$elev.scale = 0
grid.pred <- cbind(grid.pred,
                   as.data.frame(predict(novel.geo.m, newdata=grid.pred, se.fit=TRUE)))
grid.pred$raw.fit = plogis(grid.pred$fit)
grid.pred$upper <- plogis(grid.pred$fit + 1.96*grid.pred$se.fit)
grid.pred$lower <- plogis(grid.pred$fit - 1.96*grid.pred$se.fit)

# predict from model
# Model predictions ####
time.pred.df <- data.frame(bins = seq(0, 25000, len=200),
                           bin.lag.scale = 0,
                           bin.n.scale = 0,
                           bin.scale = seq(min(comb.df.geo$bins),
                                           max(comb.df.geo$bins),
                                           len=200),
                           tsLength.scale = 0,
                           tsRichness.scale = 0,
                           elev.scale = 0)

time.pred.list <- lapply(seq(5,75, 10), function(x){
  temp.df <- time.pred.df
  temp.df$abs.lat = x
  temp.df <- cbind(temp.df,
                   as.data.frame(predict(novel.geo.m, newdata=temp.df, se.fit=TRUE)))
  temp.df$raw.fit = plogis(temp.df$fit)
  return(temp.df)
})

fact.pred <- data.frame(abs.lat = rep(seq(lat.lims[1],lat.lims[2], len=length(abs.lat.seq)), 6),
                        bins = rep(seq(0,1000,200), each=length(abs.lat.seq)),
                        bin.lag.scale = 0,
                        bin.n.scale=0,
                        tsLength.scale = 0,
                        tsRichness.scale = 0,
                        elev.scale = 0)

fact.pred <- cbind(fact.pred,
                   as.data.frame(predict(novel.geo.fact.m, newdata=fact.pred, se.fit=TRUE)))
fact.pred$raw.fit <- plogis(fact.pred$fit)
grid.pred$upper <- plogis(grid.pred$fit + 1.96*grid.pred$se.fit)
grid.pred$lower <- plogis(grid.pred$fit - 1.96*grid.pred$se.fit)

# background threshold from all model
b.novel <- all.nprob.models$time.pred[all.nprob.models$time.pred$bins >=20000 & 
                                        all.nprob.models$time.pred$bins <=25000,]
b.novel$upper <- plogis(b.novel$fit + 1.96*b.novel$se.fit)
b.novel$lower <- plogis(b.novel$fit - 1.96*b.novel$se.fit)

back.low <- min(b.novel$lower)
back.up <- max(b.novel$upper)
  
print("Generating plot")

# tensor product surface plot ####
color.grad <- colorRampPalette(c("#b3cedc", "#8a87c7", 
                                 "#80388e", "#4d193c", "black"))(200)

# factor spline preds ####
fact.pred <- data.frame(abs.lat = rep(seq(lat.lims[1], lat.lims[2], len=length(abs.lat.seq)), 6),
                        bins = rep(seq(0,1000,200), each=length(abs.lat.seq)),
                        bin.lag.scale = 0,
                        bin.n.scale=0,
                        tsLength.scale = 0,
                        tsRichness.scale = 0,
                        elev.scale = 0)

fact.pred <- cbind(fact.pred,
                   as.data.frame(predict(novel.geo.fact.m, newdata=fact.pred, se.fit=TRUE)))
fact.pred$raw.fit <- plogis(fact.pred$fit)
fact.pred$upper <- plogis(fact.pred$fit + 1.96*fact.pred$se.fit)
fact.pred$lower <- plogis(fact.pred$fit - 1.96*fact.pred$se.fit)

pdf(date.wrap(paste0("./plots/latitude region tensor smooth plot (", name, ")"),".pdf"),
      height=4.75, width=8, useDingbats = FALSE)
  
  split.screen(rbind(c(0.075,0.55,0.125,0.9),
                     c(0.55,0.85,0.125,0.9),
                     c(0.90,0.93,0.3125,0.6875)))
  
  screen(1)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  
  plot(x=NULL, y=NULL, xlim = rev(major.xlims), ylim=lat.lims, xaxs="i", yaxs="i",
       xlab="", ylab="", axes=FALSE)
  
  axis(side=1, at = par("usr")[2] - seq(5000,25000,5000) + par("usr")[1], mgp=c(3,0.1,0), 
       labels=format(seq(5000,25000,5000), big.mark=","))
  axis(side=1, at=0, labels=0, mgp=c(3,0.1,0))
  axis(side=1, at=par("usr")[2] - seq(1000,25000,1000) + par("usr")[1], tcl=-0.125, labels=NA)
  mtext(side=1, line=1, text="Years before present")
  
  axis(side=2, at=seq(0,90,10), labels=parse(text=paste0(seq(0,90,10), "*degree")))
  axis(side=2, at=seq(0,90,5), labels=NA, tcl=-0.125)
  mtext(side=2, line=1.75, las=0, text="Latitude")
  
  line.col <- "grey50"
  
  axis(side=3, at=par("usr")[2] + par("usr")[1] - 19000, tcl=-1.5, labels=NA, col=line.col)
  mtext(side=3, at=par("usr")[2] + par("usr")[1] -18800, line=0.2, text="Start of\nglacial\nretreat", 
        adj=0, col=line.col, cex=0.9)

  axis(side=3, line=0.35, tcl=0.125, 
       at=par("usr")[2] + par("usr")[1] - c(14700,12700), labels=NA, col=line.col)
  mtext(side=3, at = par("usr")[2] + par("usr")[1] -mean(c(14700,12700)), line=0.4,
       text="Bølling-\nAllerød\ninterstadial",
       adj=c(0.5,0.5), col=line.col, cex=0.9)

  axis(side=3, line=0.35, tcl=0.125, at=par("usr")[2]  + par("usr")[1] - c(11000,5000), labels=NA, col=line.col)
  mtext(side=3, at =par("usr")[2] + par("usr")[1] - mean(c(11000,5000)), line=0.4,
        text="Holocene\nThermal\nMaximum",
        adj=c(0.5,0.5), col=line.col, cex=0.9)

  novel.mat <- matrix(grid.pred$raw.fit, nrow=length(bin.seq), ncol=length(abs.lat.seq), byrow=TRUE)
  novel.mat <- novel.mat[nrow(novel.mat):1,]
  
  image(x=bin.seq,
        y=abs.lat.seq,
        z=novel.mat,
        col=color.grad, add=TRUE,
        useRaster=TRUE,
        zlim=zlim)

    mtext(side=3,
          at=relative.axis.point(0, "x"),
          line=0.1,
         text = "(A)", font=2, adj=0)
    
    contour(x=bin.seq,
            y=abs.lat.seq,
            z=novel.mat,
            levels=seq(0.01,0.03,0.01), drawlabels=TRUE,
            labels=seq(0.01,0.03,0.01), lwd=0.5,
            add=TRUE, col="black", method="edge")
    
    contour(x=bin.seq,
            y=abs.lat.seq,
            z=novel.mat,
            levels=seq(0.04,max(zlim),0.01), drawlabels=TRUE,
            labels=seq(0.04,max(zlim),0.01), lwd=0.5,
            add=TRUE, col="white", method="flattest")
    
  box()
  close.screen(1)
  
  screen(2)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot(x=NULL, y=NULL, xlim=rev(c(1100,-100)), ylim=lat.lims, xaxs="i", yaxs="i",
       xlab="", ylab="", axes=FALSE)
  
  axis(side=4, at=seq(0,90,10), labels=parse(text=paste0(seq(0,90,10), "*degree")))
  axis(side=4, at=seq(0,90,5), labels=NA, tcl=-0.125)
  
  axis(side=1, mgp=c(3,0.1,0), at=0, labels=format(1000, big.mark=","))
  axis(side=1, mgp=c(3,0.1,0), at=seq(200,1000,200), labels=rev(seq(0,800,200)))
  
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
  
  temp <- do.call("rbind", lapply(1:6, function(x){
    
    bin.center <- unique(fact.pred$bins)[x]
    fact.pred[fact.pred$bins == bin.center,]$raw.fit
    
  }))
  
  rev.mat <- temp[nrow(temp):1,]
  
  image(x=c(unique(fact.pred$bins)-100, 1100),
        y=abs.lat.seq,
        z=matrix(rev.mat, nrow=6, ncol=length(abs.lat.seq)),
        col=color.grad, add=TRUE,
        useRaster=TRUE,
        zlim=zlim)
    
    nov.above <- do.call("cbind", lapply(1:6, function(x){
      
      bin.center <- unique(fact.pred$bins)[x]
      plot.bin.pos <- unique(fact.pred$bins)[(6:1)[x]]
      
      sub.pred <- fact.pred[fact.pred$bins == bin.center,]
      sub.pred$lower <- plogis(sub.pred$fit - 1.96*sub.pred$se.fit)
      
      contour(x=as.numeric(plot.bin.pos) + c(-100, 100),
              y=sub.pred$abs.lat,
              z=matrix(rbind(sub.pred$raw.fit, sub.pred$raw.fit), 
                       nrow=2, ncol=length(abs.lat.seq)),
              levels=seq(0.01,0.03,0.01), drawlabels=TRUE,
              labels=seq(0.01,0.03,0.01), lwd=0.5,
              add=TRUE, col="black", method="edge")
      
      contour(x=as.numeric(plot.bin.pos) + c(-100, 100),
              y=sub.pred$abs.lat,
              z=matrix(rbind(sub.pred$raw.fit, sub.pred$raw.fit), 
                       nrow=2, ncol=length(abs.lat.seq)),
              levels=seq(0.04,0.09,0.01), drawlabels=TRUE,
              labels=seq(0.04,0.09,0.01), lwd=0.5,
              add=TRUE, col="white", method="edge")
      
      cut.off <- sub.pred[sub.pred$lower > back.up,]
      start <- cut.off[1,]
      end <- cut.off[nrow(cut.off),]
      
      if(is.na(start$fit) | nrow(start)==0){return(NULL)}
      
      return(rbind(start$bins, start$abs.lat, end$abs.lat))
    }))
    
    sapply(1:6, function(x){
      bin.center <- unique(fact.pred$bins)[x] + c(-100, 100)
      rect(xleft=bin.center[1], xright=bin.center[2],
           ybottom=par("usr")[3], ytop=par("usr")[4], lwd=0.5)
      
    })
    
    # polygon(x=c(nov.above[1]-90, nov.above[1]+100, nov.above[1]+100, nov.above[1]-90),
    #         y=c(nov.above[2], nov.above[2], nov.above[3], nov.above[3]),
    #         border="red", lwd=2.5)
  
    mtext(side=3,
          at=relative.axis.point(0, "x"),
          line=0.1,
          text = "(B)", font=2, adj=0)

  box()
  
  close.screen(2)
  
  screen(3)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=zlim, axes=FALSE, xlab="", ylab="",
       xaxs="i", yaxs="i")
  image(x=0,
        y=seq(0, max(zlim), len=200),
        z=matrix(1:200, nrow=1),
        col=color.grad, add=TRUE, useRaster=TRUE)
  
  axis(side=4, at=seq(0,max(zlim),0.01))
  axis(side=4, at=c(seq(0,max(zlim),0.005)), tcl=-0.125, labels=NA)
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
  
  close.screen(3)
  
  close.screen(all.screens=TRUE)
  dev.off()


# surface plot with raw data and CIs ####

pdf(date.wrap(paste0("./plots/latitude region tensor smooth SUPPS plot contour (", 
                     name, ": )"), ".pdf"),
    height=13, width=9, useDingbats = FALSE)

split.screen(rbind(c(0.10,0.55,0.69,0.98),
                   c(0.55,0.85,0.69,0.98),
                   c(0.10,0.55,0.38,0.66),
                   c(0.55,0.85,0.38,0.66),
                   c(0.10,0.55,0.07,0.35),
                   c(0.55,0.85,0.07,0.35),
                   c(0.89,0.92,0.425,0.625),
                   c(0.89,0.92,0.11,0.31)))

screen(1)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(x=NULL, y=NULL, xlim=rev(major.xlims), ylim=lat.lims, xaxs="i", yaxs="i",
     xlab="", ylab="", axes=FALSE)

axis(side=1, at=par("usr")[2] - seq(5000,35000,5000) + par("usr")[1], mgp=c(3,0,0), 
     labels=format(seq(5000,35000,5000), big.mark=","))
axis(side=1, at=seq(0,35000,1000), tcl=-0.125, labels=NA)

axis(side=2, at=seq(0,90,10), labels=parse(text=paste0(seq(0,90,10), "*degree")))
axis(side=2, at=seq(0,90,5), labels=NA, tcl=-0.125)
mtext(side=2, line=1.75, las=0, text="Absolute latitude")

points(y=comb.df.geo$abs.lat, x=(par("usr")[2] - comb.df.geo$bins + par("usr")[1]), 
       pch=16, cex=0.25, col="grey60")

novel.mat <- matrix(grid.pred$raw.fit, nrow=length(bin.seq), ncol=length(abs.lat.seq), byrow=TRUE)
novel.mat <- novel.mat[nrow(novel.mat):1,]
contour(x=bin.seq,
        y=abs.lat.seq,
        z=novel.mat,
        levels=seq(0.005,0.06,0.005), drawlabels=TRUE,
        labels=seq(0.005,0.06,0.005), lwd=0.5,
        add=TRUE, col="black", method="flattest")

text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(0.965, "y"),
     labels="(A)", adj=0, font=2, cex=1.25)

box()
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(x=NULL, y=NULL, xlim=rev(c(1100,-100)), ylim=lat.lims, xaxs="i", yaxs="i",
     xlab="", ylab="", axes=FALSE)

axis(side=4, at=seq(0,90,10), labels=parse(text=paste0(seq(0,90,10), "*degree")))
axis(side=4, at=seq(0,90,5), labels=NA, tcl=-0.125)

axis(side=1, mgp=c(3,0,0), at=0, labels=format(1000, big.mark=","))
axis(side=1, mgp=c(3,0,0), at=rev(seq(200,1000,200)), labels=seq(0,800,200))

axis(side=1, mgp=c(3,0.6,0),
     at=seq(0,1000,200)[c(1,3,5)],
     labels=paste0("(",rev(c(1950,1750,1550,1350,1150,950))," AD)")[c(1,3,5)], cex.axis=0.75)

axis(side=1, mgp=c(3,0.6,0),
     at=seq(0,1000,200)[c(2,4,6)],
     labels=paste0("(",rev(c(1950,1750,1550,1350,1150,950))," AD)")[c(2,4,6)], cex.axis=0.75)

comb.geo.sub <- comb.df.full[comb.df.full$bins <= 1000,]

points(y=comb.geo.sub$abs.lat, x= jitter(1000 - comb.geo.sub$bins, amount=75),
       pch=16, cex=0.2, col="grey60")

temp <- do.call("rbind", lapply(1:6, function(x){
  
  bin.center <- unique(fact.pred$bins)[x]
  fact.pred[fact.pred$bins == bin.center,]$raw.fit
  
}))

sapply(1:6, function(x){
  
  bin.center <- unique(fact.pred$bins)[x]
  plot.center <- unique(fact.pred$bins)[(6:1)[x]]
  sub.pred <- fact.pred[fact.pred$bins == bin.center,]
  
  contour(x=as.numeric(plot.center) + c(-100, 100),
          y=sub.pred$abs.lat,
          z=matrix(rbind(sub.pred$raw.fit, sub.pred$raw.fit), nrow=2),
          levels=seq(0.01,0.06,0.01), drawlabels=TRUE,
          labels=seq(0.01,0.06,0.01), lwd=0.5,
          add=TRUE, col="black", method="edge")
  
  contour(x=as.numeric(plot.center) + c(-100, 100),
          y=sub.pred$abs.lat,
          z=matrix(rbind(sub.pred$raw.fit, sub.pred$raw.fit), nrow=2),
          levels=seq(0.05,0.09,0.01), drawlabels=TRUE,
          labels=seq(0.05,0.09,0.01), lwd=0.5,
          add=TRUE, col="white", method="edge")
})

sapply(1:6, function(x){
  bin.center <- unique(fact.pred$bins)[x] + c(-100, 100)
  rect(xleft=bin.center[1], xright=bin.center[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4], lwd=0.5)
  
})
box()

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.97, "y"),
     labels = "(B)", font=2, adj=0)
close.screen(2)

screen(3)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(x=NULL, y=NULL, xlim=rev(major.xlims), ylim=lat.lims, xaxs="i", yaxs="i",
     xlab="", ylab="", axes=FALSE)

axis(side=1, at=par('usr')[2] - seq(5000,35000,5000) + par('usr')[1], mgp=c(3,0,0), 
     labels=format(seq(5000,35000,5000), big.mark=","))
axis(side=1, at=seq(0,35000,1000), tcl=-0.125, labels=NA)

axis(side=2, at=seq(0,90,10), labels=parse(text=paste0(seq(0,90,10), "*degree")))
axis(side=2, at=seq(0,90,5), labels=NA, tcl=-0.125)
mtext(side=2, line=1.75, las=0, text="Absolute latitude")

novel.mat <- matrix(grid.pred$lower, nrow=length(bin.seq), ncol=length(abs.lat.seq), byrow=TRUE)
novel.mat <- novel.mat[nrow(novel.mat):1,]

image(x=bin.seq,
      y=abs.lat.seq,
      z=novel.mat,
      col=color.grad, add=TRUE,
      useRaster=TRUE,
      zlim=c(0,0.05))

contour(x=bin.seq,
        y=abs.lat.seq,
        z=novel.mat,
        levels=seq(0.005,0.04,0.005), drawlabels=TRUE,
        labels=seq(0.005,0.04,0.005), lwd=0.5,
        add=TRUE, col="black", method="flattest")

contour(x=bin.seq,
        y=abs.lat.seq,
        z=novel.mat,
        levels=seq(0.03,0.06,0.005), drawlabels=TRUE,
        labels=seq(0.03,0.06,0.005), lwd=0.5,
        add=TRUE, col="white", method="flattest")

text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(0.965, "y"),
     labels="(C)", adj=0, font=2, cex=1.25)

box()
close.screen(3)

screen(4)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=rev(c(1100,-100)), ylim=lat.lims, xaxs="i", yaxs="i",
     xlab="", ylab="", axes=FALSE)

axis(side=4, at=seq(0,90,10), labels=parse(text=paste0(seq(0,90,10), "*degree")))
axis(side=4, at=seq(0,90,5), labels=NA, tcl=-0.125)

axis(side=1, mgp=c(3,0,0), at=0, labels=format(1000, big.mark=","))
axis(side=1, mgp=c(3,0,0), at=rev(seq(200,1000,200)), labels=seq(0,800,200))

axis(side=1, mgp=c(3,0.6,0),
     at=seq(0,1000,200)[c(1,3,5)],
     labels=paste0("(",rev(c(1950,1750,1550,1350,1150,950))," AD)")[c(1,3,5)], cex.axis=0.75)

axis(side=1, mgp=c(3,0.6,0),
     at=seq(0,1000,200)[c(2,4,6)],
     labels=paste0("(",rev(c(1950,1750,1550,1350,1150,950))," AD)")[c(2,4,6)], cex.axis=0.75)

temp <- do.call("rbind", lapply(1:6, function(x){
  
  bin.center <- unique(fact.pred$bins)[x]
  fact.pred[fact.pred$bins == bin.center,]$lower
  
}))

image(x=c(unique(fact.pred$bins)-100, 1100),
      y=abs.lat.seq,
      z=temp[nrow(temp):1,],
      col=color.grad, add=TRUE,
      useRaster=TRUE,
      zlim=c(0,0.05))

sapply(1:6, function(x){
  
  bin.center <- unique(fact.pred$bins)[x]
  plot.center <- unique(fact.pred$bins)[(6:1)[x]]
  sub.pred <- fact.pred[fact.pred$bins == bin.center,]
  
  contour(x=as.numeric(plot.center) + c(-100, 100),
          y=sub.pred$abs.lat,
          z=rbind(sub.pred$lower, sub.pred$lower),
          levels=seq(0.01,0.02,0.01), drawlabels=TRUE,
          labels=seq(0.01,0.02,0.01), lwd=0.5,
          add=TRUE, col="black", method="edge")
  
  contour(x=as.numeric(plot.center) + c(-100, 100),
          y=sub.pred$abs.lat,
          z=rbind(sub.pred$lower, sub.pred$lower),
          levels=seq(0.03,0.09,0.01), drawlabels=TRUE,
          labels=seq(0.03,0.09,0.01), lwd=0.5,
          add=TRUE, col="white", method="edge")
})

sapply(1:6, function(x){
  bin.center <- unique(fact.pred$bins)[x] + c(-100, 100)
  rect(xleft=bin.center[1], xright=bin.center[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4], lwd=0.5)
  
})
box()

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.97, "y"),
     labels = "(D)", font=2, adj=0)

close.screen(4)

screen(5)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)

upper.color.grad <- colorRampPalette(c("#b3cedc", "#8a87c7", "#80388e", "#4d193c", "black"),
                                     bias=2)(200)

plot(x=NULL, y=NULL, xlim=rev(major.xlims), ylim=lat.lims, xaxs="i", yaxs="i",
     xlab="", ylab="", axes=FALSE)

axis(side=1, at=par("usr")[2] - seq(5000,35000,5000) + par("usr")[1], mgp=c(3,0,0), 
     labels=format(seq(5000,35000,5000), big.mark=","))
axis(side=1, at=seq(0,35000,1000), tcl=-0.125, labels=NA)
mtext(side=1, line=1, text="Years before present")

axis(side=2, at=seq(0,90,10), labels=parse(text=paste0(seq(0,90,10), "*degree")))
axis(side=2, at=seq(0,90,5), labels=NA, tcl=-0.125)
mtext(side=2, line=1.75, las=0, text="Absolute latitude")

novel.mat <- matrix(grid.pred$upper, nrow=length(bin.seq), ncol=length(abs.lat.seq), byrow=TRUE)
novel.mat <- novel.mat[nrow(novel.mat):1,]

image(x=bin.seq,
      y=abs.lat.seq,
      z=novel.mat,
      col=upper.color.grad, add=TRUE,
      useRaster=TRUE,
      zlim=c(0,0.2))

contour(x=bin.seq,
        y=abs.lat.seq,
        z=novel.mat,
        levels=seq(0,0.1,0.02), drawlabels=TRUE,
        labels=seq(0,0.1,0.02), lwd=0.5,
        add=TRUE, col="black", method="flattest")

contour(x=bin.seq,
        y=abs.lat.seq,
        z=novel.mat,
        levels=seq(0.1,0.2,0.02), drawlabels=TRUE,
        labels=seq(0.1,0.2,0.02), lwd=0.5,
        add=TRUE, col="black", method="flattest")

text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(0.965, "y"),
     labels="(E)", adj=0, font=2, cex=1.25, col="white")

box()
close.screen(5)

screen(6)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=rev(c(1100,-100)), ylim=lat.lims, xaxs="i", yaxs="i",
     xlab="", ylab="", axes=FALSE)

axis(side=4, at=seq(0,90,10), labels=parse(text=paste0(seq(0,90,10), "*degree")))
axis(side=4, at=seq(0,90,5), labels=NA, tcl=-0.125)

axis(side=1, mgp=c(3,0,0), at=0, labels=format(1000, big.mark=","))
axis(side=1, mgp=c(3,0,0), at=rev(seq(200,1000,200)), labels=seq(0,800,200))

axis(side=1, mgp=c(3,0.6,0),
     at=seq(0,1000,200)[c(1,3,5)],
     labels=paste0("(",rev(c(1950,1750,1550,1350,1150,950))," AD)")[c(1,3,5)], cex.axis=0.75)

axis(side=1, mgp=c(3,0.6,0),
     at=seq(0,1000,200)[c(2,4,6)],
     labels=paste0("(",rev(c(1950,1750,1550,1350,1150,950))," AD)")[c(2,4,6)], cex.axis=0.75)

temp <- do.call("rbind", lapply(1:6, function(x){
  
  bin.center <- unique(fact.pred$bins)[x]
  fact.pred[fact.pred$bins == bin.center,]$upper
  
}))

image(x=c(unique(fact.pred$bins)-100, 1100),
      y=abs.lat.seq,
      z=temp[nrow(temp):1,],
      col=color.grad, add=TRUE,
      useRaster=TRUE,
      zlim=c(0,0.2))

sapply(1:6, function(x){
  
  bin.center <- unique(fact.pred$bins)[x]
  plot.center <- unique(fact.pred$bins)[(6:1)[x]]
  sub.pred <- fact.pred[fact.pred$bins == bin.center,]
  
  contour(x=as.numeric(plot.center) + c(-100, 100),
          y=sub.pred$abs.lat,
          z=rbind(sub.pred$upper, sub.pred$upper),
          levels=seq(0,0.2,0.02), drawlabels=TRUE,
          labels=seq(0,0.2,0.02), lwd=0.5,
          add=TRUE, col="black", method="edge")
  
})

sapply(1:6, function(x){
  bin.center <- unique(fact.pred$bins)[x] + c(-100, 100)
  rect(xleft=bin.center[1], xright=bin.center[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4], lwd=0.5)
  
})
box()

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.97, "y"),
     labels = "(F)", font=2, adj=0)

close.screen(6)

screen(7)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,0.05), axes=FALSE, xlab="", ylab="",
     xaxs="i", yaxs="i")
image(x=0,
      y=seq(0, 0.05, len=200),
      z=matrix(1:200, nrow=1),
      col=color.grad, add=TRUE, useRaster=TRUE)

axis(side=4, at=seq(0,0.07,0.01))
axis(side=4, at=c(seq(0.005,0.075,0.005)), tcl=-0.125, labels=NA)
mtext(side=4, line=1.75, text="Probability of novel community", las=0)
box()
close.screen(7)

screen(8)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,0.2), axes=FALSE, xlab="", ylab="",
     xaxs="i", yaxs="i")
image(x=0,
      y=seq(0, 0.2, len=200),
      z=matrix(1:200, nrow=1),
      col=upper.color.grad, add=TRUE, useRaster=TRUE)

axis(side=4, at=seq(0,0.2,0.02))
axis(side=4, at=c(seq(0,0.2,0.01)), tcl=-0.125, labels=NA)
mtext(side=4, line=1.75, text="Probability of novel community", las=0)
box()
close.screen(8)

close.screen(all.screens=TRUE)
dev.off()


}