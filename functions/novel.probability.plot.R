novel.prob.plot <- function(prob.model.list, 
                            env.data,
                            mod.env.data,
                            name = "all",
                            time.age.limits,
                            factor.age.limits,
                            group.letters,
                            ylims = c(0,0.065),
                            regylims=c(0,0.15)){

comb.df <- prob.model.list$time.data
comb.df.fact <- prob.model.list$fact.data
time.m <- prob.model.list$time.m
time.fact.m <- prob.model.list$fact.model
pred.df <- prob.model.list$time.pred
reg.pred.df <- prob.model.list$time.reg.pred
fact.pred.df <- prob.model.list$fact.pred
fact.pred.df.reg <- prob.model.list$fact.reg.pred

axis.in.plot <- function(side, value=NULL, prop=NULL, at, tcl){
  
  figure.size = par("pin") * 2.54
  margin.size = par("mai") * 2.54
  
  if(side==1){
    plot.size = figure.size[2]
    target.margins = par("usr")[c(3:4)]
  } else {
    plot.size = figure.size[1]
    start.margin = par("usr")[1]
  }
  
  if(value){
    prop = (target.margins[1] - value) / (target.margins[1] - target.margins[2])
  }
  
  axis.line = (prop * plot.size) / (0.2*2.54)
  
  axis(side=side, line=-axis.line, at=at, labels=NA, tcl=tcl)
}

require(shape)
pdf(date.wrap(paste0("./plots/novel comms through time (", name, ") - no raw"), ".pdf"), 
    height=6.5, width=10, useDingbats = FALSE)

close.screen(all.screens=TRUE)
split.screen(rbind(c(0.6,0.91,0.695,0.995), # nov fact
                   c(0.085,0.6,0.695,0.995), # nov time
                   c(0.6,0.91,0.095,0.395), # temp fact
                   c(0.085,0.6,0.095,0.395), # temp time
                   c(0.085,0.6,0.095,0.995), # background
                   c(0.6,0.91,0.395,0.695), # reg nov fact
                   c(0.085,0.6,0.395,0.695))) # reg nov time

screen(5)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=rev(c(1100, time.age.limits[2])), ylim=c(0,1), 
     axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")

line.col <- rgb(0.3,0.3,0.3,0.8)

par(lheight=0.85)
abline(v=19000, lwd=1.5, col=line.col)
text(x=19250, y=relative.axis.point(0.16, "y"),
     labels="Start of\nglacial\nretreat",
     adj=c(1,1), col=line.col, cex=0.9)

rect(xleft=11000, xright=5000, ybottom=par("usr")[3], ytop=par("usr")[4],
     col="grey85", border=NA)
text(x=mean(c(5000,11000)), 
     y=relative.axis.point(0.15, "y"),
     labels="Holocene\nThermal\nMaximum",
     adj=c(0.5,0.5), col=line.col, cex=0.9)

rect(xleft=14700, xright=12700, ybottom=par("usr")[3], ytop=par("usr")[4],
     col="grey85", border=NA)
text(x=mean(c(14700,12700)), 
     y=relative.axis.point(0.23, "y"),
     labels="Bølling-\nAllerød\ninterstadial",
     adj=c(0.5,0.5), col=line.col, cex=0.9)

par(lheight=1)
close.screen(5)

screen(2)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(x=NULL, y=NULL, xlim=rev(c(1100, time.age.limits[2])),
     ylim=ylims,
     xlab="", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i")

ylims <- par("usr")[3:4]

# background novel rates
b.novel <- pred.df[pred.df$bins >=20000 & pred.df$bins <=25000,]
b.novel$upper <- plogis(b.novel$fit + 1.96*b.novel$se.fit)
b.novel$lower <- plogis(b.novel$fit - 1.96*b.novel$se.fit)

back.low <- min(b.novel$lower)
back.up = max(b.novel$upper)
back.mean = mean(plogis(b.novel$fit))

pred.df$upper <- plogis(pred.df$fit + 1.96*pred.df$se.fit)
pred.df$lower <- plogis(pred.df$fit - 1.96*pred.df$se.fit)

axis(side=1, mgp=c(3,0,0), at=seq(5000, 35000, 5000), 
     labels=NA)
axis(side=1, mgp=c(3,0,0), at=0, labels=NA)
axis(side=1, at=seq(0,40000,1000), tcl=-0.125, labels=NA)

axis(side=2)
axis(side=2, at=seq(0,0.1,0.005), labels=NA, tcl=-0.125)

mtext(side=2, line=2.25, las=0,at=par("usr")[3],
      text="Novel community emergence probability")

prop <- tapply(comb.df$novel, comb.df$bins, function(x){sum(x)/length(x)})
sample <- tapply(comb.df$novel, comb.df$bins, function(x){length(x)})
size <- cut(sample, breaks=c(1,10,100,1000,2000))

novel.col <- col2rgb("orange")/255

sub.df <- pred.df[pred.df$bins > 4000,]

post.glac = c(sub.df$bins[rev(which(sub.df$lower > back.up))[1]],
              sub.df$bins[which(sub.df$lower > back.up)[1]])

novel.col <- col2rgb("orange")/255
polygon(x=c(pred.df$bins, rev(pred.df$bins)),
        y=c(plogis(pred.df$fit + 1.96*pred.df$se.fit),
            rev(plogis(pred.df$fit - 1.96*pred.df$se.fit))),
        col=rgb(novel.col[1], novel.col[2], novel.col[3], 0.5),
        border=NA)

lines(plogis(pred.df$fit) ~ pred.df$bins, lwd=3,
      col="orange")

post.glac.max <- which.max(plogis(sub.df$fit))
post.glac.max <- sub.df[post.glac.max,]

segments(x0=rep(post.glac.max$bins, 2),
         x1=rep(par("usr")[2], 2),
         y0=c(post.glac.max$upper, post.glac.max$lower),
         y1=c(post.glac.max$upper, post.glac.max$lower),
         lty="21")

text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.925, "y"),
     labels="(A)", font=2)

box()
close.screen(2)

screen(1)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(x=NULL, y=NULL, 
     xlim=rev(c(-75,factor.age.limits[2]+100)), ylim=ylims, xaxt="n",
     xlab="", ylab="", xaxs="i", yaxs="i", yaxt="n")
  
segments(x0=rep(par("usr")[2], 2),
         x1=rep(par("usr")[1], 2),
         y0=c(post.glac.max$upper, post.glac.max$lower),
         y1=c(post.glac.max$upper, post.glac.max$lower),
         lty="21")

axis(side=4)
axis(side=4, at=seq(0,0.1,0.005), labels=NA, tcl=-0.125)
mtext(side=4, line=2.25, las=0, at=par("usr")[3],
      text="Novel community emergence probability")

axis(side=1, mgp=c(3,0,0), at=1000, labels=NA)
axis(side=1, mgp=c(3,0,0), at=seq(0,800,200), labels=NA)

axis(side=1, mgp=c(3,0.6,0),
     at=seq(0,1000,200),
     labels=NA, cex.axis=0.75)

if(!is.null(time.fact.m)){

  segments(x0=fact.pred.df$bin.num, x1=fact.pred.df$bin.num,
         y0 = fact.pred.df$upper, y1 = fact.pred.df$lower)

text(x=0,
     y=fact.pred.df$lower[fact.pred.df$bin.num==0],
     labels='Modern', pos=1, offset=0.35, adj=1)

if(sum(fact.pred.df$bin.num == 200) > 0){

text(x=200,y=fact.pred.df$lower[fact.pred.df$bin.num==200],
     labels='Pre-modern', pos=1, offset=0.35,adj=0)
}

points(fact.pred.df$fit ~ fact.pred.df$bin.num,
       pch=c(23, rep(21, nrow(fact.pred.df)-1)), bg=c(rep("orange", 2), rep("grey70", 4)),
       cex = c(1.5, rep(1, nrow(fact.pred.df)-1)))

box()

par(xpd=NA)
rect(xleft=relative.axis.point(-0.05, "x"),
     xright=relative.axis.point(0.05, "x"),
     ybottom=post.glac.max$lower + 0.00015,
     ytop=post.glac.max$upper - 0.00015,
     border=NA, col="white")
text(x=par("usr")[1],
     y=mean(unlist(post.glac.max[1,c("upper","lower")])),
     labels="Maximum post-glacial novelty")
par(xpd=FALSE)

}

text(x=relative.axis.point(0.05, "x"),
     y=relative.axis.point(0.925, "y"),
     labels="(B)", font=2)

close.screen(1)

screen(3)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL,  xlim=rev(c(-75,factor.age.limits[2]+100)), 
     ylim=c(-0.45,0.5), axes=FALSE, xlab="", ylab="",
     xaxs="i", yaxs="i")

mod.env.sub <- env.data[env.data$bin <= 1100,]

point.gap <- 20
segments(x0=seq(0,1000,200) + point.gap, 
         x1 = seq(0,1000,200) + point.gap,
         y0 = mod.env.sub$temp + 1.96 * mod.env.sub$temp.se,
         y1 = mod.env.sub$temp - 1.96 * mod.env.sub$temp.se)

points(y = mod.env.sub$temp, x= seq(0,1000,200) + point.gap, pch=21, bg="red")

mod.env.agg <- bin.env.data(mod.env.data[!is.na(mod.env.data$temp),], 
                            200, c(-100,20000), "temp")
mod.env.alt <- mod.env.agg[mod.env.agg$bin <= 1100, ]

segments(x0 = seq(0,1000,200) - point.gap, 
         x1 = seq(0,1000,200) - point.gap,
         y0 = mod.env.alt$env + 1.96 * mod.env.alt$env.se,
         y1 = mod.env.alt$env - 1.96 * mod.env.alt$env.se)

points(y = mod.env.alt$env, 
       x = seq(0,1000,200) - point.gap , pch=21, bg="gold")

axis(side=4, at=seq(-1.2,1,0.2), 
     labels=parse(text=paste0(round(seq(-1.2,1,0.2),1), "*degree")))
axis(side=4, at=seq(-1,1,0.1), tcl=-0.125, labels=NA)
mtext(side=4, line=2.25, las=0,
      text="Temperature anomaly")
mtext(side=4, line=3.25, las=0,
      text=expression("("*Delta*degree*"C from 1961-1990 mean)"))

axis(side=1, mgp=c(3,0.1,0), at=1000, labels=format(1000, big.mark=","))
axis(side=1, mgp=c(3,0.1,0), at=seq(0,800,200))

sapply(1:length(seq(0,1000,200)), function(n){
axis(side=1, mgp=c(3,0.7,0),
     at=seq(0,1000,200)[n],
     labels=paste0("(",c(1950,1750,1550,1350,1150,950)[n]," AD)"), cex.axis=0.75)
})

mtext(side=1, line=1.75, las=0,
      text="Years before present")
mtext(side=1, line=2.35, las=0,
      text="(center of 200 year bin)", cex=0.75)

text(x=relative.axis.point(0.05, "x"),
     y=relative.axis.point(0.925, "y"),
     labels="(F)", font=2)
box()
close.screen(3)

screen(4)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL,  xlim=rev(c(1100, time.age.limits[2])), 
     ylim=c(-3.5,0.9), axes=FALSE, xlab="", ylab="",
     xaxs="i", yaxs="i")
box()

env.data <- env.data[env.data$bin <= 22000,]

polygon(x=c(env.data$bin, rev(env.data$bin)),
        y=c(env.data$temp + 1.96 * env.data$temp.se,
            rev(env.data$temp - 1.96 * env.data$temp.se)),
        border=NA, col=rgb(1,0,0,0.3))

lines(y=env.data$temp, x=env.data$bin, lwd=2, col="red")

axis(side=1, mgp=c(3,0.1,0), at=seq(5000, 35000, 5000), 
     labels=format(seq(5000, 35000, 5000), big.mark=","))
axis(side=1, mgp=c(3,0.1,0), at=0)
mtext(side=1, line=1.5, text="Years before present")
axis(side=1, at=seq(0,40000,1000), tcl=-0.125, labels=NA)

axis(side=2, at=seq(-5,1,1), labels=parse(text=paste0(seq(-5,1,1), "*degree")))
axis(side=2, at=seq(-5,1,1), labels=NA)
axis(side=2, at=seq(-5,1,0.5), tcl=-0.125, labels=NA)
mtext(side=2, line=3, las=0,
      text="Temperature anomaly")
mtext(side=2, line=2, las=0,
      text=expression("("*Delta*degree*"C from 1961-1990 mean)"))

text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.925, "y"),
     labels="(E)", font=2)

box()
close.screen(4)

screen(7)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(x=NULL, y=NULL, xlim=rev(c(1100, time.age.limits[2])),
     ylim=regylims,
     xlab="", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i")

ylims <- par("usr")[3:4]

axis(side=1, mgp=c(3,0,0), at=seq(5000, 35000, 5000), 
     labels=NA)
axis(side=1, mgp=c(3,0,0), at=0, labels=NA)
axis(side=1, at=seq(0,40000,1000), tcl=-0.125, labels=NA)
axis(side=2, at=seq(0,0.2,0.02))
axis(side=2, at=seq(0,0.1,0.005), labels=NA, tcl=-0.125)

sapply(1:length(reg.pred.df), function(n){
  
  temp.pred <- reg.pred.df[[n]]
  
  if(sum(is.na(temp.pred$fit))==nrow(temp.pred)){return(NULL)}
  if(temp.pred$REGION[1] %in% c("North America", "Europe")){
    width = 3
    temp.col = ifelse(temp.pred$REGION == "North America",
                      "#629D26", "blue")
    temp.lty <- ifelse(temp.pred$REGION == "North America",
                                  "solid", "solid")
    temp.rgb <- as.vector(col2rgb(temp.col)/255)
    
    polygon(y=plogis(c(temp.pred$fit + 1.96 * temp.pred$se.fit,
                      rev(temp.pred$fit - 1.96 * temp.pred$se.fit))),
             x=c(temp.pred$bins, rev(temp.pred$bins)),
             border=NA, col=rgb(temp.rgb[1], temp.rgb[2], temp.rgb[3], 0.3))
    
  } else {
    return(NULL)
    width=0.75
    temp.col="grey70"
  }

  lines(plogis(temp.pred$fit) ~ temp.pred$bins, lwd=width,
        col=temp.col, lty=temp.lty)

})

text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.925, "y"),
     labels="(C)", font=2)

rect(xleft=relative.axis.point(0.2, "x"),
     xright=relative.axis.point(0.30, "x"),
     ybottom=relative.axis.point(0.88, "y"),
     ytop=par("usr")[4], border=NA, col="white")

text(x=relative.axis.point(0.215, "x"),
     y=relative.axis.point(0.925, "y"),
     labels=expression(bold("North America"*phantom(" & Europe"))), 
     font=2, col="#629D26")
text(x=relative.axis.point(0.215, "x"),
     y=relative.axis.point(0.925, "y"),
     labels=expression(bold(phantom("North America")*" & "*phantom("Europe"))), 
     font=2, col="black")
text(x=relative.axis.point(0.215, "x"),
     y=relative.axis.point(0.925, "y"),
     labels=expression(bold(phantom("North America & ")*"Europe")), 
     font=2, col="blue")

box()
close.screen(7)

screen(6)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(x=NULL, y=NULL, 
     xlim=rev(c(-75,factor.age.limits[2]+100)), ylim=regylims, xaxt="n",
     xlab="", ylab="", xaxs="i", yaxs="i", yaxt="n")

if(!is.null(prob.model.list$fact.reg.pred)){
  
sub.reg.df <- droplevels(fact.pred.df.reg[fact.pred.df.reg$REGION %in% c("Europe", "North America"),])

regoffset <- seq(-25,25, len=length(unique(sub.reg.df$REGION)))

axis(side=4, at=seq(0,0.06,0.02))
axis(side=4, at=seq(0,0.1,0.005), labels=NA, tcl=-0.125)

axis(side=1, mgp=c(3,0,0), at=1000, labels=NA)
axis(side=1, mgp=c(3,0,0), at=seq(0,800,200), labels=NA)

sapply(1:length(unique(sub.reg.df$REGION)), function(n){

sub.pred <- sub.reg.df[sub.reg.df$REGION == unique(sub.reg.df$REGION)[n],]
if(nrow(sub.pred)==0){return(NULL)}
suboff = regoffset[n]  
temp.col = ifelse(sub.pred$REGION == "North America",
                  "#629D26", "blue")

if(!is.null(time.fact.m)){
  
  with(sub.pred[sub.pred$upper != 1 & sub.pred$lower != 0, ],
  segments(x0=bin.num + suboff, x1=bin.num + suboff,
           y0 = upper, y1 = lower))
  
  points(y=sub.pred$fit, x=sub.pred$bin.num + suboff,
         pch=c(23, rep(21, nrow(sub.pred)-1)), 
         bg= temp.col,
         cex = c(1.5, rep(1, nrow(sub.pred)-1)))
}
})
  
  text(x=0,
       y=fact.pred.df$lower[fact.pred.df$bin.num==0],
       labels='Modern', pos=1, offset=0.35, adj=1)
  
  if(sum(fact.pred.df$bin.num == 200) > 0){
    
    text(x=200,y=fact.pred.df$lower[fact.pred.df$bin.num==200],
         labels='Pre-modern', pos=1, offset=0.35,adj=0)
  }
}
  
  box()
  
  text(x=relative.axis.point(0.05, "x"),
       y=relative.axis.point(0.925, "y"),
       labels="(D)", font=2)
  
close.screen(6)



close.screen(all.screens=TRUE)
dev.off()

# Plot with raw data ####

pdf(date.wrap(paste0("./plots/novel comms through time (", name, ") raw data"), ".pdf"), 
    height=2.75, width=7, useDingbats = FALSE)

par(mar=c(2,3,2,1), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=rev(c(-400, time.age.limits[2]+400)), 
     ylim=c(-0.002,0.09), 
     xlab="", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i")

ylims <- par("usr")[3:4]

axis(side=1, mgp=c(3,0,0), at=seq(5000, 35000, 5000), 
     labels=format(seq(5000, 35000, 5000), big.mark=","))
axis(side=1, mgp=c(3,0,0), at=0)
mtext(side=1, line=1, text="Years before present")
axis(side=1, at=seq(0,40000,1000), tcl=-0.125, labels=NA)

axis(side=2)
axis(side=2, at=seq(0,0.07,0.005), labels=NA, tcl=-0.125)

# axis(side=1, tcl=0, at=1100, labels="...", mgp=c(3,0,0))

mtext(side=2, line=2, las=0,
      text="Novel community emergence probability")

prop <- tapply(comb.df$novel, comb.df$bins, function(x){sum(x, na.rm=TRUE)/length(x[!is.na(x)])})
sample <- tapply(comb.df$novel, comb.df$bins, function(x){length(x[!is.na(x)])})
size <- cut(sample, breaks=c(1,10,100,1000,2000))

novel.col <- col2rgb("orange")/255
points(prop ~ as.numeric(names(prop)), lwd=0.5,
       pch=21, bg=rgb(novel.col[1], novel.col[2], novel.col[3], 0.3),
       cex=c(0.25,0.5,0.75,1)[cut(sample, breaks=c(1,10,50,500,5000))])

polygon(x=c(pred.df$bins, rev(pred.df$bins)),
        y=c(plogis(pred.df$fit + 1.96*pred.df$se.fit),
            rev(plogis(pred.df$fit - 1.96*pred.df$se.fit))),
        col=rgb(novel.col[1], novel.col[2], novel.col[3], 0.5),
        border=NA)

lines(plogis(pred.df$fit) ~ pred.df$bins, lwd=3,
      col="orange")

leg.pos <- relative.axis.point(seq(0.35,0.65, len=4), "x")

par(xpd=NA)
mtext(side=3, line=0.65, text="Sample size", font=2,
      at=mean(leg.pos)-1500)
points(y=rep(relative.axis.point(1.04, "y"),4),
       x=leg.pos, pch=21, cex=c(0.25,0.5,0.75,1), lwd=0.5,
       bg=rgb(novel.col[1], novel.col[2], novel.col[3], 0.3))
text(y=rep(relative.axis.point(1.035, "y"),4), adj=c(0, 0.5),
     x=leg.pos, pch=21, pos=4, offset=0.35,
     labels=c("1-10", "11-50", "51-500", "500+"))
par(xpd=FALSE)

box()
dev.off()


pdf(date.wrap(paste0("./plots/novel comms through time (", name, ") raw region data"), ".pdf"), 
    height=6.5, width=6.5, useDingbats = FALSE)

ys <- rev(seq(0.05,0.975, len=7))

split.screen(rbind(cbind(0.1, 0.65, ys[-1], ys[-length(ys)]),
             cbind(0.65, 0.99, ys[-1], ys[-length(ys)])))

colvect <- c("red", "goldenrod", "brown", "blue", "#629D26", "purple")

sapply(1:length(reg.pred.df), function(n){

  screen(n)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  
  subpred <- reg.pred.df[[n]]
  sub <- comb.df[comb.df$REGION == subpred$REGION[1],]
  
  plot(x=NULL, y=NULL, xlim=c(25000,0), ylim=c(0,0.22),
       axes=FALSE, xlab="", ylab="")
  box()
  
  if(n == 6){axis(side=1, at=seq(0,25000,5000),
                      labels=c(0,format(seq(5000,25000,5000), big.mark=",")), mgp=c(3,0,0))} else {axis(side=1, labels=NA)}
  if(n==6){mtext(side=1,line=1, text="Years before present", at=par("usr")[2])}
  axis(side=2)
  if(n == 3){mtext(side=2, line=2, text="Probability of novel community emergence", las=0)}
  axis(side=2, at=seq(0,0.2,0.01), tcl=-0.125, labels=NA)
  
  subNov <- tapply(sub$novel, sub$bins, function(x){sum(x, na.rm=TRUE)/length(x[!is.na(x)])})
  points(subNov ~ as.numeric(names(subNov)), pch=21, bg=colvect[n], lwd=0.5)
  
  col1 <- as.vector(col2rgb(colvect[n])/255)
  
  polygon(x=c(subpred$bins, rev(subpred$bins)),
          y=plogis(c(subpred$fit + 1.96 * subpred$se.fit,
                     rev(subpred$fit - 1.96 * subpred$se.fit))),
          border=NA, col=rgb(col1[1],col1[2],col1[3],0.3))
  
  lines(plogis(subpred$fit) ~ subpred$bins, col=colvect[n], lwd=1.5)

  text(x=relative.axis.point(0.025, "x"), 
       y=relative.axis.point(0.875, "y"),
       label=paste0("(", LETTERS[seq(1,11,2)[n]], ") ", subpred$REGION[1]), adj=0,
       font=2)
  close.screen(n)
  
})


sapply(1:length(reg.pred.df), function(n){
  
  screen(n+6)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  
  subpred <- split(fact.pred.df.reg, f= fact.pred.df.reg$REGION)[[n]]
  
  plot(x=NULL, y=NULL, xlim=c(1100,0), ylim=c(0,0.22),
       axes=FALSE, xlab="", ylab="")
  box()
  
  if(n == 6){axis(side=1, mgp=c(3,0,0))} else {axis(side=1, labels=NA)}
  axis(side=2, labels=NA)
  axis(side=2, at=seq(0,0.2,0.01), tcl=-0.125, labels=NA)
  
  segments(x0 = as.numeric(as.character(subpred$bin.fact)),
           x1 = as.numeric(as.character(subpred$bin.fact)),
           y0 = subpred$upper, y1 = subpred$lower)
  points(subpred$fit ~ as.numeric(as.character(subpred$bin.fact)), pch=21, bg=colvect[n], lwd=0.5)
  
  text(x=relative.axis.point(0.025, "x"), 
       y=relative.axis.point(0.875, "y"),
       label=paste0("(", LETTERS[seq(2,12,2)[n]], ")"), adj=0,
       font=2)
  
  close.screen(n+6)
  
})

close.screen(all.screens=TRUE)

dev.off()

}