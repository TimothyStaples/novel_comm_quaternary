test.richness.effect <- function(whole.novel,
                                 subsample.novel,
                                 name,
                                 sqrt.mat=FALSE){
  
div.novel <- function(novel.list){  
  
temp <- do.call("rbind", lapply(unique(novel.list$sites), function(site){
    
    raw.ssmat <- order.ssmat(novel.list$raw.ssmats[[site]])
    prop.ssmat <- order.ssmat(novel.list$prop.ssmats[[site]])
    
    if(sqrt.mat){
    raw.ssmat <- raw.ssmat^2
    }
    
    results <- data.frame(site= site,
               bins = rownames(raw.ssmat),
               count = rowSums(raw.ssmat),
               t(estimateR(round(raw.ssmat))),
               bin.lag = c(NA, abs(diff(as.numeric(rownames(raw.ssmat))))))
    
    novels <- novel.list$novel[[site]]
    
    merge(x=results, y=novels[,c("bins", "cat")],
          by.x="bins", by.y="bins",
          all.x=TRUE, all.y=FALSE, sort=FALSE)  
    
  }))

return(temp)
                      
}

all.div.novel <- div.novel(whole.novel)
sub.div.novel <- div.novel(subsample.novel)
# run models 

richness.prob.models <- lapply(list(all.div.novel, 
                                    sub.div.novel), function(novel.df){
  
  novel.df$novel = novel.df$cat == "novel"
  novel.df$log.count <- log(novel.df$count)
  novel.df$scale.count <- scale(novel.df$log.count)
  # 
  # temp.model <- glmer(novel ~ scale.count + I(scale.count^2) + (1|site), data=novel.df,
  #                     family=binomial)
  # 
  temp.model <- gam(novel ~ s(scale.count), data=novel.df,
                      family=binomial)
  
  return(list(model=temp.model,
              scale.attr = c(attr(novel.df$scale.count, "scaled:center"),
                             attr(novel.df$scale.count, "scaled:scale"))))
})

#                         Plot ####
require(merTools)

rich.prob.preds <- lapply(1:2,
                          function(n){
                            
                          temp = richness.prob.models[[n]]$model
                          scaled.pars <- richness.prob.models[[n]]$scale.attr
                            
                          pred.df <- data.frame(scale.count = seq(min(temp$model$scale.count),
                                                                 max(temp$model$scale.count),
                                                                 len=200))
                            
                          pred.df$log.count = (pred.df$scale.count * scaled.pars[2]) + scaled.pars[1]
                          pred.df$count = exp(pred.df$log.count)
                          
                          pred.df$site <- temp$model$site[1]
                            
                          # pred.df <- cbind(pred.df,
                          #                  sum.fit = predict(temp,
                          #                                    newdata=pred.df,
                          #                                    re.form=NA),
                          #                  predictInterval(temp, 
                          #                                  newdata = pred.df,
                          #                                  n.sims = 999,
                          #                                  which="fixed",
                          #                                  level=0.95,
                          #                                  include.resid.var=FALSE,
                          #                                  type="linear.prediction"))
                          
                          pred.df <- cbind(pred.df,
                                           as.data.frame(predict(temp,
                                                               newdata=pred.df,
                                                               se.fit=TRUE)))
                          pred.df$upr <- pred.df$fit + 1.96 * pred.df$se.fit
                          pred.df$lwr <- pred.df$fit - 1.96 * pred.df$se.fit
                          
                            return(pred.df)
                          })

pdf(date.wrap("./plots/prob by richness", ".pdf"), 
    height=3.5, width=7, useDingbats=FALSE)

par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(2.5,3.5,1,4), 
    tcl=-0.25, mgp=c(3,0.5,0), las=1, ps=10)

lapply(1:2, function(n){
  
  data <- richness.prob.models[[n]]$model$model
  preds <- rich.prob.preds[[n]]
  
  plot(x=NULL, y=NULL, xlim=c(min(preds$count), max(preds$count)), ylim=c(-0.01,0.5), yaxs="i",
       xlab="", ylab="", xaxt="n", log="x", yaxt="n")
  
  if(n==1){
    mtext(side=2, line=2.5, text="Probability", las=0, cex=0.8)
    axis(side=2, mgp=c(3,0.5,0))
  } else {
    axis(side=2, mgp=c(3,0.5,0), labels=NA)
    
  }
  
  mtext(side=1, line=1.25, text="Pollen count", cex=0.8)
  axis(side=1, mgp=c(3,0.2,0), at=c(1,10,100,1000,10000))
  axis(side=1, at=c(seq(1,10,1), seq(10,100,10), seq(100,1000,100),
                    seq(1000,10000,1000), seq(10000,100000,10000)), tcl=-0.125, labels=NA)
  
  per.ts.prop <- tapply(data$novel,
                        data$scale.count,
                        function(x){sum(x) / length(x)})
  
  sample.size <- tapply(data$novel,
                        data$scale.count,
                        length)
  
  
  raw.points <- exp(as.numeric(names(per.ts.prop)) * attr(data$scale.count, "scaled:scale") +
                attr(data$scale.count, "scaled:center"))
  
  novel.col <- col2rgb("orange")/255
  
  points(per.ts.prop ~ raw.points,
         cex = log(sample.size+1) / max(log(sample.size))*1.5,
         pch=16, col=rgb(novel.col[1],
                         novel.col[2],
                         novel.col[3],0.25), lwd=0.1)
  
  lines(x=preds$count,
        y=plogis(preds$upr),
        col="black", lwd=1, lty="31")
  
  lines(x=preds$count,
        y=plogis(preds$lwr),
        col="black", lwd=1, lty="31")
  
  lines(x=preds$count,
        y=plogis(preds$fit),
        col="black", lwd=2)
  
  text(x=exp(relative.axis.point(0, "x")) + c(0, 90)[n],
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(", LETTERS[n], ") ", 
                      c("All sampling bins",
                        "Only pollen counts > 99 and < 9,999")[n]), adj=0, font=2)
  
  if(n==2){
    par(xpd=NA)
    legend(x=75000, y=relative.axis.point(0.5, "y"),
           legend=c("1", "10", "100", "500"),
           pch=21, pt.bg="grey90", bty="n", yjust=0.5,
           pt.cex=log(c(1,10,100,500)+1) / max(log(sample.size))*1.5,
           x.intersp=0.75, y.intersp=0.8, title="Sample\nsize")
    par(xpd=FALSE)
  }
  
  box()
})

dev.off()

}