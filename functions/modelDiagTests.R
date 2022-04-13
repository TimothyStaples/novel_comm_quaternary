modelDiagTests <- function(model, time, data, spat.iter, 
                           spat.sub.size, sautoSubsample = FALSE){

  require(ape)
  
  modelSim <- simulateResiduals(model)
  
  unif.test = testUniformity(modelSim)
  print(ifelse(abs(unif.test$statistic) > 0.1, 
               "Uniformity not okay :(",
               "Uniformity okay :)"))
  
  disp.test = testDispersion(modelSim)
  print(ifelse(disp.test$p.value <=0.05, 
               "Dispersal not okay :(",
               "Dispersal okay :)"))
  
  time.data <- time
  
  tauto.test <- do.call("rbind", lapply(1:spat.iter, function(n){
  # 
  # while(sum(table(time.offset) > 2) > 0){
  #   time.offset <- time.data + rnorm(length(time.offset), 0, 1)
  #   } 
  subRes <- sapply(split(residuals(model), f=time), function(x){
    
    x[sample(1:length(x),1)]
    
  })

  if(class(model)[1] == "glm"){
  names(subRes) = substr(names(subRes), 1, regexpr("\\.", names(subRes))-1)
    }
  
  a <- dwtest(subRes[match(names(subRes), sort(unique(time)))] ~ sort(unique(time)))
  
  return(data.frame(statistic = a$statistic,
                    method = a$method,
                    alternative = a$alternative,
                    p = a$p.value,
                    data.name = a$data.name))
  }))
  
  comb.comp <- data[!is.na(data$novel),]
  comb.site <- comb.comp[!duplicated(comb.comp$site),]
  ids <- expand.grid(1:nrow(comb.site), 1:nrow(comb.site))
  
  sauto.test <- do.call("rbind", lapply(1:spat.iter, function(n){
  
  subRes <- sapply(split(residuals(model), f=data$site[!is.na(data$novel)]), function(x){
    x[sample(1:length(x),1)]
    })

  dists <- gc.dist(lat1=comb.site$lat[ids[,1]],
                     lat2=comb.site$lat[ids[,2]],
                     lon1 = comb.site$long[ids[,1]],
                     lon2=comb.site$long[ids[,2]])
  dists[dists==0] = 1
  
  dist.mat <- matrix(0, nrow=nrow(comb.site), ncol=nrow(comb.site),
                     dimnames=list(1:nrow(comb.site),
                                   1:nrow(comb.site)))
  dist.mat[as.matrix(ids)] = dists
  
  invDistMat <- 1/dist.mat
  diag(invDistMat) <- 0
  
  MI = ape::Moran.I(subRes, weight = invDistMat, na.rm=TRUE)
  
  return(data.frame(obs=MI$observed,
                    exp=MI$expected,
                    sd=MI$sd,
                    p=MI$p.value))
  }))
  
  return(list(unif = unif.test,
              disp = disp.test,
              tauto = tauto.test,
              sauto = sauto.test))
  
}