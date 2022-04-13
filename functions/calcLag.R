calcLag <- function(novel.list, 
                            env.data, 
                            env.var, 
                            local.env.data,
                            local.lag,
                            global.lag){
  
# Prepare modelling data ####
print("Prepping data")


# aggregate 50 year temp change estimates to 200 year sampling bins, take mean
colnames(env.data)[colnames(env.data) == env.var] = "temp"

# global data
env.agg <- bin.env.data(env.data,
                        env.var = "temp",
                        bin.width=200,
                        lims=c(-100,25000))
env.agg <- env.agg[order(env.agg$bin, decreasing=TRUE), ]
  
# site data
local.agg <- t(apply(local.env.data, 1, function(x){
  
  temp.agg <- cut(as.numeric(colnames(local.env.data)), breaks=seq(-100,25000,200))
  
  agg.char <- as.character(temp.agg)
  bin <- as.numeric(substr(agg.char, regexpr(",", agg.char)+1,
                                    nchar(agg.char)-1)) - 0.5*200
  tapply(x[-1], bin[-1], mean, na.rm=TRUE)
}))
local.agg <- local.agg[,order(as.numeric(colnames(local.agg)), decreasing=TRUE)]

# now difference temps based on required lag
env.lag <- diff(env.agg$env, differences=1, lag=unlist(global.lag))

env.agg$diff.env <- c(rep(NA, global.lag), env.lag)

local.lag.df <- t(apply(local.agg, 1, function(x){
  tempLag <- c(rep(NA, local.lag), diff(x, differences=1, lag=unlist(local.lag)))
}))
dimnames(local.lag.df) = dimnames(local.agg)

# add in local temp data
local.lag.long <- long_form(dataTable = local.lag.df,
                            data.cols = matrix(rownames(local.lag.df), ncol=1, dimnames=list(rownames(local.lag.df), NA)),
                            category.cols = local.lag.df)
colnames(local.lag.long) = c("site", "bin", "local.diff")

comb.lag <- merge(local.lag.long, env.agg[,c("bin","diff.env")],
                  by.x="bin", by.y="bin",
                  all.x=TRUE, all.y=FALSE, sort=FALSE)
head(comb.lag)
colnames(comb.lag) = c("bin", "site", paste0("localLag", local.lag),  paste0("globalLag",global.lag))

return(comb.lag)

}