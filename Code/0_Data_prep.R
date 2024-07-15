
count_k <- paste0('count_', species.no)
nstop_k <- paste0('nstop_', species.no)
doy_minus_min_k <- paste0('doy_minus_min_', species.no)

# size of pixel
size_A <- 20*20

load(paste0("./Data/all-data.RData"))

# NAO Data
dat <- read.table(paste0("./Data/nao.txt"), header = FALSE,
                  sep = "\t")
dat$V3[which(dat$V1==1+2000 & dat$V2== 3)]


load(paste0("./Data/all-data.RData"))

n.loc <- dim(df_routes)[1]
n.year <- length(unique(df_pix_year$year))
n.replicate_gev = max(df_pix_year$year)- min(df_pix_year$year)+1
n.loc <- nrow(unique(cbind(df_pix_year$x, df_pix_year$y)) )


data_list <- list()

trans_func <- function(x){-log(-(x) / 366)}
trans_func_inv <- function(x){-366*exp(-x)}

for (indx_year in 1:21){
  data_list$centroidassignment_list[[indx_year]] <- 1:n.loc
  #covariate
  #NAO DATA
  data_list$data[[indx_year]] <- rbind(trans_func(df_pix_year[,doy_minus_min_k][which(df_pix_year$year==2000+indx_year)]),
                                       rep(dat$V3[which(dat$V1==indx_year+2000 & dat$V2== 3)],n.loc) )

}


dat.mat <- matrix(ncol=n.loc, nrow=n.replicate_gev)

for (indx_year in 1:21){
  dat.mat[indx_year,] <- trans_func(df_pix_year[,doy_minus_min_k][which(df_pix_year$year==2000+indx_year)])
}


indx_nona <- colSums(is.na(dat.mat))<nrow(dat.mat)
#sum(indx_nona)

which.loc.all.na <- apply(dat.mat, 2, function(x) median(x, na.rm=TRUE))

gev_dat <- apply(dat.mat[,indx_nona], 2, function(x) mean(x, na.rm=TRUE))
gev_dat <- gev_dat[which(gev_dat>-0.5)]
#op_gev_fit <- ismev::gev.fit(gev_dat)

# hist(trans_func_inv(gev_dat), breaks=100)
# hist((gev_dat), breaks=100)


coords <- unique(cbind(df_pix_year$x, df_pix_year$y))

n.sites <- nrow(coords)

centroids_count_list <- centroids_count_list_2 <- centroids_count_list_bin <- list()

for (indx_year in 1:21){
  # set the NAs to 0!
  # covariate
  centroids_count_list[[indx_year]] <- cbind(coords,df_pix_year$count[which(df_pix_year$year==2000+indx_year)],
                                             rep(indx_year,n.loc ))
}


#sum(df_pix_year[,count_k][which(df_pix_year$year==2000+21)])/sum(df_pix_year$count[which(df_pix_year$year==2000+21)])

# Binomial component
# Change to effort later on
durations <- matrix(ncol=n.loc, nrow=21)
n_s <- matrix(ncol=n.loc, nrow=21)
n_s_i <- matrix(ncol=n.loc, nrow=21)

#min(df_pix_year$duration_median,na.rm=TRUE)
#hist(df_pix_year$duration_median, breaks=100)
#table(df_pix_year$duration_median)

df_sub <- df_pix_year[df_pix_year$count>0,]


for (indx_year in 1:21){
  duration <- df_pix_year$duration_median[which(df_pix_year$year==2000+indx_year)]
  durations[indx_year,] <- duration
  #duration[is.na(duration)] <- 0.1
  duration[which(duration==0)] <- NA

  n_s[indx_year,] <- df_pix_year$count[which(df_pix_year$year==2000+indx_year)]
  n_s_i[indx_year,] <- df_pix_year[,count_k][which(df_pix_year$year==2000+indx_year)]

  centroids_count_list_bin[[indx_year]] <- cbind(coords, df_pix_year[,count_k][which(df_pix_year$year==2000+indx_year)],
                                                 df_pix_year$count[which(df_pix_year$year==2000+indx_year)],
                                                 duration)
}


#apply( durations , 2, function(x) mean(x, na.rm=TRUE))
#apply( n_s , 2, function(x) mean(x, na.rm=TRUE))

#mean(n_s_i/n_s, na.rm=TRUE)



log_like_poissson <- function(y) {
  n <- length(y)
  function(mu) {
    log(mu) * sum(y) - n *mu - sum(lfactorial(y))
  }
}
# Optimization of the loglikelihood
op_pp <- optimise(log_like_poissson(df_pix_year$count),
         interval = c(0, 500),
         maximum = TRUE)


head(df_routes_year)
df_bbs_consider <- cbind(df_routes_year[,1:3],df_routes_year[,count_k], df_routes_year[,nstop_k])
colnames(df_bbs_consider) <- c('id_route', 'route', 'year', 'count', 'nstop')

n.route <- nrow(routes)
df_routes[is.na(df_routes)] <- 0

route_length <- apply(df_routes[,-c(1,2)], 2, sum)


indx.remove <- which(route_length==0)

dat_bbs <- matrix(nrow= n.replicate_gev, ncol=n.route)

for (indx_year in 2001:2021){
  id_lines = which(df_bbs_consider$year == indx_year)
  # merge
  route_id <- data.frame(id_route=(1:n.route))#[-indx.remove] )

  df_bbs_consider_2 <- merge(route_id, df_bbs_consider[id_lines,], by=c('id_route'), all.x=TRUE)

  # in some years, aggregate multiple routes' counts
  df_bbs_consider_2 <- aggregate(cbind(df_bbs_consider_2$count, df_bbs_consider_2$nstop), by=list(df_bbs_consider_2$id_route), sum)
  colnames(df_bbs_consider_2) <- c('id_route', 'count', 'nstop')

  ## turn into 0-1 data
  ind_cnt <- numeric(length=length(df_bbs_consider_2$count))
  ind_cnt[which(df_bbs_consider_2$count>0)] <- 1

  dat_bbs[indx_year-2000,] <- df_bbs_consider_2$count/df_bbs_consider_2$nstop

  centroids_count_list_2[[indx_year-2000]] <- cbind(df_bbs_consider_2$id_route, df_bbs_consider_2$count, df_bbs_consider_2$nstop)
  #centroids_count_list_2[[indx_year-2000]] <- cbind(df_bbs_consider_2$id_route, ind_cnt, df_bbs_consider_2$nstop)
  # set the NAs to 0!
}

