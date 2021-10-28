library(brms)

bbs_naive <- readRDS("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/brms_bbs/naive.RDS")
bbs_naive <- add_criterion(bbs_naive, "loo", moment_match = T)
# saveRDS(bbs_naive, "/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/brms_bbs/naive_moment_match.RDS")

bbs_dist <- readRDS("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/brms_bbs/dist.RDS")
bbs_dist <- add_criterion(bbs_dist, "loo", moment_match = T)
# saveRDS(bbs_dist, "/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/brms_bbs/dist_moment_match.RDS")
bbs_dist <- readRDS("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/brms_bbs/dist_moment_match.RDS")

bbs_naive2 <- readRDS("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/brms_bbs/naive2.RDS")
bbs_naive2 <- add_criterion(bbs_naive2, "loo", moment_match = F)

bbs_dist2 <- readRDS("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/brms_bbs/dist2.RDS")
bbs_dist2 <- add_criterion(bbs_dist2, "loo", moment_match = F)


bbs_naive <- readRDS("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/brms_bbs/naive_moment_match.RDS")

loo_compare(bbs_dist, bbs_dist2, criterion = "loo")
