##########################################
######### Preamble
##########################################

## A script to perform the random walks sensitivity analysis described in
## "Lifetime Mobility of an Arctic Woolly Mammoth"
## by Wooller, Bataille, ... and Willis (2021, Science)

# Authors: Clement Bataille and Amy Willis, 2020
# Software contact: Amy Willis, adwillis@uw.edu

# This script follows `isotope_guided_walk.R`. You will need to load that script first!

##########################################
######### Run the main analyses
##########################################

# This is the main analysis
# walks_d_100_iso_2 <- mclapply(X = 1:20000,
#                               FUN = run_walk_screened,
#                               dist_threshold = 100,
#                               iso_threshold = 2,
#                               mc.cores = 7)
# write_csv(x=walks_d_100_iso_2 %>% bind_rows, path="../output/walks_d_100_iso_2.csv")



##########################################
######### Run the supplementary analyses
##########################################

# walks_d_100_iso_1 <- mclapply(X = 20001:40000,
#                               FUN = run_walk_screened,
#                               dist_threshold = 100,
#                               iso_threshold = 1,
#                               mc.cores = 7)
# write_csv(x=walks_d_100_iso_1 %>% bind_rows, path="../output/walks_d_100_iso_1.csv")

# walks_d_100_iso_3 <- mclapply(X = 40001:60000,
#                               FUN = run_walk_screened,
#                               dist_threshold = 100,
#                               iso_threshold = 3,
#                               mc.cores = 7)
# write_csv(x=walks_d_100_iso_3 %>% bind_rows,
#           path="../output/walks_d_100_iso_3.csv")

# walks_d_25_iso_2 <- mclapply(X = 60001:80000,
#                              FUN = run_walk_screened,
#                              dist_threshold = 25,
#                              iso_threshold = 2,
#                              mc.cores = 7)
# write_csv(x=walks_d_25_iso_2 %>% bind_rows, path="../output/walks_d_25_iso_2.csv")


# walks_d_50_iso_2 <- mclapply(X = 80001:100000,
#                              FUN = run_walk_screened,
#                              dist_threshold = 50,
#                              iso_threshold = 2,
#                              mc.cores = 7)
# write_csv(x=walks_d_50_iso_2 %>% bind_rows, path="../output/walks_d_50_iso_2.csv")
