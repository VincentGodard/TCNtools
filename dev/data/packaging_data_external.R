library("readODS")

# data accessible to the users

tcn_depth_profiles = read_ods("dev/data/data_external/profiles.ods",sheet="data")

# reference data for eulerian calculations ----
## Production and decay parameters ----
prm = matrix(c( 4.01  , 0.012 , 0.039  , log(2)/1.387e6,
                27.93 , 0.84  , 0.081 , log(2)/0.708e6,
                12.24  , 3.31  , 0     ,log(2)/5730),
             nrow = 4,ncol=3 )
colnames(prm) <- c("Be10","Al26","C14") # we just give names to the columns of the matrix
rownames(prm) <- c("Pspal","Pstop","Pfast","lambda")  # we just give names to the rows of the matrix

## references values for production and decay ---
Lambda = c(160,1500,4320) # g/cm2
names(Lambda) <- c("Lspal","Lstop","Lfast") # we just give names to the element of the vector

#

# export -----
usethis::use_data(prm,Lambda,tcn_depth_profiles, internal = FALSE, overwrite = TRUE)

