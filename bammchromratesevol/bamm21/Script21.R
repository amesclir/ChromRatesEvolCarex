library(BAMMtools)
##library(coda)

## load data
mytree <- read.tree("./my_tree21.tree")
mcmcout <- read.csv("./mcmc_out21.txt")

## create edata
edata <- getEventData(mytree, eventdata = "./event_data21.txt", burnin=0.15, type = "trait")


#### Check convergence
burnstart <- floor(0.15 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

#effectiveSize(postburn$N_shifts)
#effectiveSize(postburn$logLik)



### Shift probabilities
shift_probs <- summary(edata)
shift_probs

### Bayes factors
bfmat <- computeBayesFactors(postburn, expectedNumberOfShifts=10, burnin=0.15)
bfmat

#### PLOT CREDIBLE SHIFTS
css <- credibleShiftSet(edata, expectedNumberOfShifts=10, threshold=5, set.limit = 0.95)
css

### PLOT BEST SHIFT
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=10)
best


MarginalBranchRateMatrix <- getMarginalBranchRateMatrix(best)
branch_matrix <- MarginalBranchRateMatrix$beta_branch_matrix
branch_matrix
write.csv(branch_matrix, file = "branch_matrix21.csv")
