rm(list=ls())
library(DevTreatRules)
source("Functions/Simulations_Chapter3.R")
source("Functions/Simulations_Chapter4.R")
Expit <- function(x) exp(x) / (1 + exp(x))

# Plots
## simulate one data-generating realization
set.seed(123)
one.realization <- FormatSimulation(n=500, scenario="chapter3",
                                          shift.for.continuous="min",
                                          unif.min=0, unif.max=2, noise.mu=0, noise.sigma=1,
                                          epsilon.mu=0, epsilon.sigma=1)
one.realization.test.df <- FormatSimulation(n=10000, scenario="chapter3",
                                     shift.for.continuous="min",
                                     unif.min=0, unif.max=2, noise.mu=0, noise.sigma=1,
                                     epsilon.mu=0, epsilon.sigma=1)$one.df
## get empirical average (incorrect) and ``population average'' (correct) of response-curves
population.optimal <- GetPopulationCurves(X.signal=one.realization.test.df[, "input_1"],
                             desirable.response=TRUE)

## make plots
png(paste0("Plots/population_averaging_", Sys.Date(), ".png"))
par(mar=c(5.1,5.1,4.1,2.1))
plot(one.realization$one.design$X[one.realization$one.design$T == 0 & one.realization$one.design$L == 0, 1],
      Expit(one.realization$one.response$S[one.realization$one.design$T == 0 & one.realization$one.design$L == 0]),
      pch=3, col=adjustcolor("blue", alpha=0.5),
      ylim=c(0, 1), xlab=expression(X[1]), ylab="P(Y | X, L, T)",
      cex.axis=1.5, cex.lab=2)
points(one.realization$one.design$X[one.realization$one.design$T == 0 & one.realization$one.design$L == 1, 1],
       Expit(one.realization$one.response$S[one.realization$one.design$T == 0 & one.realization$one.design$L == 1]),
       pch=1, col=adjustcolor("blue", alpha=0.5))
points(one.realization$one.design$X[one.realization$one.design$T == 1 & one.realization$one.design$L == 0, 1],
       Expit(one.realization$one.response$S[one.realization$one.design$T == 1 & one.realization$one.design$L == 0]),
       pch=3, col=adjustcolor("orange", alpha=0.5))
points(one.realization$one.design$X[one.realization$one.design$T == 1 & one.realization$one.design$L == 1, 1],
       Expit(one.realization$one.response$S[one.realization$one.design$T == 1 & one.realization$one.design$L == 1]),
       pch=1, col=adjustcolor("orange", alpha=0.5))
lines(one.realization.test.df[, "input_1"][population.optimal$ord.X.signal], population.optimal$correct.prob.control, lwd=5, col="blue")
lines(one.realization.test.df[, "input_1"][population.optimal$ord.X.signal], population.optimal$correct.prob.treat, lwd=5, col="orange")
dev.off()

png(paste0("Plots/sample_averaging_", Sys.Date(), ".png"))
par(mar=c(5.1,5.1,4.1,2.1))
plot(one.realization$one.design$X[one.realization$one.design$T == 0 & one.realization$one.design$L == 0, 1],
      Expit(one.realization$one.response$S[one.realization$one.design$T == 0 & one.realization$one.design$L == 0]),
      pch=3, col=adjustcolor("blue", alpha=0.5),
      ylim=c(0, 1), xlab=expression(X[1]), ylab="P(Y | X, L, T)",
      cex.axis=1.5, cex.lab=2)
points(one.realization$one.design$X[one.realization$one.design$T == 0 & one.realization$one.design$L == 1, 1],
       Expit(one.realization$one.response$S[one.realization$one.design$T == 0 & one.realization$one.design$L == 1]),
       pch=1, col=adjustcolor("blue", alpha=0.5))
points(one.realization$one.design$X[one.realization$one.design$T == 1 & one.realization$one.design$L == 0, 1],
       Expit(one.realization$one.response$S[one.realization$one.design$T == 1 & one.realization$one.design$L == 0]),
       pch=3, col=adjustcolor("orange", alpha=0.5))
points(one.realization$one.design$X[one.realization$one.design$T == 1 & one.realization$one.design$L == 1, 1],
       Expit(one.realization$one.response$S[one.realization$one.design$T == 1 & one.realization$one.design$L == 1]),
       pch=1, col=adjustcolor("orange", alpha=0.5))
lines(one.realization.test.df[, "input_1"][population.optimal$ord.X.signal], population.optimal$incorrect.prob.control, lwd=5, col="blue")
lines(one.realization.test.df[, "input_1"][population.optimal$ord.X.signal], population.optimal$incorrect.prob.treat, lwd=5, col="orange")
dev.off()

## make common legend
png(paste0("Plots/common_legend_", Sys.Date(), ".png"), width=480, height=480)
plot(0, 0, type = "n", ann = F, axes = F)
legend("center", cex=1.5, legend=c("standard of care, L=1", "treatment, L=1", "standard of care, L=0", "treatment, L=0", "averaging, standard of care", "averaging, treatment"),
                         col=c("blue", "orange", "blue", "orange", "blue", "orange"), pch=c(1, 1, 3, 3, NA, NA), lty=c(NA, NA, NA, NA, "solid", "solid"),
                         lwd=c(NA, NA, NA, NA, 3, 3))
dev.off()

# Run simulations (shown here for smaller n.reps than used in paper to speed up the run-time; just change n.reps argument to 10,000 to match results from paper)
# set up
n.reps <- 10
test.n <- 10000
vec.training.n <- c(50, 100, 200, 500, 1000)
mat.simulation.means <- matrix(NA, nrow=5, ncol=length(vec.training.n))
colnames(mat.simulation.means) <- paste0("n_", vec.training.n)
rownames(mat.simulation.means) <- c("split.regression.IPW", "split.regression.naive", "optimal.rule", "treating.all", "treating.none")
mat.simulation.sds <- matrix(NA, nrow=5, ncol=length(vec.training.n))
colnames(mat.simulation.sds) <- paste0("n_", vec.training.n)
rownames(mat.simulation.sds) <- c("split.regression.IPW", "split.regression.naive", "optimal.rule", "treating.all", "treating.none")

for (i in 1:length(vec.training.n)) {
    set.seed(i)
    one.training.n <- vec.training.n[i]
    print(paste0("training set sample size - ", one.training.n))
    mat.mean.response <- matrix(NA, nrow=nrow(mat.simulation.means), ncol=n.reps)
    rownames(mat.mean.response) <- rownames(mat.simulation.means)
    for (b in 1:n.reps) {
        print(paste0("replication - ", b))
       # generate training dataset
        one.training.df <- FormatSimulation(n=one.training.n,
                                                          scenario="chapter3",
                                                          shift.for.continuous="min",
                                                          unif.min=0, unif.max=2, noise.mu=0, noise.sigma=1,
                                                          epsilon.mu=0, epsilon.sigma=1)$one.df
        # build rules
        ## split regression (correct IPW weights)
        rule.split.regression <- BuildRule(data=one.training.df,
                                           study.design="observational",
                                           prediction.approach="split.regression",
                                           name.outcome="no_relapse",
                                           type.outcome="binary",
                                           name.treatment="intervention",
                                           names.influencing.treatment="good_prognosis",
                                           names.influencing.rule="input_1",
                                           desirable.outcome=TRUE,
                                           rule.method="glm.regression",
                                           propensity.method="logistic.regression")
        ## split-regression (incorrect uniform weights)
        rule.split.regression.naive <- BuildRule(data=one.training.df,
                                                 study.design="naive",
                                                 prediction.approach="split.regression",
                                                 name.outcome="no_relapse",
                                                 type.outcome="binary",
                                                 name.treatment="intervention",
                                                 names.influencing.treatment="good_prognosis",
                                                 names.influencing.rule="input_1",
                                                 desirable.outcome=TRUE,
                                                 rule.method="glm.regression",
                                                 propensity.method="logistic.regression")
        # evaluate rules
        ## generate test set
        one.test.df <- FormatSimulation(n=test.n, scenario="chapter3",
                                                     shift.for.continuous="min",
                                                     unif.min=0, unif.max=2, noise.mu=0, noise.sigma=1,
                                                     epsilon.mu=0, epsilon.sigma=1)$one.df
        ## apply split-regression rule (correct IPW weights) to test set and evaluate mean outcome probability
        predictions.split.regression  <- PredictRule(rule.split.regression,
                                                                new.X=one.test.df,
                                                                desirable.outcome=TRUE)
        population.split.regression <- GetPopulationCurves(X.signal=one.test.df[, "input_1"],
                                                                            idx.test.positive=as.logical(predictions.split.regression),
                                                                             desirable.response=TRUE)
        ## apply naive version of split-regression rule (incorrect uniform weights) to test set and evaluate mean outcome probability
        predictions.split.regression.naive  <- PredictRule(rule.split.regression.naive,
                                                                new.X=one.test.df,
                                                                desirable.outcome=TRUE)
        population.split.regression.naive <- GetPopulationCurves(X.signal=one.test.df[, "input_1"],
                                                                            idx.test.positive=as.logical(predictions.split.regression.naive),
                                                                             desirable.response=TRUE)
        ## compute mean response probability for optimal rule (just leave idx.test.positive argument empty)
        population.optimal <- GetPopulationCurves(X.signal=one.test.df[, "input_1"],
                                                                   desirable.response=TRUE)
        # store results
        mat.mean.response["split.regression.IPW", b] <- population.split.regression$mean.response.probability.under.specified.rule
        mat.mean.response["split.regression.naive", b] <- population.split.regression.naive$mean.response.probability.under.specified.rule
        mat.mean.response["optimal.rule", b] <- population.optimal$mean.response.probability.under.optimal.rule
        mat.mean.response["treating.all", b] <- population.optimal$mean.response.probability.treating.all
        mat.mean.response["treating.none", b] <- population.optimal$mean.response.probability.treating.none
    }
    mat.simulation.means[, i] <- apply(mat.mean.response, 1, mean)
    mat.simulation.sds[, i] <- apply(mat.mean.response, 1, sd)
}
write.csv(mat.simulation.means, paste0("Results/simulation_means_", n.reps, "_reps_", Sys.Date(), ".csv"))
write.csv(mat.simulation.sds, paste0("Results/simulation_sds_", n.reps, "_reps_", Sys.Date(), ".csv"))
