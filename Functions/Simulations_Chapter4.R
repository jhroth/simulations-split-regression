DetailedSummaries <- function(y.test, treatment.test, treatment.rule, true.positive.test, true.negative.test) {
    sensitivity <- mean(treatment.rule[true.positive.test == 1])
    specificity <- mean(treatment.rule[true.negative.test == 1] == 0)
    misclassification.rate <- 1 - mean(treatment.rule == true.positive.test)
    proportion.marker.negative <- mean(treatment.rule == 0)
    return(list("sensitivity"=sensitivity, "specificity"=specificity,
                  "misclassification.rate"=misclassification.rate,
                  "proportion.marker.negative"=proportion.marker.negative))
}

DoAllSimulations <- function(n.reps, scenario, n.training, n.test,
                                                 shift.for.continuous="min",
                                                 unif.min=0, unif.max=2, noise.mu=0, noise.sigma=1, epsilon.mu=0, epsilon.sigma=1, clinical.threshold=0,
                                                 vec.approaches,
                                                 bootstrap.CI=FALSE,
                                                 debugging.DI=FALSE,
                                                 OWL.framework.shift.by.min=TRUE,
                                                 direct.interactions.center.continuous.Y=TRUE,
                                                 direct.interactions.exclude.A.from.penalty) {
    # create empty matrix shells
    mat.prop.pos <- matrix(NA, nrow=length(vec.approaches) + 3, ncol=n.reps)
    mat.ATE.pos <- matrix(NA, nrow=length(vec.approaches) + 3, ncol=n.reps)
    mat.ATE.neg <- matrix(NA, nrow=length(vec.approaches) + 3, ncol=n.reps)
    mat.ABR <- matrix(NA, nrow=length(vec.approaches) + 3, ncol=n.reps)
    mat.exact.ATE.pos <- matrix(NA, nrow=length(vec.approaches) + 3, ncol=n.reps)
    mat.exact.ATE.neg <- matrix(NA, nrow=length(vec.approaches) + 3, ncol=n.reps)
    mat.exact.ABR <- matrix(NA, nrow=length(vec.approaches) + 3, ncol=n.reps)
    mat.misclassification.rate <- matrix(NA, nrow=length(vec.approaches) + 3, ncol=n.reps)
    mat.FPR <- matrix(NA, nrow=length(vec.approaches) + 3, ncol=n.reps)
    mat.TPR <- matrix(NA, nrow=length(vec.approaches) + 3, ncol=n.reps)
    # do n.reps iterations
    for (rep in 1:n.reps) {
        set.seed(rep)
        one.result <- DoOneSimulation(scenario=scenario, n.training=n.training, n.test=n.test,
                                              shift.for.continuous=shift.for.continuous,
                                              unif.min=unif.min, unif.max=unif.max, noise.mu=noise.mu, noise.sigma=noise.sigma, epsilon.mu=epsilon.mu, epsilon.sigma=epsilon.sigma, clinical.threshold=clinical.threshold,
                                              vec.approaches=vec.approaches,
                                              bootstrap.CI=bootstrap.CI,
                                              debugging.DI=debugging.DI,
                                              OWL.framework.shift.by.min=OWL.framework.shift.by.min,
                                              direct.interactions.center.continuous.Y=direct.interactions.center.continuous.Y,
                                              direct.interactions.exclude.A.from.penalty=direct.interactions.exclude.A.from.penalty)
        one.mat.summary <- one.result$mat.summary
        mat.prop.pos[, rep] <- one.mat.summary[, "prop.pos"]
        mat.ATE.pos[, rep] <- one.mat.summary[, "ATE.pos"]
        mat.ATE.neg[, rep] <- one.mat.summary[, "ATE.neg"]
        mat.ABR[, rep] <- one.mat.summary[, "ABR"]
        mat.exact.ATE.pos[, rep] <- one.mat.summary[, "exact.ATE.pos"]
        mat.exact.ATE.neg[, rep] <- one.mat.summary[, "exact.ATE.neg"]
        mat.exact.ABR[, rep] <- one.mat.summary[, "exact.ABR"]
        mat.misclassification.rate[, rep] <- one.mat.summary[, "misclassification.rate"]
        mat.TPR[, rep] <- one.mat.summary[, "TPR"]
        mat.FPR[, rep] <- one.mat.summary[, "FPR"]
    }
    rownames(mat.ATE.pos) <- rownames(one.mat.summary)
    rownames(mat.ATE.neg) <- rownames(one.mat.summary)
    rownames(mat.ABR) <- rownames(one.mat.summary)
    rownames(mat.exact.ATE.pos) <- rownames(one.mat.summary)
    rownames(mat.exact.ATE.neg) <- rownames(one.mat.summary)
    rownames(mat.exact.ABR) <- rownames(one.mat.summary)
    rownames(mat.misclassification.rate) <- rownames(one.mat.summary)
    rownames(mat.TPR) <- rownames(one.mat.summary)
    rownames(mat.FPR) <- rownames(one.mat.summary)
    # compute sample mean and SD for each statistic
    ## mean
    mat.means <- matrix(NA, nrow=length(vec.approaches) + 3, ncol=10)
    rownames(mat.means) <- rownames(one.mat.summary)
    colnames(mat.means) <- c("prop.pos", "ATE.pos", "ATE.neg", "ABR", "exact.ATE.pos", "exact.ATE.neg", "exact.ABR", "misclassification.rate", "TPR", "FPR")
    mat.means[rownames(mat.prop.pos), "prop.pos"] <- apply(mat.prop.pos, 1, mean)
    mat.means[rownames(mat.ATE.pos), "ATE.pos"] <- apply(mat.ATE.pos, 1, function(x) mean(x, na.rm=TRUE))
    mat.means[rownames(mat.ATE.neg), "ATE.neg"] <- apply(mat.ATE.neg, 1, function(x) mean(x, na.rm=TRUE))
    mat.means[rownames(mat.ABR), "ABR"] <- apply(mat.ABR, 1, function(x) mean(x, na.rm=TRUE))
    mat.means[rownames(mat.exact.ATE.pos), "exact.ATE.pos"] <- apply(mat.exact.ATE.pos, 1, function(x) mean(x, na.rm=TRUE))
    mat.means[rownames(mat.exact.ATE.neg), "exact.ATE.neg"] <- apply(mat.exact.ATE.neg, 1, function(x) mean(x, na.rm=TRUE))
    mat.means[rownames(mat.exact.ABR), "exact.ABR"] <- apply(mat.exact.ABR, 1, function(x) mean(x, na.rm=TRUE))
    mat.means[rownames(mat.misclassification.rate), "misclassification.rate"] <- apply(mat.misclassification.rate, 1, mean)
    mat.means[rownames(mat.TPR), "TPR"] <- apply(mat.TPR, 1, mean)
    mat.means[rownames(mat.FPR), "FPR"] <- apply(mat.FPR, 1, mean)

    ## SD
    mat.sds <- matrix(NA, nrow=length(vec.approaches) + 3, ncol=10)
    rownames(mat.sds) <- rownames(one.mat.summary)
    colnames(mat.sds) <- c("prop.pos", "ATE.pos", "ATE.neg", "ABR", "exact.ATE.pos", "exact.ATE.neg", "exact.ABR", "misclassification.rate", "TPR", "FPR")
    mat.sds[rownames(mat.prop.pos), "prop.pos"] <- apply(mat.prop.pos, 1, sd)
    mat.sds[rownames(mat.ATE.pos), "ATE.pos"] <- apply(mat.ATE.pos, 1, function(x) sd(x, na.rm=TRUE))
    mat.sds[rownames(mat.ATE.neg), "ATE.neg"] <- apply(mat.ATE.neg, 1, function(x) sd(x, na.rm=TRUE))
    mat.sds[rownames(mat.ABR), "ABR"] <- apply(mat.ABR, 1, function(x) sd(x, na.rm=TRUE))
    mat.sds[rownames(mat.exact.ATE.pos), "exact.ATE.pos"] <- apply(mat.exact.ATE.pos, 1, function(x) sd(x, na.rm=TRUE))
    mat.sds[rownames(mat.exact.ATE.neg), "exact.ATE.neg"] <- apply(mat.exact.ATE.neg, 1, function(x) sd(x, na.rm=TRUE))
    mat.sds[rownames(mat.exact.ABR), "exact.ABR"] <- apply(mat.exact.ABR, 1, function(x) sd(x, na.rm=TRUE))
    mat.sds[rownames(mat.misclassification.rate), "misclassification.rate"] <- apply(mat.misclassification.rate, 1, sd)
    mat.sds[rownames(mat.TPR), "TPR"] <- apply(mat.TPR, 1, sd)
    mat.sds[rownames(mat.FPR), "FPR"] <- apply(mat.FPR, 1, sd)
    
    ## count NA
    mat.counts <- matrix(NA, nrow=length(vec.approaches) + 3, ncol=10+1)
    rownames(mat.counts) <- rownames(one.mat.summary)
    colnames(mat.counts) <- c("prop.pos.0", "prop.pos.1", "ATE.pos", "ATE.neg", "ABR", "exact.ATE.pos", "exact.ATE.neg", "exact.ABR", "misclassification.rate", "TPR", "FPR")
    mat.counts[rownames(mat.prop.pos), "prop.pos.0"] <- apply(mat.prop.pos, 1, function(x) sum(x == 0))
    mat.counts[rownames(mat.prop.pos), "prop.pos.1"] <- apply(mat.prop.pos, 1, function(x) sum(x == 1))
    mat.counts[rownames(mat.ATE.pos), "ATE.pos"] <- apply(mat.ATE.pos, 1, function(x) sum(is.na(x)))
    mat.counts[rownames(mat.ATE.neg), "ATE.neg"] <- apply(mat.ATE.neg, 1, function(x) sum(is.na(x)))
    mat.counts[rownames(mat.ABR), "ABR"] <- apply(mat.ABR, 1, function(x) sum(is.na(x)))
    mat.counts[rownames(mat.exact.ATE.pos), "exact.ATE.pos"] <- apply(mat.exact.ATE.pos, 1, function(x) sum(is.na(x)))
    mat.counts[rownames(mat.exact.ATE.neg), "exact.ATE.neg"] <- apply(mat.exact.ATE.neg, 1, function(x) sum(is.na(x)))
    mat.counts[rownames(mat.exact.ABR), "exact.ABR"] <- apply(mat.exact.ABR, 1, function(x) sum(is.na(x)))
    mat.counts[rownames(mat.misclassification.rate), "misclassification.rate"] <- apply(mat.misclassification.rate, 1, function(x) sum(is.na(x)))
    mat.counts[rownames(mat.TPR), "TPR"] <- apply(mat.TPR, 1, function(x) sum(is.na(x)))
    mat.counts[rownames(mat.FPR), "FPR"] <- apply(mat.FPR, 1, function(x) sum(is.na(x)))
    
    return(list("one.mat.summary"=one.mat.summary,
                 "mat.prop.pos"=mat.prop.pos,
                "mat.ATE.pos"=mat.ATE.pos,
                "mat.ATE.neg"=mat.ATE.neg,
                "mat.ABR"=mat.ABR,
                "mat.exact.ATE.pos"=mat.exact.ATE.pos,
                "mat.exact.ATE.neg"=mat.exact.ATE.neg,
                "mat.exact.ABR"=mat.exact.ABR,
                "mat.misclassification.rate"=mat.misclassification.rate,
                "mat.FPR"=mat.FPR,
                "mat.TPR"=mat.TPR,
                "mat.means"=mat.means,
                "mat.sds"=mat.sds,
                "mat.counts"=mat.counts))
}

DoOneSimulation <- function(scenario, n.training, n.test,
                                                 shift.for.continuous="min",
                                                 unif.min, unif.max, noise.mu, noise.sigma, epsilon.mu, epsilon.sigma, clinical.threshold=0,
                                                 vec.approaches,
                                                 bootstrap.CI,
                                                 debugging.DI=FALSE,
                                                 OWL.framework.shift.by.min=TRUE,
                                                 direct.interactions.center.continuous.Y=TRUE,
                                                 direct.interactions.exclude.A.from.penalty) {
    # set scenario-specific parameters
    if (scenario %in% c("default", "main.effect.added", "quadratic.interaction.added", "piecewise.cubic.interaction.added", "chapter3")) {
        names.influencing.rule <- "input_1"
        names.influencing.treatment <- "good_prognosis"
        rule.method <- "glm.regression"
        name.outcome <- "no_relapse"
        type.outcome <- "binary"
        name.treatment <- "intervention"
    } else if (scenario %in% c("default.piecewise.constant", "two.way.interactions.piecewise.constant", "three.way.interactions.piecewise.constant", "nonlinear.interactions.piecewise.constant")) {
        names.influencing.rule <- c("input_1", "input_2", "input_3")
        names.influencing.treatment <- "good_prognosis"
        rule.method <- "glm.regression"
        name.outcome <- "no_relapse"
        type.outcome <- "binary"
        name.treatment <- "intervention"
    } else if (scenario %in% c("default.piecewise.constant.continuous", "two.way.interactions.piecewise.constant.continuous",
                                       "three.way.interactions.piecewise.constant.continuous", "nonlinear.interactions.piecewise.constant.continuous")) {
        names.influencing.rule <- c("input_1", "input_2", "input_3")
        names.influencing.treatment <- "good_prognosis"
        rule.method <- "glm.regression"
        name.outcome <- "time_to_relapse"
        type.outcome <- "continuous"
        name.treatment <- "intervention"
    } else if (scenario %in% c("noise.added")) {
        names.influencing.rule <- c("input_1", paste0("noise_", 1:40))
        names.influencing.treatment <- "good_prognosis"
        rule.method <- "lasso" 
        name.outcome <- "no_relapse"
        type.outcome <- "binary"
        name.treatment <- "intervention"
    } else if (scenario %in% c("high.dimensional.noise.added")) {
        names.influencing.rule <- c("input_1",
                                             paste0("noise_", 1:500))
        names.influencing.treatment <- "good_prognosis"
        rule.method <- "lasso" 
        name.outcome <- "no_relapse"
        type.outcome <- "binary"
        name.treatment <- "intervention"
    } else if (scenario %in% c("high.dimensional.noise.and.prognostic.added")) {
        names.influencing.rule <- c("input_1", "input_2", "input_3", "input_4", "input_5", "input_6",
                                              paste0("noise_", 1:500))
        names.influencing.treatment <- "good_prognosis"
        rule.method <- "lasso" 
        name.outcome <- "no_relapse"
        type.outcome <- "binary"
        name.treatment <- "intervention"
    } else if (scenario %in% c("default.continuous", "main.effect.added.continuous", "quadratic.interaction.added.continuous", "piecewise.cubic.interaction.added.continuous")) {
        names.influencing.rule <- "input_1"
        names.influencing.treatment <- "good_prognosis"
        rule.method <- "glm.regression"
        name.outcome <- "time_to_relapse"
        type.outcome <- "continuous"
        name.treatment <- "intervention"
    } else if (scenario %in% c("noise.added.continuous")) {
        names.influencing.rule <- c("input_1", paste0("noise_", 1:40))
        names.influencing.treatment <- "good_prognosis"
        rule.method <- "lasso" 
        name.outcome <- "time_to_relapse"
        type.outcome <- "continuous"
        name.treatment <- "intervention"
    } else if (scenario %in% c("high.dimensional.noise.added.continuous")) {
        names.influencing.rule <- c("input_1",
                                             paste0("noise_", 1:500))
        names.influencing.treatment <- "good_prognosis"
        rule.method <- "lasso" 
        name.outcome <- "time_to_relapse"
        type.outcome <- "continuous"
        name.treatment <- "intervention"
    } else if (scenario %in% c("high.dimensional.noise.and.prognostic.added.continuous")) {
        names.influencing.rule <- c("input_1", "input_2", "input_3", "input_4", "input_5", "input_6",
                                              paste0("noise_", 1:500))
        names.influencing.treatment <- "good_prognosis"
        rule.method <- "lasso" 
        name.outcome <- "time_to_relapse"
        type.outcome <- "continuous"
        name.treatment <- "intervention"
    }
    # generate data
    training.set <- FormatSimulation(scenario=scenario,
                                          n=n.training,
                                          shift.for.continuous=shift.for.continuous,
                                          unif.min=unif.min, unif.max=unif.max, noise.mu=noise.mu, noise.sigma=noise.sigma, epsilon.mu=epsilon.mu, epsilon.sigma=epsilon.sigma,
                                          debugging.DI=debugging.DI)
    test.set <- FormatSimulation(scenario=scenario,
                                          n=n.test,
                                          shift.for.continuous=shift.for.continuous,
                                          unif.min=unif.min, unif.max=unif.max, noise.mu=noise.mu, noise.sigma=noise.sigma, epsilon.mu=epsilon.mu, epsilon.sigma=epsilon.sigma,
                                          debugging.DI=debugging.DI)
    # create matrix to store one set of simulation summaries (prop+, ATE+, ABR, MCR, FPR, TPR) for the different approaches
    mat.summary <- matrix(NA, nrow=length(vec.approaches)+3, ncol=10)
    rownames(mat.summary) <- c(vec.approaches, "optimal.rule", "treat.all", "treat.none")
    colnames(mat.summary) <- c("prop.pos", "ATE.pos", "ATE.neg", "ABR", "exact.ATE.pos", "exact.ATE.neg", "exact.ABR", "misclassification.rate", "FPR", "TPR")
    # fill in matrix
    for (a in 1:length(vec.approaches)) {
        one.approach <- vec.approaches[a]
        if (one.approach == "split.regression.naive") {
            my.study.design <- "naive"
            my.approach <- "split.regression"
        } else {
            my.study.design <- "observational"
            my.approach <- one.approach
        }
        get.one.rule <- BuildRule(data=training.set$one.df,
                                   study.design=my.study.design,
                                   prediction.approach=my.approach,
                                   name.outcome=name.outcome,
                                   type.outcome=type.outcome,
                                   name.treatment=name.treatment,
                                   names.influencing.treatment=names.influencing.treatment,
                                   names.influencing.rule=names.influencing.rule,
                                   desirable.outcome=TRUE,
                                   rule.method=rule.method,
                                   propensity.method="logistic.regression",
                                   type.observation.weights=NULL,
                                   truncate.propensity.score=TRUE,
                                   OWL.framework.shift.by.min=OWL.framework.shift.by.min, 
                                   direct.interactions.center.continuous.Y=direct.interactions.center.continuous.Y,
                                   direct.interactions.exclude.A.from.penalty=direct.interactions.exclude.A.from.penalty)
        predicted.one.rule <- PredictRule(get.one.rule,
                                           new.X=test.set$one.df[, names.influencing.rule, drop=FALSE],
                                           desirable.outcome=TRUE,
                                           return.predicted.response=TRUE)
        predicting.with.oracle <- GetOptimalRule(scenario=scenario, B=predicted.one.rule$recommended.treatment, use.prespecified.data=TRUE,
                                                       X=test.set$one.design$X, T=test.set$one.design$T, L=test.set$one.design$L,
                                                       unif.min=unif.min, unif.max=unif.max, noise.mu=noise.mu, noise.sigma=noise.sigma)
        two.by.two.table <- table(predicted.one.rule$recommended.treatment, predicting.with.oracle$one.optimal.rule.evaluation$recommended.treatment)
        evaluate.one.rule <- EvaluateRule(data=test.set$one.df,
                                           BuildRule.object=get.one.rule,
                                           study.design="observational",
                                           name.outcome=name.outcome,
                                           type.outcome=type.outcome,
                                           desirable.outcome=TRUE,
                                           clinical.threshold=clinical.threshold,
                                           name.treatment=name.treatment,
                                           names.influencing.treatment=names.influencing.treatment,
                                           names.influencing.rule=names.influencing.rule,
                                           propensity.method="logistic.regression",
                                           bootstrap.CI=bootstrap.CI)
        mat.summary[one.approach, "prop.pos"] <- evaluate.one.rule$n.test.positives / nrow(test.set$one.df)
        mat.summary[one.approach, "ATE.pos"] <- evaluate.one.rule$ATE.test.positives
        mat.summary[one.approach, "ATE.neg"] <- evaluate.one.rule$ATE.test.negatives
        mat.summary[one.approach, "ABR"] <- evaluate.one.rule$ABR
        one.result.two.by.two <- SummarizeTable(two.by.two.table)
        mat.summary[one.approach, "exact.ATE.pos"] <- predicting.with.oracle$one.estimated.rule.evaluation$ATE.positives
        mat.summary[one.approach, "exact.ATE.neg"] <- predicting.with.oracle$one.estimated.rule.evaluation$ATE.negatives
        mat.summary[one.approach, "exact.ABR"] <- predicting.with.oracle$one.estimated.rule.evaluation$ABR
        mat.summary[one.approach, "misclassification.rate"] <- one.result.two.by.two$misclassification.rate
        mat.summary[one.approach, "FPR"] <- one.result.two.by.two$FPR
        mat.summary[one.approach, "TPR"] <- one.result.two.by.two$TPR
    }
    #comparison rules
    ## optimal rule
    one.comparison.rules <- GetComparisonRules(scenario=scenario,
                                                     B.optimal=predicting.with.oracle$one.optimal.rule.evaluation$recommended.treatment,
                                                     test.df=test.set$one.df,
                                                     study.design="observational",
                                                     name.outcome=name.outcome, type.outcome=type.outcome, desirable.outcome=TRUE,
                                                     name.treatment=name.treatment,
                                                     names.influencing.treatment=names.influencing.treatment,
                                                     names.influencing.rule=names.influencing.rule,
                                                     propensity.method="logistic.regression")
    
    mat.summary["optimal.rule", "ATE.pos"] <- one.comparison.rules$evaluate.optimal.rule$ATE.test.positives
    mat.summary["optimal.rule", "ATE.neg"] <- one.comparison.rules$evaluate.optimal.rule$ATE.test.negatives
    mat.summary["optimal.rule", "ABR"] <- one.comparison.rules$evaluate.optimal.rule$ABR
    mat.summary["optimal.rule", "prop.pos"] <- predicting.with.oracle$one.optimal.rule.evaluation$prop.positives
    mat.summary["optimal.rule", "exact.ATE.pos"] <- predicting.with.oracle$one.optimal.rule.evaluation$ATE.positives
    mat.summary["optimal.rule", "exact.ATE.neg"] <- predicting.with.oracle$one.optimal.rule.evaluation$ATE.negatives
    mat.summary["optimal.rule", "exact.ABR"] <- predicting.with.oracle$one.optimal.rule.evaluation$ABR
    mat.summary["optimal.rule", "misclassification.rate"] <- 0
    mat.summary["optimal.rule", "FPR"] <- 0
    mat.summary["optimal.rule", "TPR"] <- 1

    ## treating none
    B.treat.none <- rep(0, n.test)
    predicting.with.oracle.treat.none <- GetOptimalRule(scenario=scenario, B=B.treat.none, use.prespecified.data=TRUE,
                                                   X=test.set$one.design$X, T=test.set$one.design$T, L=test.set$one.design$L,
                                                   unif.min=0, unif.max=2, noise.mu=0, noise.sigma=1)
    mat.summary["treat.none", "ATE.pos"] <- one.comparison.rules$evaluate.treat.none$ATE.test.positives
    mat.summary["treat.none", "ATE.neg"] <- one.comparison.rules$evaluate.treat.none$ATE.test.negatives
    mat.summary["treat.none", "ABR"] <- one.comparison.rules$evaluate.treat.none$ABR
    mat.summary["treat.none", "prop.pos"] <- predicting.with.oracle.treat.none$one.estimated.rule.evaluation$prop.positives
    mat.summary["treat.none", "exact.ATE.pos"] <- predicting.with.oracle.treat.none$one.estimated.rule.evaluation$ATE.positives
    mat.summary["treat.none", "exact.ATE.neg"] <- predicting.with.oracle.treat.none$one.estimated.rule.evaluation$ATE.negatives
    mat.summary["treat.none", "exact.ABR"] <- predicting.with.oracle.treat.none$one.estimated.rule.evaluation$ABR
    two.by.two.table.treat.none <- table(B.treat.none, predicting.with.oracle$one.optimal.rule.evaluation$recommended.treatment)
    one.result.two.by.two.treat.none <- SummarizeTable(two.by.two.table.treat.none)
    mat.summary["treat.none", "misclassification.rate"] <- one.result.two.by.two.treat.none$misclassification.rate
    mat.summary["treat.none", "FPR"] <- 0
    mat.summary["treat.none", "TPR"] <- 0

    ## treating none
    B.treat.all <- rep(1, n.test)
    predicting.with.oracle.treat.all <- GetOptimalRule(scenario=scenario, B=B.treat.all, use.prespecified.data=TRUE,
                                                   X=test.set$one.design$X, T=test.set$one.design$T, L=test.set$one.design$L,
                                                   unif.min=0, unif.max=2, noise.mu=0, noise.sigma=1)
    mat.summary["treat.all", "ATE.pos"] <- one.comparison.rules$evaluate.treat.all$ATE.test.positives
    mat.summary["treat.all", "ATE.neg"] <- one.comparison.rules$evaluate.treat.all$ATE.test.negatives
    mat.summary["treat.all", "ABR"] <- one.comparison.rules$evaluate.treat.all$ABR
    mat.summary["treat.all", "prop.pos"] <- predicting.with.oracle.treat.all$one.estimated.rule.evaluation$prop.positives
    mat.summary["treat.all", "exact.ATE.pos"] <- predicting.with.oracle.treat.all$one.estimated.rule.evaluation$ATE.positives
    mat.summary["treat.all", "exact.ATE.neg"] <- predicting.with.oracle.treat.all$one.estimated.rule.evaluation$ATE.negatives
    mat.summary["treat.all", "exact.ABR"] <- predicting.with.oracle.treat.all$one.estimated.rule.evaluation$ABR
    two.by.two.table.treat.all <- table(B.treat.all, predicting.with.oracle$one.optimal.rule.evaluation$recommended.treatment)
    one.result.two.by.two.treat.all <- SummarizeTable(two.by.two.table.treat.all)
    mat.summary["treat.all", "misclassification.rate"] <- one.result.two.by.two.treat.all$misclassification.rate
    mat.summary["treat.all", "FPR"] <- 1
    mat.summary["treat.all", "TPR"] <- 1
    return(list("training.set"=training.set,
                  "test.set"=test.set,
                   "mat.summary"=mat.summary))
}

FormatSimulation <- function(scenario,
                                               shift.for.continuous=0,
                                               n, unif.min, unif.max, noise.mu, noise.sigma,
                                               epsilon.mu, epsilon.sigma,
                                               debugging.DI=FALSE) {
    one.design <- GetDesign(n=n, unif.min=unif.min, unif.max=unif.max, noise.mu=noise.mu, noise.sigma=noise.sigma, debugging.DI=debugging.DI)
    one.response <- GetResponse(scenario=scenario,
                                                 n=n,
                                                 T=one.design$T, X=one.design$X, L=one.design$L,
                                                 shift.for.continuous=shift.for.continuous,
                                                 epsilon.mu=epsilon.mu, epsilon.sigma=epsilon.sigma)
    if (scenario %in% c("default", "main.effect.added", "quadratic.interaction.added", "piecewise.cubic.interaction.added", "noise.added",
                               "default.piecewise.constant", "two.way.interactions.piecewise.constant", "three.way.interactions.piecewise.constant", "nonlinear.interactions.piecewise.constant",
                                "high.dimensional.noise.added", "high.dimensional.noise.and.prognostic.added")) {
        one.df <- data.frame(one.design$L, one.design$X, one.design$Z,
                                     "intervention"=one.design$T,
                                     "no_relapse"=one.response$Y.binary)
    } else if (scenario %in% "chapter3") {
        one.df <- data.frame(one.design$L, one.design$X, 
                                     "intervention"=one.design$T,
                                     "no_relapse"=one.response$Y.binary)
    } else if (scenario %in% c("default.continuous", "main.effect.added.continuous", "quadratic.interaction.added.continuous", "piecewise.cubic.interaction.added.continuous", "noise.added.continuous",
                                       "default.piecewise.constant.continuous", "two.way.interactions.piecewise.constant.continuous",
                                       "three.way.interactions.piecewise.constant.continuous", "nonlinear.interactions.piecewise.constant.continuous",
                                       "high.dimensional.noise.added.continuous", "high.dimensional.noise.and.prognostic.added.continuous")) {
        one.df <- data.frame(one.design$L, one.design$X, one.design$Z,
                             "intervention"=one.design$T,
                             "time_to_relapse"=one.response$Y.continuous)
    }
    return(list("one.design"=one.design, "one.response"=one.response, "one.df"=one.df))
}

GetDesign <- function(n, p=40, q=500,
                                    unif.min, unif.max,
                                    noise.mu, noise.sigma,
                                    debugging.DI=FALSE) {
    # generate (X, Z), the patient characteristics that aren't associated with treatment or prognosis
    X <- matrix(NA, nrow=n, ncol=p)
    for (j in 1:p) {
        X[, j] <- runif(n=n, unif.min, unif.max)
    }
    colnames(X) <- paste0("input_", 1:p)
    Z <- matrix(NA, nrow=n, ncol=q)
    for (j in 1:q) {
        Z[, j] <- rnorm(n=n, mean=noise.mu, sd=noise.sigma)
    }
    colnames(Z) <- paste0("noise_", 1:q)
    # generate indicator of good prognosis
    if (debugging.DI == TRUE) {
        L <- as.matrix(rbinom(n=n, size=1, prob=1))
    } else {
        L <- as.matrix(rbinom(n=n, size=1, prob=0.5))
    }
    colnames(L) <- "good_prognosis"
    # generate treatment indicator based on prognosis
    T <- rep(NA, n)
    T[L[, 1] == 0] <- rbinom(sum(L[, 1] == 0), size=1, prob=0.75)
    T[L[, 1] == 1] <- rbinom(sum(L[, 1] == 1), size=1, prob=0.25)
    return(list("X"=X, "Z"=Z, "L"=L, "T"=T))
}

GetOptimalRule <- function(scenario, n.for.optimal=1e+05, B,
                                            use.prespecified.data, T=NULL, X=NULL, L=NULL,
                                            unif.min, unif.max, noise.mu, noise.sigma) {
    stopifnot(is.logical(use.prespecified.data))
    if (scenario %in% c("default", "main.effect.added", "quadratic.interaction.added", "piecewise.cubic.interaction.added", "noise.added",
                               "default.piecewise.constant", "two.way.interactions.piecewise.constant", "three.way.interactions.piecewise.constant", "nonlinear.interactions.piecewise.constant",
                               "high.dimensional.noise.added", "high.dimensional.noise.and.prognostic.added", "chapter3")) {
        name.object.to.return <- "prob.binary.response"
    } else if (scenario %in% c("default.continuous", "main.effect.added.continuous", "quadratic.interaction.added.continuous", "piecewise.cubic.interaction.added.continuous", "noise.added.continuous",
                                       "default.piecewise.constant.continuous", "two.way.interactions.piecewise.constant.continuous",
                                       "three.way.interactions.piecewise.constant.continuous", "nonlinear.interactions.piecewise.constant.continuous",
                                       "high.dimensional.noise.added.continuous", "high.dimensional.noise.and.prognostic.added.continuous")) {
        name.object.to.return <- "S"
    } else {
        print("invalid scenario specified")
    }
    if (use.prespecified.data == TRUE) {
        stopifnot(!is.null(T) | !is.null(X) | !is.null(L))
        stopifnot(!is.null(B))
        my.n <- nrow(X)
        one.response.all.control <- GetResponse(scenario=scenario,
                                                       n=my.n, T=rep(0, my.n),
                                                       X=X, L=L)[[name.object.to.return]]
        one.response.all.treatment <- GetResponse(scenario=scenario,
                                                       n=my.n, T=rep(1, my.n),
                                                       X=X, L=L)[[name.object.to.return]]
        one.optimal.rule.evaluation <- OneEvaluationForOptimalRule(B=NULL,
                                                                                           my.response.all.treatment=one.response.all.treatment,
                                                                                           my.response.all.control=one.response.all.control)
        one.estimated.rule.evaluation <- OneEvaluationForOptimalRule(B=B,
                                                                                              my.response.all.treatment=one.response.all.treatment,
                                                                                              my.response.all.control=one.response.all.control)
        return(list("one.optimal.rule.evaluation"=one.optimal.rule.evaluation,
                      "one.estimated.rule.evaluation"=one.estimated.rule.evaluation))
    } else {
        one.design.for.optimal <- GetDesign(n=n.for.optimal, unif.min=unif.min, unif.max=unif.max, noise.mu=noise.mu, noise.sigma=noise.sigma, debugging.DI=debugging.DI)
        one.response.all.control <- GetResponse(scenario=scenario, n=n.for.optimal,
                                                                       T=rep(0, n.for.optimal), X=one.design.for.optimal$X, L=one.design.for.optimal$L)[[name.object.to.return]]
        one.response.all.treatment <- GetResponse(scenario=scenario, n=n.for.optimal,
                                                                             T=rep(1, n.for.optimal), X=one.design.for.optimal$X, L=one.design.for.optimal$L)[[name.object.to.return]]
        one.optimal.rule.evaluation <- OneEvaluationForOptimalRule(B=NULL,
                                                                                           my.response.all.treatment=one.response.all.treatment,
                                                                                           my.response.all.control=one.response.all.control)
        return(list("one.optimal.rule.evaluation"=one.optimal.rule.evaluation))
    }
}

GetResponse <- function(scenario,
                                       n, T, X, L,
                                       shift.for.continuous=c("none", "min.and.unit.sd", "min", "10", "50"), 
                                       epsilon.mu=0, epsilon.sigma=1) {
    my.shift.for.continuous <- match.arg(shift.for.continuous)
    stopifnot(n == nrow(X))
    # compute ``score'' underlying response
    S <- rep(NA, n)
    S.mat <- matrix(nrow=n, ncol=3)
    if (scenario %in% c("default", "default.continuous", "noise.added", "noise.added.continuous")) {
        S[T == 0] <- 0 + X[T==0, 1] * 0 + L[T==0, 1] * 1.5 + X[T==0, 1] * -0.55 + 0 * 0 + 0 * 0
        S[T == 1] <- -1 + X[T==1, 1] * 0 + L[T==1, 1] * 1.5 + X[T==1, 1] * 0.55 + 0 * 0 + 0 * 0
    } else if (scenario %in% "chapter3") {
        S[T == 0] <- 0 + X[T==0, 1] * 0 + L[T==0, 1] * 1.5 + X[T==0, 1] * -0.55 + 0 * 0 + 0 * 0
        S[T == 1] <- -1.4 + X[T==1, 1] * 0 + L[T==1, 1] * 1.5 + X[T==1, 1] * 0.55 + 0 * 0 + 0 * 0
    } else if (scenario %in% c("main.effect.added", "main.effect.added.continuous")) {
        S[T == 0] <- 0 + X[T==0, 1] * 1 + L[T==0, 1] * 1.5 + X[T==0, 1] * -0.7 + 0 * 0 + 0 * 0
        S[T == 1] <- -1.4 + X[T==1, 1] * 1 + L[T==1, 1] * 1.5 + X[T==1, 1] * 0.4 + 0 * 0 + 0 * 0
    } else if (scenario %in% c("quadratic.interaction.added", "quadratic.interaction.added.continuous")) {
        S[T == 0] <- -0.2 + X[T==0, 1] * 0 + L[T==0, 1] * 1.5 + X[T==0, 1] * 1 + (X[T==0, 1]^2) * -1 + 0 * 0
        S[T == 1] <- -0.5 + X[T==1, 1] * 0 + L[T==1, 1] * 1.5 + X[T==1, 1] * -0.3 + (X[T==1, 1]^2) * 0.8 + 0 * 0
    } else if (scenario %in% c("piecewise.cubic.interaction.added", "piecewise.cubic.interaction.added.continuous")) {
        S[T == 0] <- 0 + X[T==0, 1] * 0 + L[T==0, 1] * 1.5 + X[T==0, 1] * -0.7 + (X[T==0, 1]^2) * -0.3 + 0
        S[T == 1] <- -1 + X[T==1, 1] * 0 + L[T==1, 1] * 1.5 + (X[T==1, 1] < 0.8) * (X[T==1, 1]  * 0.4  + X[T==1, 1]^3 * 2) +
                                                                                              (X[T==1, 1] > 0.8) * (X[T==1, 1]  * 2.3 + X[T==1, 1]^3 * -1)
    } else if (scenario %in% c("high.dimensional.noise.added", "high.dimensional.noise.added.continuous")) {
        S[T == 0] <- 0 + X[T==0, 1] * 0 + L[T==0, 1] * 1.5 + X[T==0, 1] * -0.7 + 0 * 0 + 0 * 0
        S[T == 1] <- -1.4 + X[T==1, 1] * 0 + L[T==1, 1] * 1.5 + X[T==1, 1] * 0.4 + 0 * 0 + 0 * 0
    } else if (scenario %in% c("high.dimensional.noise.and.prognostic.added", "high.dimensional.noise.and.prognostic.added.continuous")) {
        S[T == 0] <- 0 + X[T==0, 1] * 0 + L[T==0, 1] * 1.5 + X[T==0, 1] * -0.7 + 0 * 0 + 0 * 0 +
                            X[T==0, 2] * 0.2 + X[T==0, 3] * 0.3 + X[T==0, 4] * -0.3 + X[T==0, 5] * -0.2 + X[T==0, 6] * 0.4
        S[T == 1] <- -1.4 + X[T==1, 1] * 0 + L[T==1, 1] * 1.5 + X[T==1, 1] * 0.4 + 0 * 0 + 0 * 0 +
                            X[T==1, 2] * 0.2 + X[T==1, 3] * 0.3 + X[T==1, 4] * -0.3 + X[T==1, 5] * -0.2 + X[T==1, 6] * 0.4
    } else if (scenario %in% c("default.piecewise.constant", "default.piecewise.constant.continuous",
                                       "two.way.interactions.piecewise.constant", "two.way.interactions.piecewise.constant.continuous",
                                       "three.way.interactions.piecewise.constant", "three.way.interactions.piecewise.constant.continuous",
                                       "nonlinear.interactions.piecewise.constant", "nonlinear.interactions.piecewise.constant.continuous")) {
        S.mat[T==0, 1] <- 0.8 * (X[T==0, 1] > 0 & X[T==0, 1] < 0.7) + -1 * (X[T==0, 1] > 0.7 & X[T==0, 1] < 1.5) + 1 * (X[T==0, 1] > 1.5 & X[T==0, 1] < 2)
        S.mat[T==0, 2] <- 0.3 * (X[T==0, 2] > 0 & X[T==0, 2] < 0.5) + -0.7 * (X[T==0, 2] > 0.5 & X[T==0, 2] < 1) + 0 * (X[T==0, 2] > 1 & X[T==0, 2] < 1.3) + -0.5 * (X[T==0, 2] > 1.3 & X[T==0, 2] < 2)
        S.mat[T==0, 3] <- -0.3 * (X[T==0, 3] > 0 & X[T==0, 3] < 0.8) + 0.4 * (X[T==0, 3] > 0.8 & X[T==0, 3] < 1.3) + -0 * (X[T==0, 3] > 1.3 & X[T==0, 3] < 2)
                          
        S.mat[T==1, 1] <- 0.5 * (X[T==1, 1] > 0 & X[T==1, 1] < 0.5) + 0 * (X[T==1, 1] > 0.5 & X[T==1, 1] < 1.3) + -0.5 * (X[T==1, 1] > 1.3 & X[T==1, 1] < 2)
        S.mat[T==1, 2] <- -0.3 * (X[T==1, 2] > 0 & X[T==1, 2] < 0.6) + -1 * (X[T==1, 2] > 0.6 & X[T==1, 2] < 1.1) + 0.6 * (X[T==1, 2] > 1.1 & X[T==1, 2] < 1.4) + 0.2 * (X[T==1, 2] > 1.4 & X[T==1, 2] < 2)
        S.mat[T==1, 3] <- 0.5 * (X[T==1, 3] > 0 & X[T==1, 3] < 0.4) + -1 * (X[T==1, 3] > 0.4 & X[T==1, 3] < 0.9) + 1.1 * (X[T==1, 3] > 0.9 & X[T==1, 3] < 2)
        if (scenario %in% c("default.piecewise.constant", "default.piecewise.constant.continuous")) {
            S[T==0] <- apply(S.mat[T==0, ], 1, sum) + 0 + L[T==0, 1] * 0.5
            S[T==1] <- apply(S.mat[T==1, ], 1, sum) + 0 + L[T==1, 1] * 0.5
        } else if (scenario %in% c("two.way.interactions.piecewise.constant", "two.way.interactions.piecewise.constant.continuous")) {
            S[T==0] <- S.mat[T==0, 1] + S.mat[T==0, 2] + S.mat[T==0, 3] +
                             -0.2 * (S.mat[T==0, 1] * S.mat[T==0, 2]) +
                              0.2 * (S.mat[T==0, 1] * S.mat[T==0, 3]) +
                             -0.3 * (S.mat[T==0, 2] * S.mat[T==0, 3]) +
                             0 + L[T==0, 1] * 0.5
            S[T==1] <- S.mat[T==1, 1] + S.mat[T==1, 2] + S.mat[T==1, 3] +
                              0.1 * (S.mat[T==1, 1] * S.mat[T==1, 2]) +
                              0.3 * (S.mat[T==1, 1] * S.mat[T==1, 3]) +
                             -0.1 * (S.mat[T==1, 2] * S.mat[T==1, 3]) +
                             0 + L[T==1, 1] * 0.5
        } else if (scenario %in% c("three.way.interactions.piecewise.constant", "three.way.interactions.piecewise.constant.continuous")) {
            S[T==0] <- S.mat[T==0, 1] + S.mat[T==0, 2] + S.mat[T==0, 3] +
                             -0.5 * (S.mat[T==0, 1] * S.mat[T==0, 2] * S.mat[T==0, 3]) +
                              0 + L[T==0, 1] * 0.5
            S[T==1] <- S.mat[T==1, 1] + S.mat[T==1, 2] + S.mat[T==1, 3] +
                             0.5 * (S.mat[T==1, 1] * S.mat[T==1, 2] * S.mat[T==1, 3]) +
                              0 + L[T==1, 1] * 0.5
        } else if (scenario %in% c("nonlinear.interactions.piecewise.constant", "nonlinear.interactions.piecewise.constant.continuous")) {
            S[T==0] <- S.mat[T==0, 1] + S.mat[T==0, 2] + S.mat[T==0, 3] +
                              -0.3 * (S.mat[T==0, 1] > 1) * (S.mat[T==0, 2]) ^ 2 +
                              0.4 * (S.mat[T==0, 2] > 0.5) * (S.mat[T==0, 3]) ^ 2 +
                              -0.1 * (S.mat[T==0, 1] < 1.5) * S.mat[T==0, 3] ^3
                              0 + L[T==0, 1] * 0.5
            S[T==1] <- S.mat[T==1, 1] + S.mat[T==1, 2] + S.mat[T==1, 3] +
                              0.4 * (S.mat[T==1, 1] > 1) * (S.mat[T==1, 2]) ^ 2 +
                              -0.5 * (S.mat[T==1, 2] > 0.5) * (S.mat[T==1, 3]) ^ 2 +
                              -0.6 * (S.mat[T==1, 1] < 1.5) * S.mat[T==1, 3] ^3
                              0 + L[T==0, 1] * 0.5
        }
    }
    # convert score into probability of response (for binary outcome)
    prob.binary.response <- Expit(S)
    # get response variables
    Y.binary <- rbinom(n, size=1, prob=prob.binary.response)
    if (my.shift.for.continuous %in% c(10, 50)) {
        Y.continuous <- S + rnorm(n, mean=epsilon.mu, sd=epsilon.sigma) + as.numeric(my.shift.for.continuous)
    } else if (my.shift.for.continuous %in% c("none")) {
        Y.continuous <- S + rnorm(n, mean=epsilon.mu, sd=epsilon.sigma)
    } else if (my.shift.for.continuous %in% c("min")) {
        Y.continuous.first <- S + rnorm(n, mean=epsilon.mu, sd=epsilon.sigma)
        Y.continuous <- Y.continuous.first + abs(min(Y.continuous.first)) * 1.01
    } else if (my.shift.for.continuous %in% c("min.and.unit.sd")) {
        Y.continuous.first <- S + rnorm(n, mean=epsilon.mu, sd=epsilon.sigma)
        Y.continuous <- Y.continuous.first + abs(min(Y.continuous.first)) * 1.01
        Y.continuous <- Y.continuous / sd(Y.continuous)
    }
    return(list("S"=S, "S.mat"=S.mat,
                   "prob.binary.response"=prob.binary.response,
                   "Y.binary"=Y.binary, "Y.continuous"=Y.continuous))
}

OneEvaluationForOptimalRule <- function(my.response.all.treatment, my.response.all.control, B) {
    one.diff <- my.response.all.treatment - my.response.all.control
    if (is.null(B)) {
        one.idx.positives <- one.diff > 0
        one.idx.negatives <- one.diff <= 0
        prop.positives <- mean(one.idx.positives)
        prop.negatives <- mean(one.idx.negatives)
    } else {
        one.idx.positives <- B == 1
        one.idx.negatives <- B == 0
        prop.positives <- mean(one.idx.positives)
        prop.negatives <- mean(one.idx.negatives)
    }
    if (all(one.idx.positives == FALSE)) {
        my.ATE.positives <- NA
        my.ATE.negatives <- mean(one.diff[one.idx.negatives])
        my.ABR <- -1 * my.ATE.negatives
    } else if (all(one.idx.negatives == FALSE)) {
        my.ATE.positives <- mean(one.diff[one.idx.positives])
        my.ATE.negatives <- NA
        my.ABR <- my.ATE.positives
    } else {
        my.ATE.positives <- mean(one.diff[one.idx.positives])
        my.ATE.negatives <- mean(one.diff[one.idx.negatives])
        my.ABR <- prop.positives * mean(one.diff[one.idx.positives]) + prop.negatives * -1 * mean(one.diff[one.idx.negatives])
    }
    return(list("recommended.treatment"=as.numeric(one.diff > 0),
                "prop.positives"=prop.positives,
                "ATE.positives"=my.ATE.positives,
                "prop.negatives"=prop.negatives,
                "ATE.negatives"=my.ATE.negatives,
                "ABR"=my.ABR,
                "ABR.2"=mean(abs(my.response.all.treatment - my.response.all.control))))
}

PlotLegendMeanSummary <- function(mat.from.simulation) {
    my.n.vec <- colnames(mat.from.simulation)
    my.vec.approaches <- rownames(mat.from.simulation)
    plot(0, 0, type = "n", ann = F, axes = F)
    legend("center", cex=2.5, legend=my.vec.approaches, col=1:length(my.vec.approaches), pch=1:length(my.vec.approaches), lwd=rep(3, length(my.vec.approaches)))
}

PlotMeanSummary <- function(mat.from.simulation,
                                           my.ylab="", my.xlab="Sample Size", my.title="", my.lwd=4, my.cex.axis=1.7, my.cex.lab=1.7, my.cex.main=1.5) {
    my.n.vec <- colnames(mat.from.simulation)
    my.vec.approaches <- rownames(mat.from.simulation)
    my.ylim <- c(-0.1, max(mat.from.simulation))
    plot(my.n.vec, mat.from.simulation[1, ],
          ylim=my.ylim,
          type="b", xlab=my.xlab, ylab=my.ylab, main=my.title, 
          col=1, pch=1, lwd=my.lwd, cex.axis=my.cex.axis, cex.lab=my.cex.lab, cex.main=my.cex.main)
    for (i in 2:length(my.vec.approaches)) {
        lines(my.n.vec, mat.from.simulation[i, ], type="b", col=i, pch=i,
               lwd=my.lwd, cex.axis=my.cex.axis)
    }
}

SummarizeTable <- function(my.two.by.two.table) {
    if (nrow(my.two.by.two.table) == 1 & all(rownames(my.two.by.two.table) == "0")) {
        FPR <- 0
        TPR <- 0
        misclassification.rate <- 1 - (my.two.by.two.table["0", "0"] / sum(my.two.by.two.table))
    } else if (nrow(my.two.by.two.table) == 1 & all(rownames(my.two.by.two.table) == "1")) {
        misclassification.rate <- 1 - (my.two.by.two.table["1", "1"] / sum(my.two.by.two.table))
        FPR <- 1
        TPR <- 1
    } else {
        FPR <- my.two.by.two.table["1", "0"] / sum(my.two.by.two.table[, "0"])
        TPR <- my.two.by.two.table["1", "1"] / sum(my.two.by.two.table[, "1"])
        misclassification.rate <- 1 - ((my.two.by.two.table["0", "0"] + my.two.by.two.table["1", "1"]) / sum(my.two.by.two.table))
    }
    return(list("misclassification.rate"=misclassification.rate,
                   "TPR"=TPR,
                   "FPR"=FPR))
}
