GetPopulationCurves <- function(X.signal,
                                 prob.L1=0.5, prob.T1.given.L0=0.75, prob.T1.given.L1=0.25,
                                 intercept.control=0, slope.control=-0.55, gamma.control=1.5,
                                 intercept.treat=-1.4, slope.treat=0.55, gamma.treat=1.5,
                                 desirable.response, idx.test.positive=NULL) {
    if (!is.null(idx.test.positive)) {
        stopifnot(is.logical(idx.test.positive))
        user.rule <- TRUE
    } else {
        user.rule <- FALSE
    }
    # compute response.curves
    prob.T1 <- prob.T1.given.L0 * (1 - prob.L1) + prob.T1.given.L1 * prob.L1
    prob.L0.given.T1 <- prob.T1.given.L0 * (1 - prob.L1) / prob.T1
    prob.L0.given.T0 <- (1 - prob.T1.given.L0) * (1 - prob.L1) / (1 - prob.T1)
    group.weights <- GetTrueWeights(prob.L0.given.T0=prob.L0.given.T0,
                                                  prob.L0.given.T1=prob.L0.given.T1,
                                                  prob.T1.given.L0=prob.T1.given.L0,
                                                  prob.T1.given.L1=prob.T1.given.L1)
    ## incorrect (sample averaging)
    ord.X.signal <- order(X.signal) ## sorting like this only makes sense for univariate X, but that's all we're considering in this small simulation study
    X.signal.sorted <- X.signal[ord.X.signal]
    incorrect.prob.control <-  group.weights$incorrect.weight.control.L0 * Expit(intercept.control + slope.control * X.signal.sorted) + group.weights$incorrect.weight.control.L1 * Expit(intercept.control + slope.control * X.signal.sorted + gamma.control)
    incorrect.prob.treat <- group.weights$incorrect.weight.treat.L0 * Expit(intercept.treat + slope.treat * X.signal.sorted) + group.weights$incorrect.weight.treat.L1 * Expit(intercept.treat + slope.treat * X.signal.sorted + gamma.treat)
    ## correct (population averaging)
    correct.prob.control <-  (group.weights$correct.weight.control.L0 * Expit(intercept.control + slope.control * X.signal.sorted) + group.weights$correct.weight.control.L1 * Expit(intercept.control + slope.control * X.signal.sorted + gamma.control)) / (group.weights$correct.weight.control.L0 + group.weights$correct.weight.control.L1)
    correct.prob.treat <- (group.weights$correct.weight.treat.L0 * Expit(intercept.treat + slope.treat * X.signal.sorted) + group.weights$correct.weight.treat.L1 * Expit(intercept.treat + slope.treat * X.signal.sorted + gamma.treat)) / (group.weights$correct.weight.treat.L0 + group.weights$correct.weight.treat.L1)
    if (desirable.response == TRUE) {
        crossing.point <- X.signal.sorted[min(which(correct.prob.treat > correct.prob.control))]
        if (is.null(idx.test.positive) == TRUE) {
            idx.test.positive <- correct.prob.treat > correct.prob.control
        } else {
            idx.test.positive <- idx.test.positive[ord.X.signal]
        }
        benefit.under.rule <- rep(NA, length(X.signal))
        benefit.under.rule[idx.test.positive] <- (correct.prob.treat - correct.prob.control)[idx.test.positive]
        benefit.under.rule[!idx.test.positive] <- (correct.prob.control - correct.prob.treat)[!idx.test.positive]
    } else {
        crossing.point <- X.signal.sorted[min(which(correct.prob.control < correct.prob.treat))]
        if (is.null(idx.test.positive) == TRUE) {
            idx.test.positive <- correct.prob.treat < correct.prob.control
        }        
        benefit.under.rule <- rep(NA, length(X.signal))
        benefit.under.rule[idx.test.positive] <- (correct.prob.control - correct.prob.treat)[idx.test.positive]
        benefit.under.rule[!idx.test.positive] <- (correct.prob.treat - correct.prob.control)[!idx.test.positive]
    }
    response.probability.under.optimal.rule <- rep(NA, length(X.signal))
    response.probability.under.optimal.rule[idx.test.positive] <- correct.prob.treat[idx.test.positive]
    response.probability.under.optimal.rule[!idx.test.positive] <- correct.prob.control[!idx.test.positive]
    if (user.rule == FALSE) {
        return(list("ord.X.signal"=ord.X.signal,
                    "incorrect.prob.control"=incorrect.prob.control, "incorrect.prob.treat"=incorrect.prob.treat,
                    "correct.prob.control"=correct.prob.control, "correct.prob.treat"=correct.prob.treat, "crossing.point"=crossing.point,
                    "mean.response.probability.under.optimal.rule"=mean(response.probability.under.optimal.rule),
                    "mean.response.probability.treating.all"=mean(correct.prob.treat), "mean.response.probability.treating.none"=mean(correct.prob.control)))
    } else {
        return(list("ord.X.signal"=ord.X.signal,
                    "incorrect.prob.control"=incorrect.prob.control, "incorrect.prob.treat"=incorrect.prob.treat,
                    "correct.prob.control"=correct.prob.control, "correct.prob.treat"=correct.prob.treat, "crossing.point"=crossing.point,
                    "mean.response.probability.under.specified.rule"=mean(response.probability.under.optimal.rule),
                    "mean.response.probability.treating.all"=mean(correct.prob.treat), "mean.response.probability.treating.none"=mean(correct.prob.control)))
    }
}

GetTrueWeights <- function(prob.L0.given.T0, prob.L0.given.T1,
                                      prob.T1.given.L0, prob.T1.given.L1) {
    # compute weights
    ## incorrect weights (sample averaging)
    incorrect.weight.control.L0 <- prob.L0.given.T0
    incorrect.weight.control.L1 <- 1 - prob.L0.given.T0
    incorrect.weight.treat.L0 <- prob.L0.given.T1
    incorrect.weight.treat.L1 <- 1 - prob.L0.given.T1
    ## correct weights (population averaging)
    correct.weight.control.L0 <- prob.L0.given.T0 * (1 / (1 - prob.T1.given.L0))
    correct.weight.control.L1 <- (1 - prob.L0.given.T0) * (1 / (1 - prob.T1.given.L1))
    correct.weight.treat.L0 <- prob.L0.given.T1 * (1 / prob.T1.given.L0)
    correct.weight.treat.L1 <- (1 - prob.L0.given.T1) * (1 / prob.T1.given.L1)
    return(list("incorrect.weight.control.L0"=incorrect.weight.control.L0,
                   "incorrect.weight.control.L1"=incorrect.weight.control.L1,
                   "incorrect.weight.treat.L0"=incorrect.weight.treat.L0,
                   "incorrect.weight.treat.L1"=incorrect.weight.treat.L1,
                   "correct.weight.control.L0"=correct.weight.control.L0,
                   "correct.weight.control.L1"=correct.weight.control.L1,
                   "correct.weight.treat.L0"=correct.weight.treat.L0,
                   "correct.weight.treat.L1"=correct.weight.treat.L1))
}
