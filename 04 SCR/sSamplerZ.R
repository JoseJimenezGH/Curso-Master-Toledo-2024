custom_zs_sampler <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        zNode <- target
        sNodes <- control$sNodes
        birthProb <- control$birthProb
        deathProb <- control$deathProb
        calcNodes <- model$getDependencies(c(zNode, sNodes))
    },
    run = function() {
        currentIndicatorValue <- model[[zNode]]
        currentLogProb <- getLogProb(model, calcNodes)
        if(currentIndicatorValue == 0) {
            if(runif(1,0,1) > birthProb) return()
            ## propose birth
            logProbReverseProposalValues <- getLogProb(model, sNodes)
            simulate(model, sNodes)
            logProbProposalValues <- calculate(model, sNodes)
            model[[zNode]] <<- 1
            proposalLogProb <- calculate(model, calcNodes)
            logMHR <- proposalLogProb - currentLogProb + log(deathProb) + logProbReverseProposalValues - log(birthProb) - logProbProposalValues
        } else {
            if(runif(1,0,1) > deathProb) return()
            ## propose death
            logProbReverseProposalValues <- getLogProb(model, sNodes)
            simulate(model, sNodes)
            logProbProposalValues <- calculate(model, sNodes)
            model[[zNode]] <<- 0
            proposalLogProb <- calculate(model, calcNodes)
            logMHR <- proposalLogProb - currentLogProb + log(birthProb) + logProbReverseProposalValues - log(deathProb) - logProbProposalValues
        }

        jump <- decide(logMHR)
        if(jump)
            copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else
            copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list(
        reset = function () {}
    )
)



custom_s_with_indicator_sampler <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        scaleFactor <- control$scaleFactor
        scaleNode <- control$scaleNode
        calcNodes <- model$getDependencies(target)
        target1 <- target[1]
        target2 <- target[2]
        indicatorNode <- control$indicatorNode
    },
    run = function() {
        if(model[[indicatorNode]]==0) return()
        propSD <- model[[scaleNode]]*scaleFactor
        model[[target1]] <<- rnorm(1, mean = model[[target1]], sd = propSD)
        model[[target2]] <<- rnorm(1, mean = model[[target2]], sd = propSD)
        logMHR <- calculateDiff(model, calcNodes)
        jump <- decide(logMHR)
        if(jump)
            copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else
            copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list(
        reset = function () {}
    )
)

