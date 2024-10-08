# TODO: Currently assumes a single outcome!!
# TODO: Assumes intervention produces a table
# TODO: Needs better input processing

function apply_intervention(parents, intervention)
    intervened = intervention(parents) 

    if DataAPI.ncol(intervened) != DataAPI.ncol(treatment(parents))
        throw(ArgumentError("Table of intervened treatments does not contain the same number of treatment vectors as are specified in the StructuralCausalModel"))
    end

    unintervened = CausalTables.reject(parents, keys(intervened))
    return(CausalTables.replace(unintervened, data = merge(unintervened.data, intervened)))
end

draw_counterfactual(counterfactual_covariates) = rand.(condensity(scm, counterfactual_covariates, scm.response[1]))

function cfmean(scm::StructuralCausalModel, intervention::Function; samples = 10^6)

    if(length(scm.response) != 1)
        throw(ArgumentError("More than one response not allowed"))
    end
    ct = rand(scm, samples)
    parents = responseparents(ct)
    counterfactual_covariates = apply_intervention(parents, intervention)
    Ystar = draw_counterfactual(counterfactual_covariates)
    return((μ = mean(Ystar), eff_bound = var(Ystar)))
end

function cfdiff(scm::StructuralCausalModel, intervention1::Function, intervention2::Function; samples = 10^6)

    if(length(scm.response) != 1)
        throw(ArgumentError("More than one response not allowed"))
    end

    ct = rand(scm, samples)
    parents = responseparents(ct)
    counterfactual_covariates1 = apply_intervention(parents, intervention1)
    counterfactual_covariates2 = apply_intervention(parents, intervention2)

    Ystar1 = draw_counterfactual(counterfactual_covariates1)
    Ystar2 = draw_counterfactual(counterfactual_covariates2)
    diffstar = Ystar1 .- Ystar2

    return((μ = mean(diffstar), eff_bound = var(diffstar)))
end

# TODO: Assume only a single treatment!!!
treat_all(ct) = NamedTuple{Tuple(ct.treatment)}((ones(DataAPI.nrow(ct)),))
treat_none(ct) = NamedTuple{Tuple(ct.treatment)}((zeros(DataAPI.nrow(ct)),))

function ate(scm::StructuralCausalModel; samples = 10^6)
    return(cfdiff(scm, treat_all, treat_none))
end



