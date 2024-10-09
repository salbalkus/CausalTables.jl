# TODO: Currently assumes a single outcome!!
# TODO: Assumes intervention produces a table
# TODO: Needs better input processing

check_response(scm) = (length(scm.response) != 1 && throw(ArgumentError("More than one response not allowed")))
check_treatment(scm) = (length(scm.treatment) != 1 && throw(ArgumentError("More than one treatment not allowed")))


"""
    draw_counterfactual(scm::StructuralCausalModel, parents::CausalTable, intervention::Function) -> Vector

Generate counterfactual responses based on a given structural causal model (SCM), a table of response parents, and an intervention function. That is, sample the responses that would have occurred had some intervention been applied to the treatment specified by the structural causal model.

# Arguments
- `scm::StructuralCausalModel`: The structural causal model used to generate counterfactual outcomes.
- `parents::CausalTable`: A table containing the variables causally preceding the response variable.
- `intervention::Function`: A function that defines the intervention to be applied to the parent variables.

# Returns
A vector of counterfactual responses.

"""
function draw_counterfactual(scm::StructuralCausalModel, parents::CausalTable, intervention::Function)
    intervened = intervention(parents) 

    if DataAPI.ncol(intervened) != DataAPI.ncol(treatment(parents))
        throw(ArgumentError("Table of intervened treatments does not contain the same number of treatment vectors as are specified in the StructuralCausalModel"))
    end

    unintervened = CausalTables.reject(parents, keys(intervened))
    counterfactual_covariates = CausalTables.replace(unintervened, data = merge(unintervened.data, intervened))
    Ystar = rand.(condensity(scm, counterfactual_covariates, scm.response[1]))
    return(Ystar)
end


"""
    cfmean(scm::StructuralCausalModel, intervention::Function; samples = 10^6)

Approximate the counterfactual mean of the response had `intervention` been applied to the treatment, along with its efficiency bound, for a given structural causal model (SCM). These statistical quantities are approximated using Monte Carlo sampling. 

# Arguments
- `scm::StructuralCausalModel`: The SCM from which data is to be simulated.
- `intervention::Function`: The intervention function to apply to the SCM.
- `samples`: The number of samples to draw from `scm` for Monte Carlo approximation (default is 10^6). This controls the precision of the approximation.

# Returns
A named tuple containing:
- `μ`: The mean of the counterfactual outcomes.
- `eff_bound`: The variance of the counterfactual response, which is equal to the efficiency bound for IID data. If observations are correlated, this may not have a meaningful interpretation.

# Example
```@example
using Distributions
dgp = @dgp(
    L ~ Beta(2, 4),
    A ~ @.(Bernoulli(L)),
    Y ~ @.(Normal(A + L))
)
scm = StructuralCausalModel(dgp, [:A], [:Y], [:L])
cfmean(scm, treat_all)
cfmean(scm, treat_none)
```
"""
function cfmean(scm::StructuralCausalModel, intervention::Function; samples = 10^6)
    check_response(scm)
    ct = rand(scm, samples)
    parents = responseparents(ct)
    Y_counterfactual = draw_counterfactual(scm, parents, intervention)
    return((μ = mean(Y_counterfactual), eff_bound = var(Y_counterfactual)))
end

"""
    cfdiff(scm::StructuralCausalModel, intervention1::Function, intervention2::Function; samples = 10^6)

Approximate the difference between two counterfactual response means -- that under `intervention1` having been applied to the treatment, and that under `intervention2` -- for a given structural causal model (SCM), along with its efficiency bound. These statistical quantities are approximated using Monte Carlo sampling. 

# Arguments
- `scm::StructuralCausalModel`: The SCM from which data is to be simulated.
- `intervention1::Function`: The first intervention function to be contrasted.
- `intervention2::Function`: The second intervention function to be contrasted.
- `samples`: The number of samples to draw from `scm` for Monte Carlo approximation (default is 10^6). This controls the precision of the approximation.

# Returns
A named tuple containing:
- `μ`: The mean difference in counterfactual outcomes.
- `eff_bound`: The variance of the difference in counterfactual responses, which is equal to the efficiency bound for IID data. If observations are correlated, this may not have a meaningful interpretation.

# Example
```@example
using Distributions
dgp = @dgp(
    L ~ Beta(2, 4),
    A ~ @.(Bernoulli(L)),
    Y ~ @.(Normal(A + L))
)
scm = StructuralCausalModel(dgp, [:A], [:Y], [:L])
cfdiff(scm, treat_all, treat_none)
```
"""
function cfdiff(scm::StructuralCausalModel, intervention1::Function, intervention2::Function; samples = 10^6)
    check_response(scm)
    ct = rand(scm, samples)
    parents = responseparents(ct)
    Y_cf1 = draw_counterfactual(scm, parents, intervention1)
    Y_cf2 = draw_counterfactual(scm, parents, intervention2)
    diff_cf = Y_cf1 .- Y_cf2
    return((μ = mean(diff_cf), eff_bound = var(diff_cf)))
end

# TODO: Assume only a single treatment!!!
treat_all(ct) = NamedTuple{Tuple(ct.treatment)}((ones(DataAPI.nrow(ct)),))
treat_none(ct) = NamedTuple{Tuple(ct.treatment)}((zeros(DataAPI.nrow(ct)),))

"""
    ate(scm::StructuralCausalModel; samples = 10^6)

Approximate the average treatment effect (ATE) for a given structural causal model (SCM), along with its efficiency bound. This statistical quantity is approximated using Monte Carlo sampling.

# Arguments
- `scm::StructuralCausalModel`: The SCM from which data is to be simulated.
- `samples`: The number of samples to draw from `scm` for Monte Carlo approximation (default is 10^6). This controls the precision of the approximation.

# Returns
A named tuple containing:
- `μ`: The ATE approximation.
- `eff_bound`: The variance of the counterfactual response, which is equal to the efficiency bound for IID data. If observations are correlated, this may not have a meaningful interpretation.

# Example
```@example
using Distributions
dgp = @dgp(
    L ~ Beta(2, 4),
    A ~ @.(Bernoulli(L)),
    Y ~ @.(Normal(A + L))
)
scm = StructuralCausalModel(dgp, [:A], [:Y], [:L])
ate(scm, treat_all, treat_none)
```
"""
function ate(scm::StructuralCausalModel; samples = 10^6)
    return(cfdiff(scm, treat_all, treat_none; samples = samples))
end

"""
    att(scm::StructuralCausalModel; samples = 10^6)

Approximate the average treatment effect among the treated (ATT) for a given structural causal model (SCM), along with its efficiency bound. This statistical quantity is approximated using Monte Carlo sampling.

# Arguments
- `scm::StructuralCausalModel`: The SCM from which data is to be simulated.
- `samples`: The number of samples to draw from `scm` for Monte Carlo approximation (default is 10^6). This controls the precision of the approximation.

# Returns
A named tuple containing:
- `μ`: The ATT approximation.
- `eff_bound`: The variance of the counterfactual response, which is equal to the efficiency bound for IID data. If observations are correlated, this may not have a meaningful interpretation.

# Example
```@example
using Distributions
dgp = @dgp(
    L ~ Beta(2, 4),
    A ~ @.(Bernoulli(L)),
    Y ~ @.(Normal(A + L))
)
scm = StructuralCausalModel(dgp, [:A], [:Y], [:L])
att(scm, treat_all, treat_none)
```
"""
function att(scm::StructuralCausalModel; samples = 10^6)
    check_response(scm)
    check_treatment(scm)

    ct = rand(scm, samples)
    ct_treated = Tables.subset(ct, Tables.getcolumn(treatment(ct), 1))
    parents = responseparents(ct_treated)
    Y_cf1 = draw_counterfactual(scm, parents, treat_all)
    Y_cf2 = draw_counterfactual(scm, parents, treat_none)
    diff_cf = Y_cf1 .- Y_cf2
    return((μ = mean(diff_cf), eff_bound = var(diff_cf)))
end

"""
    atu(scm::StructuralCausalModel; samples = 10^6)

Approximate the average treatment effect among the untreated (ATU) for a given structural causal model (SCM), along with its efficiency bound. This statistical quantity is approximated using Monte Carlo sampling.

# Arguments
- `scm::StructuralCausalModel`: The SCM from which data is to be simulated.
- `samples`: The number of samples to draw from `scm` for Monte Carlo approximation (default is 10^6). This controls the precision of the approximation.

# Returns
A named tuple containing:
- `μ`: The ATU approximation.
- `eff_bound`: The variance of the counterfactual response, which is equal to the efficiency bound for IID data. If observations are correlated, this may not have a meaningful interpretation.

# Example
```@example
using Distributions
dgp = @dgp(
    L ~ Beta(2, 4),
    A ~ @.(Bernoulli(L)),
    Y ~ @.(Normal(A + L))
)
scm = StructuralCausalModel(dgp, [:A], [:Y], [:L])
atu(scm, treat_all, treat_none)
```
"""
function atu(scm::StructuralCausalModel; samples = 10^6)
    check_response(scm)
    check_treatment(scm)

    ct = rand(scm, samples)
    ct_treated = Tables.subset(ct, .!(Tables.getcolumn(treatment(ct), 1)))
    parents = responseparents(ct_treated)
    Y_cf1 = draw_counterfactual(scm, parents, treat_all)
    Y_cf2 = draw_counterfactual(scm, parents, treat_none)
    diff_cf = Y_cf1 .- Y_cf2
    return((μ = mean(diff_cf), eff_bound = var(diff_cf)))
end


treatment_identity(ct) = columntable(treatment(ct))

"""
    cast_matrix_to_table_function(func::Function)

Wraps a given function `func` that operates on a matrix and returns a new function that operates on a `CausalTable` object. The returned function converts the `CausalTable`'s treatment matrix to a table, applies `func` to this matrix, and then converts the result back to a column table with the same header as the original treatment matrix.

# Arguments
- `func::Function`: A function that takes a matrix as input and returns a matrix.

# Returns
- A function that takes a `CausalTable` object as input and returns a column table.

# Example
```@example
custom_intervention = cast_matrix_to_table_function(x -> exp.(x))
```
"""
function cast_matrix_to_table_function(func::Function) 
    ct -> Tables.columntable(Tables.table(func(Tables.matrix(treatment(ct))); header = ct.treatment))
end

"""
    additive_mtp(δ)

Constructs a function that adds a constant (or constant vector) δ to the treatment variable(s) in a `CausalTable` object. This function is intended to be used as an argument to `ape`.

# Arguments
- δ: The "additive shift" to be applied to the treatment variable of a `CausalTable`.

# Returns
- A function that takes a `CausalTable` object as input and returns a column table of treatments that have been shifted by δ units.

# Example
```@example
using Distributions
dgp = @dgp(
    L ~ Beta(2, 4),
    A ~ @.(Normal(L)),
    Y ~ @.(Normal(A + 2 * L + 1))
)
scm = StructuralCausalModel(dgp, [:A], [:Y], [:L])
ape(scm, additive_mtp(0.5))
```
"""
function additive_mtp(δ)
    mat_func = (x -> x .+ δ)
    cast_matrix_to_table_function(mat_func)
end

"""
    multiplicative_mtp(δ)

Constructs a function that scales the treatment variable(s) in a `CausalTable` object by a constant δ. This function is intended to be used as an argument to `ape`.

# Arguments
- δ: The "multiplicative shift" to be applied to the treatment variable of a `CausalTable`.

# Returns
- A function that takes a `CausalTable` object as input and returns a column table of treatments that have been scaled by δ units.

# Example
```@example
using Distributions
dgp = CausalTables.@dgp(
    L ~ Beta(2, 4),
    A ~ @.(Normal(L)),
    Y ~ @.(Normal(A + 2 * L + 1))
)
scm = CausalTables.StructuralCausalModel(dgp, [:A], [:Y], [:L])
ape(scm, multiplicative_mtp(2.0))
```
"""
function multiplicative_mtp(δ)
    mat_func = (x -> x .* δ)
    cast_matrix_to_table_function(mat_func)
end

"""
    ape(scm::StructuralCausalModel, intervention::Function; samples = 10^6)

Approximate the average policy effect for a given structural causal model (SCM), along with its efficiency bound. This is also known as the causal effect of a modified treatment policy, and is approximated using Monte Carlo sampling. Note that unless `intervention` is piecewise smooth invertible, the estimated statistical quantity may not have a causal interpretation; see [Haneuse and Rotnizky (2013)](https://pubmed.ncbi.nlm.nih.gov/23913589/).

Convenience functions for generating `intervention` functions include `additive_mtp` and `multiplicative_mtp`, which construct functions that respectively add or multiply a constant (or constant vector) to the treatment variable(s). One can also implement their own intervention function; this function must take as input a `CausalTable` object and return a NamedTuple object with each key indexing a treatment variable that has been modified according to the intervention. Also see `cast_matrix_to_table_function` for a convenience function for constructing interventions.

# Arguments
- `scm::StructuralCausalModel`: The SCM from which data is to be simulated.
- `intervention::Function`: The intervention function to apply to the SCM.
- `samples`: The number of samples to draw from `scm` for Monte Carlo approximation (default is 10^6). This controls the precision of the approximation.

# Returns
A named tuple containing:
- `μ`: The ATU approximation.
- `eff_bound`: The variance of the difference between the natural and counterfactual responses, which is equal to the efficiency bound for IID data. If observations are correlated, this may not have a meaningful interpretation.

# Example
```@example
using Distributions
dgp = CausalTables.@dgp(
    L ~ Beta(2, 4),
    A ~ @.(Normal(L)),
    Y ~ @.(Normal(A + 2 * L + 1))
)
scm = CausalTables.StructuralCausalModel(dgp, [:A], [:Y], [:L])
ape(scm, additive_mtp(0.5))
ape(scm, multiplicative_mtp(2.0))

# example of a custom intervention function
custom_intervention = cast_matrix_to_table_function(x -> exp.(x))
ape(scm, custom_intervention)

```
"""
function ape(scm::StructuralCausalModel, intervention::Function; samples = 10^6)
    return(cfdiff(scm, intervention, treatment_identity))
end




