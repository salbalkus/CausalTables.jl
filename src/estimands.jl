# TODO: Currently assumes a single outcome
# TODO: Assumes intervention produces a table
# TODO: Needs better input processing

check_response(scm) = (length(scm.response) != 1 && throw(ArgumentError("More than one response not allowed")))
check_treatment(scm) = (length(scm.treatment) != 1 && throw(ArgumentError("More than one treatment not allowed")))
function check_treatment_binary(ct)
    treatcol = Tables.getcolumn(ct, ct.treatment[1])
    all((treatcol .== 1) .| (treatcol .== 0)) || throw(ArgumentError("Treatment variable must be binary (only either 0 or 1 valued)"))
end

@doc raw"""
    intervene(ct::CausalTable, intervention::Function)

Applies `intervention` to the treatment vector(s) within a CausalTable, and outputs a new CausalTable with the intervened treatment.

# Arguments
- `ct::CausalTable`: The data on which treatment should be intervened
- `intervention::Function`: A function that defines the intervention to be applied to the parent variables. Use `cast_matrix_to_table_function` to convert a function acting on a treatment vector or matrix to a function that acts on a `CausalTable`.

# Returns
A `CausalTable` containing the same data as `ct`, but with the treatment variable(s) modified accoding to `intervention`

# Example
```@example
using Distributions
dgp = @dgp(
    L ~ Beta(2, 4),
    A ~ @.(Bernoulli(L)),
    Y ~ @.(Normal(A + L))
)
scm = StructuralCausalModel(dgp, :A, :Y, [:L])
ct = rand(scm, 100)
intervene(ct, treat_all)
```
"""
function intervene(ct::CausalTable, intervention::Function)
    # apply intervention
    intervened = intervention(ct) 

    # check for inconsistencies in function output
    if DataAPI.ncol(intervened) != DataAPI.ncol(treatment(ct))
        throw(ArgumentError("Table of intervened treatments does not contain the same number of treatment vectors as are specified in the CausalTable"))
    end

    # merge intervened treatments with the rest of the data
    unintervened = CausalTables.reject(ct, keys(intervened))
    newtbl = CausalTables.replace(ct, data = merge(unintervened.data, intervened) |> Select(Tables.columnnames(ct))) 
    return(newtbl)
end


@doc raw"""
    draw_counterfactual(scm::StructuralCausalModel, parents::CausalTable, intervention::Function) -> Vector

Generate counterfactual responses based on a given structural causal model (SCM), a table of response parents, and an intervention function. That is, sample the responses that would have occurred had some intervention been applied to the treatment specified by the structural causal model.

# Arguments
- `scm::StructuralCausalModel`: The structural causal model used to generate counterfactual outcomes.
- `parents::CausalTable`: A table containing the variables causally preceding the response variable.
- `intervention::Function`: A function that defines the intervention to be applied to the parent variables. Use `cast_matrix_to_table_function` to convert a function acting on a treatment vector or matrix to a function that acts on a `CausalTable`.

# Returns
A vector of counterfactual responses.

"""
function draw_counterfactual(scm::StructuralCausalModel, parents::CausalTable, intervention::Function)
    counterfactual_covariates = intervene(parents, intervention)
    Ystar = rand.(condensity(scm, counterfactual_covariates, scm.response[1]))
    return(Ystar)
end


@doc raw"""
    cfmean(scm::StructuralCausalModel, intervention::Function; samples = 10^6)

Approximate the counterfactual mean of the response had `intervention` been applied to the treatment, along with its efficiency bound, for a given structural causal model (SCM). Mathematically, this estimand is

```math
E(Y(d(a)))
```

where ``d(a)`` represents an intervention on the treatment variable(s) ``A``. This statistical quantity is approximated using Monte Carlo sampling. 

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
scm = StructuralCausalModel(dgp, :A, :Y, [:L])
cfmean(scm, treat_all)
cfmean(scm, treat_none)
```
"""
function cfmean(scm::StructuralCausalModel, intervention::Function; samples = 10^6)
    check_response(scm)
    ct = rand(scm, samples)
    parents = responseparents(ct)
    Y_counterfactual = draw_counterfactual(scm, parents, intervention)
    return((μ = mean(Y_counterfactual),))
end

@doc raw"""
    cfdiff(scm::StructuralCausalModel, intervention1::Function, intervention2::Function; samples = 10^6)

Approximate the difference between two counterfactual response means -- that under `intervention1` having been applied to the treatment, and that under `intervention2` -- for a given structural causal model (SCM), along with its efficiency bound. Mathematically, this is

```math
E(Y(d_1(a)) - Y(d_2(a)))
```

where ``d_1`` and ``d_2`` represent `intervention1` and `intervention2` being applied on the treatment variable(s) ``A``. This statistical quantity is approximated using Monte Carlo sampling. 

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
scm = StructuralCausalModel(dgp, :A, :Y, [:L])
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
    return((μ = mean(diff_cf),))
end

# TODO: Assume only a single treatment!!!


"""
    treat_all(ct::CausalTable)

Intervenes on a `CausalTable` object by setting all treatment variables to 1.

# Arguments
- `ct::CausalTable`: A `CausalTable` object with a univariate binary treatment.

# Returns
A `NamedTuple` object with the same header as the treatment matrix in `ct`, where each treatment variable is set to 1.

# Example
```@example
using Distributions
dgp = @dgp(
    L ~ Beta(2, 4),
    A ~ @.(Bernoulli(L)),
    Y ~ @.(Normal(A + L))
)
scm = StructuralCausalModel(dgp, :A, :Y, [:L])
data = rand(scm, 100)
treat_all(data)
```
"""
treat_all(ct::CausalTable) = NamedTuple{Tuple(ct.treatment)}((ones(DataAPI.nrow(ct)),))

"""
    treat_all(ct::CausalTable)

Intervenes on a `CausalTable` object by setting all treatment variables to 0.

# Arguments
- `ct::CausalTable`: A `CausalTable` object with a univariate binary treatment.

# Returns
A `NamedTuple` object with the same header as the treatment matrix in `ct`, where each treatment variable is set to 0.

# Example
```@example
using Distributions
dgp = @dgp(
    L ~ Beta(2, 4),
    A ~ @.(Bernoulli(L)),
    Y ~ @.(Normal(A + L))
)
scm = StructuralCausalModel(dgp, :A, :Y, [:L])
data = rand(scm, 100)
treat_none(data)
```
"""
treat_none(ct) = NamedTuple{Tuple(ct.treatment)}((zeros(DataAPI.nrow(ct)),))

@doc raw"""
    ate(scm::StructuralCausalModel; samples = 10^6)

Approximate the average treatment effect (ATE) for a given structural causal model (SCM), along with its efficiency bound, for a univariate binary treatment. Mathematically, this is

```math
E(Y(1) - Y(0))
```
    
where ``Y(a)`` represents the counterfactual ``Y`` had the treatment ``A`` been set to ``a``. This statistical quantity is approximated using Monte Carlo sampling.

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
scm = StructuralCausalModel(dgp, :A, :Y, [:L])
ate(scm)
```
"""
function ate(scm::StructuralCausalModel; samples = 10^6)
    check_response(scm)
    check_treatment(scm)

    ct = rand(scm, samples)
    check_treatment_binary(ct)

    ct1 = intervene(ct, treat_all)
    ct0 = intervene(ct, treat_none)

    responsesymb = scm.response[1]
    diff = (conmean(scm, ct1, responsesymb) - conmean(scm, ct0, responsesymb))

    Y = vec(responsematrix(ct))
    A = vec(treatmentmatrix(ct))
    p = propensity(scm, ct, scm.treatment[1])
    μ = conmean(scm, ct, responsesymb)
    eif = (@. diff  + (((2 * A) -1) / p) * (Y - μ))

    return((μ = mean(diff), eff_bound = var(eif)))
end

@doc raw"""
    att(scm::StructuralCausalModel; samples = 10^6)

Approximate the average treatment effect among the treated (ATT) for a given structural causal model (SCM), along with its efficiency bound, for a univariate binary treatment. Mathematically, this is

```math
E(Y(1) - Y(0) \mid A = 1)
```
    
where ``Y(a)`` represents the counterfactual ``Y`` had the treatment ``A`` been set to ``a``. This statistical quantity is approximated using Monte Carlo sampling.

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
scm = StructuralCausalModel(dgp, :A, :Y, [:L])
att(scm)
```
"""
function att(scm::StructuralCausalModel; samples = 10^6)
    check_response(scm)
    check_treatment(scm)

    ct = rand(scm, samples)
    check_treatment_binary(ct)

    ct1 = intervene(ct, treat_all)
    ct0 = intervene(ct, treat_none)

    rs = scm.response[1]
    Y = vec(responsematrix(ct))
    A = vec(treatmentmatrix(ct))
    p = propensity(scm, ct1, scm.treatment[1])
    q = mean(A)

    # Compute the EIF from Theorem 1 of Kennedy, Sjolander, and Small (2015): Semiparametric causal inference in matched cohort studies.
    # Or, check against the code from npcausal: https://rdrr.io/github/ehkennedy/npcausal/src/R/att.R
    μ0 = conmean(scm, ct0, rs)

    eif = (@. (A/q) * (Y - μ0) - ((1-A)/(1-q)) * (p / (1-p)) * (Y - μ0))

    return((μ = mean(eif), eff_bound = var(eif)))
end

@doc raw"""
    atu(scm::StructuralCausalModel; samples = 10^6)

Approximate the average treatment effect among the untreated (ATU) for a given structural causal model (SCM), along with its efficiency bound, for a univariate binary treatment. Mathematically, this is

```math
E(Y(1) - Y(0) \mid A = 0)
```
    
where ``Y(a)`` represents the counterfactual ``Y`` had the treatment ``A`` been set to ``a``. This statistical quantity is approximated using Monte Carlo sampling.

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
scm = StructuralCausalModel(dgp, :A, :Y, [:L])
atu(scm)
```
"""
function atu(scm::StructuralCausalModel; samples = 10^6)
    check_response(scm)
    check_treatment(scm)

    ct = rand(scm, samples)
    check_treatment_binary(ct)

    ct1 = intervene(ct, treat_all)
    ct0 = intervene(ct, treat_none)

    rs = scm.response[1]
    Y = vec(responsematrix(ct))
    # Flip the treatment variable
    A = vec(treatmentmatrix(ct)) .== 0
    p = propensity(scm, ct0, scm.treatment[1])
    q = mean(A)

    # Invert the mean from above to estimate E{Y(1)-Y|A=0}
    μ1 = conmean(scm, ct1, rs)
    eif = (@. ((1-A) / (1-q)) * (p / (1-p)) * (Y - μ1) - (A / q) * (Y - μ1))

    return((μ = mean(eif), eff_bound = var(eif)))
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

@doc raw"""
    ape(scm::StructuralCausalModel, intervention::Function; samples = 10^6)

Approximate the average policy effect for a given structural causal model (SCM), along with its efficiency bound. This is also known as the causal effect of a modified treatment policy, and is approximated using Monte Carlo sampling. Note that unless `intervention` is piecewise smooth invertible, the estimated statistical quantity may not have a causal interpretation; see [Haneuse and Rotnizky (2013)](https://pubmed.ncbi.nlm.nih.gov/23913589/). Mathematically, this is

```math
E(Y(d(a) - Y(a))
```

where ``d(a)`` represents the intervention on the treatment variable(s) ``A``, ``Y(d(a))`` represents the counterfactual ``Y`` under treatment ``d(a)``, and ``Y(a)`` represents the counterfactual outcome under the naturally observed value of treatment. This statistical quantity is approximated using Monte Carlo sampling.

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