# TODO: Currently assumes a single outcome!!
# TODO: Assumes intervention produces a table
# TODO: Needs better input processing

check_response(scm) = (length(scm.response) != 1 && throw(ArgumentError("More than one response not allowed")))
check_treatment(scm) = (length(scm.treatment) != 1 && throw(ArgumentError("More than one treatment not allowed")))


function draw_counterfactual(parents, intervention)
    intervened = intervention(parents) 

    if DataAPI.ncol(intervened) != DataAPI.ncol(treatment(parents))
        throw(ArgumentError("Table of intervened treatments does not contain the same number of treatment vectors as are specified in the StructuralCausalModel"))
    end

    unintervened = CausalTables.reject(parents, keys(intervened))
    counterfactual_covariates = CausalTables.replace(unintervened, data = merge(unintervened.data, intervened))
    Ystar = rand.(condensity(scm, counterfactual_covariates, scm.response[1]))
    return(Ystar)
end


function cfmean(scm::StructuralCausalModel, intervention::Function; samples = 10^6)
    check_input(scm)
    ct = rand(scm, samples)
    parents = responseparents(ct)
    Y_counterfactual = draw_counterfactual(parents, intervention)
    return((μ = mean(Y_counterfactual), eff_bound = var(Y_counterfactual)))
end

function cfdiff(scm::StructuralCausalModel, intervention1::Function, intervention2::Function; samples = 10^6)
    check_input(scm)
    ct = rand(scm, samples)
    parents = responseparents(ct)
    Y_cf1 = draw_counterfactual(parents, intervention1)
    Y_cf2 = draw_counterfactual(parents, intervention2)
    diff_cf = Y_cf1 .- Y_cf2
    return((μ = mean(diff_cf), eff_bound = var(diff_cf)))
end

# TODO: Assume only a single treatment!!!
treat_all(ct) = NamedTuple{Tuple(ct.treatment)}((ones(DataAPI.nrow(ct)),))
treat_none(ct) = NamedTuple{Tuple(ct.treatment)}((zeros(DataAPI.nrow(ct)),))

function ate(scm::StructuralCausalModel; samples = 10^6)
    return(cfdiff(scm, treat_all, treat_none; samples = samples))
end

function att(scm::StructuralCausalModel; samples = 10^6)
    check_response(scm)
    check_treatment(scm)

    ct = rand(scm, samples)
    ct_treated = Tables.subset(ct, Tables.getcolumn(treatment(ct), 1))
    parents = responseparents(ct_treated)
    Y_cf1 = draw_counterfactual(parents, treat_all)
    Y_cf2 = draw_counterfactual(parents, treat_none)
    diff_cf = Y_cf1 .- Y_cf2
    return((μ = mean(diff_cf), eff_bound = var(diff_cf)))
end

function atu(scm::StructuralCausalModel; samples = 10^6)
    check_response(scm)
    check_treatment(scm)

    ct = rand(scm, samples)
    ct_treated = Tables.subset(ct, .!(Tables.getcolumn(treatment(ct), 1)))
    parents = responseparents(ct_treated)
    Y_cf1 = draw_counterfactual(parents, treat_all)
    Y_cf2 = draw_counterfactual(parents, treat_none)
    diff_cf = Y_cf1 .- Y_cf2
    return((μ = mean(diff_cf), eff_bound = var(diff_cf)))
end


treatment_identity(ct) = columntable(treatment(ct))


function cast_matrix_to_table_function(func::Function) 
    ct -> Tables.columntable(Tables.table(func(Tables.matrix(treatment(ct))); header = ct.treatment))
end

function additive_mtp(δ)
    mat_func = (x -> x .+ δ)
    cast_matrix_to_table_function(mat_func)
end

function multiplicative_mtp(δ)
    mat_func = (x -> x .* δ)
    cast_matrix_to_table_function(mat_func)
end

function ape(scm::StructuralCausalModel, intervention; samples = 10^6)
    return(cfdiff(scm, intervention, treatment_identity))
end




