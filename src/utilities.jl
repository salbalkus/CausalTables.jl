"""
    Distributions.convolve(ds::Vector{T}) where {T <: UnivariateDistribution}

Overload the `convolve` function to work on a vector of `UnivariateDistribution`.

# Arguments
- `ds::Vector{T}`: A vector of `UnivariateDistribution` objects.

# Returns
- `output`: The result of convolving all the distributions in `ds`. If `ds` is empty, will return `Binomial(0, 0.5)` denoting a point mass at 0.

"""
function Distributions.convolve(ds::Vector{T}) where {T <: UnivariateDistribution}
    try
        if length(ds) == 0
            return Binomial(0, 0.5)
        else
            return reduce(Distributions.convolve, ds)
        end
    catch
        throw(ArgumentError("Attempted to convolve a vector of distributions that may not have a closed-form convolution formula"))
    end
end

function Distributions.convolve(ds::Vector{T}, m) where {T <: UnivariateDistribution}
    try
        result = Vector{UnivariateDistribution}(undef, length(ds))
        result .= Binomial(0, 0.5)
        indices = findall(!iszero, m)
        for index in indices
            result[index[2]] = (result[index[2]] == Binomial(0, 0.5) ? ds[index[1]] : convolve(result[index[2]], ds[index[1]]))
        end
        return result
    catch
        throw(ArgumentError("Attempted to convolve a vector of distributions that may not have a closed-form convolution formula"))
    end
end


