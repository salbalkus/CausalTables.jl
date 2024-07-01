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
            output = ds[1]
            for d in ds[2:end]
                output = Distributions.convolve(output, d)
            end
            return output
        end
    catch
        throw(ArgumentError("Attempted to convolve a vector of distributions that may not have a closed-form convolution formula"))
    end
end