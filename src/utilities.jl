"""
    Distributions.convolve(ds::Vector{T}) where {T <: UnivariateDistribution}

Overload the `convolve` function to work on a vector of `UnivariateDistribution`.

# Arguments
- `ds::Vector{T}`: A vector of `UnivariateDistribution` objects.

# Returns
- `output`: The result of convolving all the distributions in `ds`.

"""
function Distributions.convolve(ds::Vector{T}) where {T <: UnivariateDistribution}
    try
        output = ds[1]
        for d in ds[2:end]
            output = Distributions.convolve(output, d)
        end
        return output
    catch
        if length(ds) == 0
            error("Attempted to convolve an empty vector")
        else
            error("Attempted to convolve a vector of distributions that may not have a closed-form convolution formula")
        end
    end
end