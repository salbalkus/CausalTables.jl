
mutable struct Summary
    target::Symbol
    matrix::Symbol
end

summarize(o, x) = o.arrays[x.matrix] * o.data[x.target]