__precompile__()

module SumOfSquares

export SOSModel

import Base.show, Base.length, Base.getindex, Base.vect, Base.isless, Base.isempty, Base.start, Base.done, Base.next, Base.convert, Base.dot
import Base: *, +, -, /, ^, ==,
    promote_rule, convert, show, isless, size, getindex,
    one, zero, transpose, isapprox, @pure, dot, copy


using MultivariatePolynomials
const MP = MultivariatePolynomials
using MultivariateMoments
using SemialgebraicSets

include("matpoly.jl")
include("sosdec.jl")

include("certificate.jl")

using PolyJuMP, JuMP
import JuMP: validmodel, addtoexpr_reorder

include("variable.jl")
include("constraint.jl")

function SOSModel(; solver=()->nothing, kwargs...)
    m = Model(; kwargs...)
    PolyJuMP.getpolydata(m).solver = solver
    setpolymodule!(m, SumOfSquares)
    m
end

export getslack
function PolyJuMP.getslack(c::SOSConstraint)
    JuMP.resultvalue(c.slack)
end

end # module
