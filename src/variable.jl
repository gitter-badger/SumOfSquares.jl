export polytype, createpoly, DSOSPoly, SDSOSPoly, SOSPoly

function JuMP.resultvalue(p::MatPolynomial{JuMP.Variable})
    MatPolynomial(map(JuMP.resultvalue, p.Q), p.x)
end

for poly in (:DSOSPoly, :SDSOSPoly, :SOSPoly)
    @eval begin
        struct $poly{MT, MV} <: PolyJuMP.AbstractPoly
            x::MV
        end
        $poly{MT}(x::MV) where {MT, MV} = $poly{MT, MV}(x)
        $poly(x::MV) where MV = $poly{:Default}(x)
    end
end

const PosPoly{MT, MV} = Union{DSOSPoly{MT, MV}, SDSOSPoly{MT, MV}, SOSPoly{MT, MV}, Poly{true, MT, MV}}

polytype(m::JuMP.Model, p) = _polytype(m, p, p.x)
polytype(m::JuMP.Model, p, X::AbstractVector) = _polytype(m, p, monovec(X))

# Free polynomial

_polytype(m::JuMP.Model, ::Poly{false}, x::AbstractVector{MT}) where MT<:AbstractMonomial = polynomialtype(MT, JuMP.Variable)

# x should be sorted and without duplicates
function _createpoly(m::JuMP.Model, ::Poly{false}, x::AbstractVector{<:AbstractMonomial}, binary::Bool, integer::Bool)
    function _newvar(i)
        v = Variable(m)
        if binary
            setbinary(v)
        end
        if integer
            setinteger(v)
        end
        v
    end
    polynomial(_newvar, x)
end
function createpoly(m::JuMP.Model, p::Union{Poly{false, :Default}, Poly{false, :Classic}}, binary::Bool, integer::Bool)
    _createpoly(m, p, monovec(p.x), binary, integer)
end
function createpoly(m::JuMP.Model, p::Poly{false, :Gram}, binary::Bool, integer::Bool)
    _createpoly(m, p, monomials(sum(p.x)^2), binary, integer)
end

# Sum-of-Squares polynomial

_polytype(m::JuMP.Model, ::PosPoly, x::MVT) where {MT<:AbstractMonomial, MVT<:AbstractVector{MT}} = MatPolynomial{JuMP.Variable, MT, MVT}

function _constraintmatpoly!(m, p, ::Union{SOSPoly, Poly{true}})
    JuMP.addconstraint(m, JuMP.SDVariableConstraint(p.Q))
end
function _constraintmatpoly!(m, p, ::DSOSPoly)
    n = length(p.x)
    Q = Matrix{JuMP.Variable}(n, n)
    for i in 1:n
        for j in 1:n
            if i == j
                Q[i, j] = p[i, j]
            else
                Q[j, i] = Q[i, j] = Variable(m)
                @constraint m Q[i, j] >= p[i, j]
                @constraint m Q[i, j] >= -p[i, j]
            end
        end
    end
    for i in 1:n
        @constraint m 2Q[i, i] >= sum(Q[i, :])
    end
    # If n > 1, this is implied by the constraint but it doesn't hurt to add the variable cone
    # Adding things on varCones makes JuMP think that it is SDP
    # push!(m.varCones, (:NonNeg, map(i -> p[i, i].col, 1:n)))
end
function _matpolynomial(m, x::AbstractVector{<:AbstractMonomial}, binary::Bool, integer::Bool)
    if isempty(x)
        zero(polytype(m, SOSPoly(x)))
    else
        function _newvar(i, j)
            v = Variable(m)
            if length(x) == 1
                # 1x1 matrix is SDP iff its only entry is nonnegative
                # We handle this case here and do not create any SDP constraint
                setlowerbound(v, 0)
            end
            if binary
                setbinary(v)
            end
            if integer
                setinteger(v)
            end
            v
        end
        MatPolynomial{JuMP.Variable}(_newvar, x)
    end
end
function _createpoly(m::JuMP.Model, set::PosPoly, x::AbstractVector{<:AbstractMonomial}, binary::Bool, integer::Bool)
    p = _matpolynomial(m, x, binary, integer)
    if length(x) > 1
        _constraintmatpoly!(m, p, set)
    end
    p
end
function createpoly(m::JuMP.Model, p::Union{PosPoly{:Default}, PosPoly{:Gram}}, binary::Bool, integer::Bool)
    _createpoly(m, p, monovec(p.x), binary, integer)
end
function createpoly(m::JuMP.Model, pp::PosPoly{:Classic}, binary::Bool, integer::Bool)
    p = _createpoly(m, pp, getmonomialsforcertificate(pp.x), binary, integer)
    # The coefficients of a monomial not in Z do not all have to be zero, only their sum
    addpolyconstraint!(m, removemonomials(Polynomial(p), p.x), ZeroPoly(), FullSpace())
    p
end
