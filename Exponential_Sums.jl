using Oscar

include("Finite_Field_Utilities.jl")

"""
Fundamental additive character over a finite field.
Expects a field element as input.
"""
function fundamental_additive_character(x::FqFieldElem)
    K = parent(x)
    p = characteristic(K)
    return exp(2 * im * pi * BigInt(lift(ZZ, absolute_tr(x))) / BigInt(p))
end

"""
Character sum with a univariate polynomial argument.
Expects a univariate polynomial over a finite field as input.
"""
function character_sum_polynomial_argument(f::FqPolyRingElem)
    K = base_ring(f)
    s = Complex{BigFloat}(0)
    for el in K
        s += fundamental_additive_character(f(el))
    end
    return s
end

"""
Walsh transform of a univariate polynomial.
Expects as input a univariate polynomial over a finite field, and two masks as
elements of the finite field.
"""
function walsh_transform(f::FqPolyRingElem, a::FqFieldElem, b::FqFieldElem)
    P = parent(f)
    x = gen(P)
    return character_sum_polynomial_argument(a * x + b * f)
end

"""
Computes the maximal Walsh transform of a univariate polynomial.
Expects a univariate polynomial over a finite field.
Also takes an optional parameter "constant_a" to fix one of the masks to one.
Default value of "constant_a" is false, if set to "true" the result will only
be corrrect for special functions like power functions.
"""
function maximal_walsh_transform(f::FqPolyRingElem; constant_a=false)
    K = base_ring(f)
    s = BigFloat(0)
    if constant_a
        a = one(K)
        for b in K
            if !iszero(b)
                s_tmp = abs(walsh_transform(f, a, b))
                if s_tmp > s
                    s = s_tmp
                end
            end
        end
    else
        for a in K
            for b in K
                if !iszero(a * b)
                    s_tmp = abs(walsh_transform(f, a, b))
                    if s_tmp > s
                        s = s_tmp
                    end
                end
            end
        end
    end
    return s
end

"""
Computes the Kloosterman sum over a multiplicative subgroup over a finie field of order q.
Takes as input two masks "a" and "b" as field elements, and an integer "m"
to generate the subgroup of order (q - 1) / gcd (m, q - 1).
Also takes an optional parameter "generator" to call the function with an
explicit generator of the multiplicative group.
Default value is "nothing", if no generator is provided a random one is used.
"""
function kloosterman_sum_over_subgroup(a::FqFieldElem, b::FqFieldElem, m::Int64; generator=nothing)
    K = parent(a)
    if isnothing(generator) || !is_multiplicative_generator(K(generator))
        g = multiplicative_generator(K)
    else
        g = K(deepcopy(generator))
    end
    g_m = g^m
    g_k_m = deepcopy(g_m)
    s = Complex{BigFloat}(0)
    s += fundamental_additive_character(a * one(K) + b * one(K))
    while !isone(g_k_m)
        s += fundamental_additive_character(a * g_k_m + b * inv(g_k_m))
        g_k_m *= g_m
    end
    return s
end

"""
Computes the classical Kloosterman sum over the multiplicative group over a finite field.
Takes as input two masks "a" and "b" as field elements.
"""
function kloosterman_sum(a::FqFieldElem, b::FqFieldElem)
    return kloosterman_sum_over_subgroup(a, b, 1)
end

"""
Computes the maximal Kloosterman sum over a multiplicative subgroup of a finite field.
Takes as input a finite field of order q and an integer "m" to generate the subgroup
of order (q - 1) / gcd (m, q - 1).
"""
function maximal_kloosterman_sum_over_subgroup(K::FqField, m::Int64)
    s = BigFloat(0)
    g = multiplicative_generator(K)
    for r in 0:(m - 1)
        a = g^r
        for b in K
            if !iszero(b)
                s_tmp = abs(kloosterman_sum_over_subgroup(a, b, m, generator=g))
                if s_tmp > s
                    s = s_tmp
                end
            end
        end
    end
    return s
end
