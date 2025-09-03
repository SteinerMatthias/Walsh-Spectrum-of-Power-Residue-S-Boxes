using Oscar

"""
Takes as input a field element "a", and checks whether it
is a generator of the multiplicative group.
"""
function is_multiplicative_generator(a::FqFieldElem)
    if iszero(a)
        return false
    end
    K = parent(a)
    n = order(K) - 1
    f = factor(n)
    for p = keys(f.fac)
        if a^divexact(n, p) == one(K)
            return false
        end
    end
    return true
end

"""
Takes as input a finite field, and returns a random
multiplicative generator.
"""
function multiplicative_generator(K::FqField)
    g = rand(K)
    while !is_multiplicative_generator(g)
        g = rand(K)
    end
    return g
end

"""
Takes as input a finite field of order q, and an integer "m",
and returns the elements of the multiplicative subroup of
order (q - 1) / gcd (m, q - 1) as vector.
"""
function multiplicative_subgroup(K::FqField, m::Int64)
    q = order(K)
    @assert(mod(q - 1, m) == 0, "m must divide q - 1.")
    g = multiplicative_generator(K)
    g_m = g^m
    g_k_m = deepcopy(g_m)
    G = FqFieldElem[one(K)]
    while !isone(g_k_m)
        push!(G, g_k_m)
        g_k_m *= g_m
    end
    return G
end
