"""
Script to compute the maximal value of the Walsh spectrum of S-Boxes
"S (x) = x^d * x^((q - 1) / m)".

The script takes the following command line arguments as input:

d -- Exponent of power function.
m -- To generate subgroup of order (q - 1) / gcd (m, q - 1), where q is prime.
p_min -- Minimal value of prime interval.
p_max -- Maximal value of prime interval.
"""

using Oscar

include("Exponential_Sums.jl")

d = parse(Int64, ARGS[1])
m = parse(Int64, ARGS[2])
p_min = parse(Int64, ARGS[3])
if p_min % 2 == 0
    p_min += 1
end
p_max = parse(Int64, ARGS[4])

log_file_name = "PR_Sbox_d_" * string(d) * "_" * "m_" * string(m) * ".log"

function main()
    open(log_file_name, "w") do file
        text = "Exhastive Walsh Transform Computation for Power Residue S-Box."
        text *= "\n" * "Parameters:" * "\n"
        text *= "d: " * string(d) * "\n"
        text *= "m: " * string(m) * "\n"
        text *= "p_min: " * string(p_min) * "\n"
        text *= "p_max: " * string(p_max) * "\n"
        text *= repeat("-", 100)
        text *= "\n" * "p" *  "\t" * "(Walsh Transform) / p^(1 / 2)" * "\n"
        text *= repeat("-", 100)
        write(file, text * "\n")
        println(text)
    end

    for p in p_min:2:p_max
        if is_prime(p) && (p - 1) % m == 0
            K = GF(p)
            _, x = polynomial_ring(K)
            f = x^d * x^Int64((p - 1) / m)
            s = maximal_walsh_transform(f, constant_a=true) / sqrt(p)
            open(log_file_name, "a") do file
                text = string(p) * "\t" * string(s)
                write(file, text * "\n")
                println(text)
            end
        end
    end
end

main()
