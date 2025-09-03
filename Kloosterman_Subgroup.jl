"""
Script to compute the maximal value of a Kloosterman over a subgroup
for all prime fields within a specified interval.

The script takes the following command line arguments as input:

m -- To generate subgroup of order (q - 1) / gcd (m, q - 1), where q is prime.
p_min -- Minimal value of prime interval.
p_max -- Maximal value of prime interval.
"""

using Oscar

include("Exponential_Sums.jl")

m = parse(Int64, ARGS[1])
p_min = parse(Int64, ARGS[2])
if p_min % 2 == 0
    p_min += 1
end
p_max = parse(Int64, ARGS[3])

log_file_name = "Kloosterman_" * "m_" * string(m) * ".log"

function main()
    open(log_file_name, "w") do file
        text = "Kloosterman Sum Over Multiplicative Subgroup."
        text *= "\n" * "Parameters:" * "\n"
        text *= "m: " * string(m) * "\n"
        text *= "p_min: " * string(p_min) * "\n"
        text *= "p_max: " * string(p_max) * "\n"
        text *= repeat("-", 100)
        text *= "\n" * "p" *  "\t" * "(Kloosterman Sum) / p^(1 / 2)" * "\n"
        text *= repeat("-", 100)
        write(file, text * "\n")
        println(text)
    end

    for p in p_min:2:p_max
        if is_prime(p) && (p - 1) % m == 0
            K = GF(p)
            s = maximal_kloosterman_sum_over_subgroup(K, m) / sqrt(p)
            open(log_file_name, "a") do file
                text = string(p) * "\t" * string(s)
                write(file, text * "\n")
                println(text)
            end
        end
    end
end

main()
