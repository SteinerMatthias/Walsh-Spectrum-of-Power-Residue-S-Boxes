"""
Script to compute the maximal value of the Walsh spectrum for the GKRS
S-Box.

The script takes the following command line arguments as input:

d_p -- Exponent d-plus.
d_m -- Exponent d-minus.
p_min -- Minimal value of prime interval.
p_max -- Maximal value of prime interval.
"""

using Oscar

include("Exponential_Sums.jl")

d_p = parse(Int64, ARGS[1])
d_m = parse(Int64, ARGS[2])
p_min = parse(Int64, ARGS[3])
if p_min % 2 == 0
    p_min += 1
end
p_max = parse(Int64, ARGS[4])

log_file_name = "GKRS_Sbox_d_+_" * string(d_p) * "_" * "d_-_" * string(d_m) * ".log"

function main()
    open(log_file_name, "w") do file
        text = "Exhastive Walsh Transform Computation for GKRS S-Box."
        text *= "\n" * "Parameters:" * "\n"
        text *= "d_+: " * string(d_p) * "\n"
        text *= "d_-: " * string(d_m) * "\n"
        text *= "p_min: " * string(p_min) * "\n"
        text *= "p_max: " * string(p_max) * "\n"
        text *= repeat("-", 100)
        text *= "\n" * "p" *  "\t" * "(Walsh Transform) / p^(1 / 2)" * "\n"
        text *= repeat("-", 100)
        write(file, text * "\n")
        println(text)
    end

    for p in p_min:2:p_max
        if is_prime(p)
            K = GF(p)
            _, x = polynomial_ring(K)
            f = (x^d_p * (1 + x^Int64((p - 1) / 2)) + x^d_m * (1 - x^Int64((p - 1) / 2))) / 2
            s = maximal_walsh_transform(f) / sqrt(p)
            open(log_file_name, "a") do file
                text = string(p) * "\t" * string(s)
                write(file, text * "\n")
                println(text)
            end
        end
    end
end

main()
