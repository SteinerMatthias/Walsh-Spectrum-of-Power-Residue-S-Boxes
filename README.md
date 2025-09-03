# Walsh-Spectrum-of-Power-Residue-S-Boxes
This repository contains the code for the numerical Kloosterman and Walsh spectrum experiments in the paper:

> **"A Note on the Walsh Spectrum of Power Residue S-Boxes"**<br>
> Authors: Matthias Johann Steiner<br>
> Link: https://arxiv.org/abs/2507.06808


## Requirements
To run the scripts one requires [Julia](https://julialang.org/) and the [OSCAR](https://www.oscar-system.org/) package.
In Julia, Oscar can be installed with
```Julia
julia> using Pkg; Pkg.add("Oscar")
```

## Usage
The repository contains the following scripts to compute Kloosterman/Walsh spectrum experiments:
- [Kloosterman_Subgroup.jl](./Kloosterman_Subgroup.jl):
Computes for all prime numbers $p$ in an interval $[p_{min}, p_{max}] \subset \mathbb{Z}$ the maximal value of the Kloosterman sum over a multiplicative subgroup of order $\frac{p - 1}{m}$.
Takes as input the command line parameters `m p_min p_max`, for example:
```shell
julia Kloosterman_Subgroup.jl 4 13 257
```
- [PR_Sbox_Inverse_Walsh_Transformjl](./PR_Sbox_Inverse_Walsh_Transform.jl): Computes for all prime numbers $p$ in an interval $[p_{min}, p_{max}] \subset \mathbb{Z}$ the maximal value of the Walsh spectrum of power residue S-Boxes $S (x) = x^{p - 2} \cdot x^\frac{p - 1}{m}$.
Takes as input the command line parameters `m p_min p_max`, for example:
```shell
julia PR_Sbox_Inverse_Walsh_Transform.jl 4 13 257
```
- [PR_Sbox_Walsh_Transform.jl](./PR_Sbox_Walsh_Transform.jl): Computes for all prime numbers $p$ in an interval $[p_{min}, p_{max}] \subset \mathbb{Z}$ the maximal value of the Walsh spectrum of power residue S-Boxes $S (x) = x^d \cdot x^\frac{p - 1}{m}$.
Takes as input the command line parameters `d m p_min p_max`, for example:
```shell
julia PR_Sbox_Walsh_Transform.jl 3 4 13 257
```
- [GKRS_Sbox_Walsh_Transform.jl](./GKRS_Sbox_Walsh_Transform.jl): Computes for all prime numbers $p$ in an interval $[p_{min}, p_{max}] \subset \mathbb{Z}$ the maximal value of the Walsh spectrum of power residue S-Boxes $S (x) = \frac{x^{d_+} \cdot \left(1 + x^\frac{p - 1}{2} \right) + x^{d_-} \cdot \left( 1 - x^\frac{p - 1}{2} \right)}{2}$.
Takes as input the command line parameters `d_plus d_mins p_min p_max`, for example:
```shell
julia GKRS_Sbox_Walsh_Transform.jl 3 5 13 127
```

The results of all scripts are printed to console and written in a log file.
