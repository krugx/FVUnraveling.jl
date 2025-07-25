# FVUnraveling

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://krugx.github.io/FVUnraveling.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://krugx.github.io/FVUnraveling.jl/dev/)
[![Build Status](https://github.com/krugx/FVUnraveling.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/krugx/FVUnraveling.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Current use

```julia
using FVUnraveling
using SparseArrays
using Plots
## Define Constants
const GHz = 2π
const MHz = 0.001GHz
const ns = 1.0
const μs = 1000ns

function main()
  # System setup
  omega_q = 1.0GHz
  H = omega_q * [0 1; 1 0]
  q = [1 0; 0 -1]

  ## Create a time grid 
  N = Int(100)
  tf = 1 * pi / omega_q
  dt = tf / N
  tgrid = [0:dt:tf;]


  ## Reservoir parameters
  beta = 1.0 / omega_q # inverse reservoir temperature
  omega_c = 2.0 * omega_q # cutoff frequency
  alpha = 0.1 * omega_q  # coupling constant
  rtol = 1e-1

  d, gamma, K = bary_fit(beta, alpha, omega_c, rtol)

  display(d)

  T = 2
  hild = size(H, 1)

  structure = build_heom_structure(K, T, hild; useLtrunc=true)

  prop = HEOMPropagator(structure, H, q, d, gamma)

  ado = spzeros(ComplexF64, length(structure.ados) * hild^2)
  rho = ComplexF64[1 0; 0 0]
  ado[1:hild^2] = reshape(rho, (hild^2,))

  obs = zeros(N + 1)
  obs[1] = 1

  for i in 1:N
    ado = rk4_step(ado, dt, prop)
    display(ado)
    obs[i+1] = real(ado[1])
  end

  P = plot(tgrid, obs)
  display(P)

  return 0
end
main()
```
