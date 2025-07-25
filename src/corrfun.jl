function bary_fit(beta, alpha, omega_c, rtol)
  K_max = 100
  J(omega) = alpha * omega / (1 + (omega / omega_c)^2)^2
  S(omega) = J(omega) * 1 / (1 - exp(-beta * omega))
  d = ComplexF64[]
  gamma = ComplexF64[]
  if alpha == 0.0
    return [0.0], [0.0], 1
  end

  npoints = Int(1e4)
  omega = collect(range(-2 * omega_c, 2 * omega_c, length=npoints))

  Ga = aaa(omega, real(S.(omega)), tol=rtol, mmax=2 * K_max + 1, verbose=false)
  poles = prz(Ga)[1]
  residuals = prz(Ga)[2]
  for i in eachindex(poles)
    if imag(poles[i]) <= 0
      push!(d, -2 * im * residuals[i])
      push!(gamma, im * poles[i])
    end
  end
  return d, gamma, length(d)
end
