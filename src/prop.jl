function rk4_step(phi, dt, prop::HEOMPropagator; backwards=false)
  Hfwd = prop.mat
  Hrev = prop.mat'

  H = backwards ? Hrev : Hfwd

  k1 = H * phi
  k2 = H * (phi + k1 * dt / 2)
  k3 = H * (phi + k2 * dt / 2)
  k4 = H * (phi + k3 * dt)

  return phi + dt / 6 * (k1 + 2k2 + 2k3 + k4)
end
