function [f, df_dX, df_dpsi] = kernel_B(n, psi0, dX, par)
% integral kernel and derivatives for
% random contribution to valley-splitting (without modulus)

% exponential function
  G0  = compute_G0(par);
  arg = (n*G0 + 4*par.k0)*par.z;

  exp_iarg = exp(1i*arg);

% prefactor
  prefactor = (par.dEc)^2 * par.Omega_a/(2*pi*par.l_x * par.l_y);

% squared wave function
  psi0_sq = psi0.*psi0;

% total profile
  X_tot = par.X_QW + dX;

% integral kernel  
  f = prefactor .* exp_iarg .* X_tot .* (1 - X_tot) .* psi0_sq .* psi0_sq /(par.dz^2);

% integral kernel derivatives  
  if nargout > 1
    df_dX   = prefactor .* exp_iarg .* (1 - 2*X_tot) .* psi0_sq .* psi0_sq /(par.dz^2);
    df_dpsi = prefactor .* exp_iarg .* X_tot .* (1 - X_tot) .* 4 .* psi0_sq .* psi0 /(par.dz^2);
  end


end