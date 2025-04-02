function [f, df_dX, df_dpsi] = kernel_I_det_t(n, psi0, x, par)
% integral kernel and derivatives for
% deterministic contribution to valley-splitting (transversal / QD)

% exponential function
  G0  = compute_G0(par);
  arg = (n*G0(3) + 2*par.k0)*par.z;

  exp_iarg = exp(-1i*arg);

% gaussians
  gauss_x = exp( -(n*G0(1)*par.l_x/2)^2 );
  gauss_y = exp( -(n*G0(2)*par.l_y/2)^2 );
  gauss   = gauss_x * gauss_y;

% prefactor
  prefactor =   0.25 * par.const.hbar * par.omega_x * (1 - (n*G0(1)*par.l_x/sqrt(2))^2 ) ...
              + 0.25 * par.const.hbar * par.omega_y * (1 - (n*G0(2)*par.l_y/sqrt(2))^2 );

% squared wave function (renormalize w.r.t. dz)
  psi    = psi0/sqrt(par.dz);
  psi_sq = psi.*psi;

% integral kernel  
  f = gauss * prefactor .* exp_iarg .* psi_sq;

% integral kernel derivatives
  if nargout > 1
    df_dX   = zeros(size(f));
    df_dpsi = gauss * prefactor .* exp_iarg .* (2 * psi);
  end


end