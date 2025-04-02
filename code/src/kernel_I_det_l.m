function [f, df_dX, df_dpsi] = kernel_I_det_l(n, psi0, x, par)
% integral kernel and derivatives for
% deterministic contribution to valley-splitting (longitudinal / QW)

% exponential function
  G0  = compute_G0(par);
  arg = (n*G0(3) + 2*par.k0)*par.z;

  exp_iarg = exp(-1i*arg);

% gaussians
  gauss_x = exp( -(n*G0(1)*par.l_x/2)^2 );
  gauss_y = exp( -(n*G0(2)*par.l_y/2)^2 );
  gauss   = gauss_x * gauss_y;

% potential including profile variation  
  U_x = potential_modification(x, par);
  U   = par.U_F + par.U_QW + U_x;

% squared wave function (renormalize w.r.t. dz)
  psi    = psi0/sqrt(par.dz);
  psi_sq = psi.*psi;

% integral kernel
  exp_iarg_psi_sq = exp_iarg .* psi_sq;
  f   = gauss * U .* exp_iarg_psi_sq;

% integral kernel derivatives  
  if nargout > 1
    df_dX   = gauss * par.dEc .* exp_iarg_psi_sq;
    df_dpsi = gauss * U .* exp_iarg .* (2*psi);
  end


end