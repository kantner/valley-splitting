function [f, df_dX, df_dpsi] = kernel_A(n, j, psi0, x, par)
% integral kernel and derivatives for
% random contribution to valley-splitting
% with j=0 (with Delta_rand with modulus)
%   or j=4 (with Delta_rand without modulus)

% check if j = 0 or 4
  if j~=0 && j~=4
    error('this is only allowed for j=0 or j=4')
  end

% exponential function
  G0  = compute_G0(par);
  arg = (n*G0(3) + j*par.k0)*par.z;

  exp_iarg = exp(-1i*arg);

% gaussians
  gauss_x = exp( -0.5 * (n*G0(1)*par.l_x/2)^2 );
  gauss_y = exp( -0.5 * (n*G0(2)*par.l_y/2)^2 );
  gauss   = gauss_x * gauss_y;  

% prefactor
  prefactor = (par.dEc)^2 * par.Omega_a/(2*pi*par.l_x * par.l_y);

% squared wave function (renormalize w.r.t. dz)
  psi    = psi0/sqrt(par.dz);
  psi_sq = psi.*psi;
  
% total profile
  X = par.X_QW + x;

  %{ 
  if sum(X > 1) > 0
    warning('X can not be larger than 1')
    disp('X can not be larger than 1')
    %pause
  end

  if sum(X < 0) > 0
    warning('X can not be smaller than 0')
    disp('X can not be smaller than 0')
    %pause
  end
  %}


% integral kernel 
  gauss_prefactor = gauss * prefactor;
  X_1mX           = X .* (1 - X);
  exp_iarg_psi_sq = exp_iarg .* psi_sq;
  exp_iarg_psi3   = exp_iarg_psi_sq .* psi;
  exp_iarg_psi4   = exp_iarg_psi_sq .* psi_sq;


  f = gauss_prefactor * X_1mX .* exp_iarg_psi4;
  
% integral kernel derivatives  
  if nargout > 1
    df_dX   = gauss_prefactor * (1 - 2*X) .* exp_iarg_psi4;
    df_dpsi = gauss_prefactor * 4 * X_1mX .* exp_iarg_psi3;
  end


end