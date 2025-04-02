function [Delta_det, D_Delta_det_dX, D_Delta_det_dpsi, Delta_det_n] = compute_Delta_det(psi0, x, par)
% compute deterministic part of complex-valued valley-splitting parameter

  int_QW = zeros(length(par.n_range), 1);
  int_QD = zeros(length(par.n_range), 1);

  if nargout > 1
    D_Delta_det_dX   = zeros(par.N, 1);
    D_Delta_det_dpsi = zeros(par.N, 1);
  end
         
  for in = 1 : length(par.n_range)
    [fQW,dfQW_dX,dfQW_dpsi] = kernel_I_det_l(par.n_range(in), psi0, x, par);
    [fQD,dfQD_dX,dfQD_dpsi] = kernel_I_det_t(par.n_range(in), psi0, x, par);

    int_QW(in) = par.dz * sum(fQW);
    int_QD(in) = par.dz * sum(fQD);

    if nargout > 1
      % compute functional derivates
        D_Delta_det_dX   = D_Delta_det_dX ...
                           + par.C2(in) * (dfQW_dX + dfQD_dX);
        D_Delta_det_dpsi = D_Delta_det_dpsi ...
                           + par.C2(in) * (dfQW_dpsi + dfQD_dpsi);

    end
  end

% print  
%  [par.n_range', par.C2, int_QW/par.units.meV, int_QD/par.units.meV]
    
  Delta_QW = sum(par.C2 .* int_QW);
  Delta_QD = sum(par.C2 .* int_QD);

  Delta_det = Delta_QW + Delta_QD;

% export details about contributions from resonances index by n
  if nargout >= 4
    integral_J_sum = int_QW + int_QD;
    Delta_det_n    = par.C2 .* integral_J_sum;
  end

end

