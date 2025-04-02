function [abs_Delta_rand_sq, dabs_Delta_rand_sq_dX, dabs_Delta_rand_sq_dpsi] = compute_abs_Delta_rand_sq(psi0, x, par)
% compute second moment of random component of complex-valued valley-splitting parameter
% and their functional derivatives

    int_rand_abs    = zeros(length(par.n_range),1);
    
    if nargout > 1
      dabs_Delta_rand_sq_dX   = zeros(par.N, 1);
      dabs_Delta_rand_sq_dpsi = zeros(par.N, 1);
    end


  % set j-parameter (multiple of j*k0 in z-integral exponent)
    j = 0;

  % compute for all n  
    for in = 1 : length(par.n_range)
      [f,df_dX,df_dpsi] = kernel_A(par.n_range(in), j, psi0, x, par);      
      int_rand_abs(in)  = par.dz * sum(f);

      % derivatives
      if nargout > 1
        dabs_Delta_rand_sq_dX   = dabs_Delta_rand_sq_dX ...
                                  + par.C4(in) * df_dX;
        dabs_Delta_rand_sq_dpsi = dabs_Delta_rand_sq_dpsi ...
                                  + par.C4(in) * df_dpsi;
      end
    end     

    abs_Delta_rand_sq  = sum(par.C4 .* int_rand_abs);

    %{
     [par.n_range, par.C4m, int_rand_abs/par.units.meV^2 ]
     abs_Delta_rand_sq/par.units.meV^2
      pause
    %}

    if abs_Delta_rand_sq < 0
      warning(['<|Delta_rand|^2> is negative: ',num2str(abs_Delta_rand_sq,'%.4e')])
      abs_Delta_rand_sq = par.sigma_failure;
      warning(['----> replace by dummy value: ',num2str(abs_Delta_rand_sq,'%.4e')])
    end


  %%%%%%%%%%%  
  % everthing should be real (filter out residual imaginary parts from inexact arithmetics)
    abs_Delta_rand_sq = real(abs_Delta_rand_sq);
    
    if nargout > 1
      dabs_Delta_rand_sq_dX   = real(dabs_Delta_rand_sq_dX);
      dabs_Delta_rand_sq_dpsi = real(dabs_Delta_rand_sq_dpsi);
    end




end