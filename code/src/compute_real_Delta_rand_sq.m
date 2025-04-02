function [real_Delta_rand_sq, dreal_Delta_rand_sq_dX, dreal_Delta_rand_sq_dpsi] = compute_real_Delta_rand_sq(psi0, dX, par)
% compute second moment of random component of complex-valued valley-splitting parameter
% and their functional derivatives
% This function computes < Delta_rand^2 > (without modulus), which is
% neglected in the Rician

    int_rand_real    = zeros(length(par.n_range),1);
    
    if nargout > 1
      dreal_Delta_rand_sq_dX   = zeros(par.N, 1);
      dreal_Delta_rand_sq_dpsi = zeros(par.N, 1);
    end


    for in = 1 : length(par.n_range)
      [f,df_dX,df_dpsi] = kernel_A(par.n_range(in), 4, psi0, dX, par);      
      int_rand_real(in) = par.dz * sum(f);

      % derivatives
      if nargout > 1
        dreal_Delta_rand_sq_dX   = dreal_Delta_rand_sq_dX ...
                                  + par.C4p(in) * df_dX;
        dreal_Delta_rand_sq_dpsi = dreal_Delta_rand_sq_dpsi ...
                                  + par.C4p(in) * df_dpsi;
      end
    end     
  
    real_Delta_rand_sq  = sum(par.C4p .* int_rand_real);

  


end