function [mean_abs_Delta_rand_sq, mean_real_Delta_rand_sq] = compute_Delta_rand(psi0, dX, par)
% compute second moments of random component of complex-valued valley-splitting parameter
% and their functional derivatives

    int_fluct_abs    = zeros(length(par.n_range),1);
    int_fluct_real   = zeros(length(par.n_range),1);
   
    for in = 1 : length(par.n_range)
      [fA,dfA_dX,dfA_dpsi] = kernel_A(par.n_range(in), psi0, dX, par);
      [fB,dfB_dX,dfB_dpsi] = kernel_B(par.n_range(in), psi0, dX, par);

      int_fluct_abs(in)   = par.dz * sum(fA);
      int_fluct_real(in)  = par.dz * sum(fB);
    end     
  
   mean_abs_Delta_rand_sq  = real(sum(par.C4 .* int_fluct_abs));
   mean_real_Delta_rand_sq =      sum(par.D4 .* int_fluct_real);



end