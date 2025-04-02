function [X_ww] = compute_wiggle_well_amplitude(X_budget, q,phi,par)
% compute amplitude of wiggle well such that Ge-budget is met

% initialize
  h = par.h_QW;
  %X_ww_0 = X_budget;
  X_ww_0 = X_budget /(0.5*( 1 + (sin(q*h-phi) +sin(phi))/(q*h) ));

% function  
  mean_mode = 2;
 
  f = @(x) X_budget - mean_Ge_budget(mean_mode, x_wiggle_well(x, q, phi, par), zeros(par.N,1), par);

  opts = optimset('TolX',1E-36,'Display','iter');
  X_ww = fzero(f, X_ww_0,opts);

end
